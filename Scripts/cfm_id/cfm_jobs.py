"""Create MS2 spectra using CFM-ID."""

import os
import re
import sys
from copy import deepcopy
from subprocess import call

from minedatabase.databases import MINE

# pylint: disable=redefined-outer-name


def start_cfm_jobs(file_dir, db, spec_type, job_template=None,
                   job_comp_number=100, max_mass=1500):
    """Run CFM-predict to predict MS2 spectra for MINE compounds.

    Parameters
    ----------
    file_dir : str
        Path to output directory to contain cfm-predict spectra.
    db : MINE
        MINE DB object containing compounds to get spectra for.
    spec_type : str
        'pos', 'neg', or 'ei' spectra type.
    job_template : str, optional (default: None)
        Path to jobfile to run cfm-predict on Quest.
    job_comp_number : int, optional (default: 100)
        Number of compounds to run per job.
    max_mass : numeric, optional (default: 1500)
        Upper limit on mass for compounds (in Da) to get spectra for.
    """
    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
    os.chdir(file_dir)

    out_file = open(db.name + "_1.txt", "w")
    i = 1
    j = 1
    inchi_set = set()
    if spec_type == 'pos':
        search = {"_id": {"$regex": "^C"}, "Mass": {"$lt": max_mass},
                  "Charge": 0, "Pos_CFM_spectra": {"$exists": 0}}
    elif spec_type == 'neg':
        search = {"_id": {"$regex": "^C"}, "Mass": {"$lt": max_mass},
                  "Charge": 0, "Neg_CFM_spectra": {"$exists": 0}}
    elif spec_type == 'ei':
        search = {"_id": {"$regex": "^C"}, "Mass": {"$lt": max_mass},
                  "Charge": 0, "EI_CFM_spectra": {"$exists": 0}}
    else:
        raise ValueError('invalid spectrum spec_type')

    for compound in db.compounds.find(search, {"SMILES": 1, "Inchikey": 1}):
        con_block = compound['Inchikey'].split('-')[0]
        if con_block not in inchi_set:
            out_file.write("%s %s\n" % (con_block, compound['SMILES']))
            i += 1
            inchi_set.add(con_block)
            if not i % job_comp_number:
                out_file.close()
                j += 1
                out_file = open(db_name + "_" + str(j) + ".txt", 'w')
    out_file.close()

    if job_template:
        # for Quest?
        # fetch the information from the job template
        try:
            infile = open(job_template)
        except IOError:
            sys.exit("Job template file not found")
        temp_data = infile.readlines()
        infile.close()

        for comp_file in os.listdir("."):
            if ".txt" not in comp_file:
                continue
            job_data = deepcopy(temp_data)
            # edit the job template data for the specific job
            job_data[1] = "#$ -e logs/cfm_%s.err\n" % (comp_file.strip('.txt'))
            job_data[2] = "#$ -o logs/cfm_%s.log\n" % (comp_file.strip('.txt'))
            new_file = file_dir + "/" + comp_file.strip()
            job_data[11] = job_data[11].replace("{smile_file}", new_file)

            # write the specific job file
            job_file = "cfm_%s.job" % comp_file.strip('.txt')
            outfile = open(job_file, 'w')
            outfile.writelines(job_data)
            outfile.close()
            rc = call(["qsub", job_file], shell=True)
            if not rc:
                os.remove(job_file)


def load_cfm_results(result_dir, db_list, spec_type, cleanup=False):
    """Insert CFM-predict results into MINE DB.

    Parameters
    ----------
    result_dir : str
        Path to directory containing spectra predictions from cfm-predict.
    db_list : list
        List of MINE DBs to update.
    spec_type : str
        Mode for spectra prediction ('pos', 'neg' or 'ei').
    cleanup : bool, optional (default: False)
        Whether to delete all spectra (in result_dir) after loading them into
        the MongoDB.
    """
    for spec_file in os.listdir(result_dir):
        if ('param' in spec_file) or (spec_file[-4:] != ".log"):
            continue
        with open(os.path.join(result_dir, spec_file)) as infile:
            data = []
            for x in re.split(r"energy\d\n", infile.read())[1:]:
                split_data = []
                for y in x.strip().split('\n'):
                    split_data.append([float(z) for z in y.split()])
                data.append(split_data)
            try:
                for db in db_list:
                    if spec_type == 'pos':
                        db.compounds.update_many(
                            {"Inchikey": {"$regex": "^" + spec_file[:-4]}},
                            {"$set": {"Pos_CFM_spectra": {"10 V": data[0],
                                                          "20 V": data[1],
                                                          "40 V": data[2]}}})
                    elif spec_type == 'neg':
                        db.compounds.update_many(
                            {"Inchikey": {"$regex": "^" + spec_file[:-4]}},
                            {"$set": {"Neg_CFM_spectra": {"10 V": data[0],
                                                          "20 V": data[1],
                                                          "40 V": data[2]}}})
                    elif spec_type == 'ei':
                        db.compounds.update_many(
                            {"Inchikey": {"$regex": "^" + spec_file[:-4]}},
                            {"$set": {"EI_CFM_spectra": {"70 V": data[0]}}})
                    else:
                        raise ValueError('invalid spectrum spec_type')
            except IndexError:
                print(spec_file)
        if cleanup:
            os.remove(os.path.join(result_dir, spec_file))


if __name__ == "__main__":

    # pylint: disable=invalid-name
    # collect user input
    if sys.argv[1] == "calculate":
        db_name = sys.argv[2]
        file_dir = sys.argv[3]
        job_comp_number = int(sys.argv[4])
        spec_type = sys.argv[5]
        if len(sys.argv) == 7:
            job_template = sys.argv[6]
        else:
            job_template = None

        db = MINE(db_name)
        start_cfm_jobs(file_dir, db, spec_type, job_template=job_template,
                       job_comp_number=job_comp_number)

    if sys.argv[1] == 'load':
        result_dir = sys.argv[2]
        spec_type = sys.argv[3]
        dbs = [MINE(x) for x in sys.argv[4:]]
        load_cfm_results(result_dir, dbs, spec_type=spec_type)
