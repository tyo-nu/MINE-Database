"""

"""
import sys
import os
from copy import deepcopy
from subprocess import call
import re
from minedatabase.databases import MINE

def start_cfm_jobs(type=True, max_mass=1500):
    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
    os.chdir(file_dir)

    out_file = open(db_name+"_1.txt", "w")
    i = 1
    j = 1
    inchi_set = set()
    if type == 'pos':
        search = {"_id": {"$regex": "^C"}, "Mass": {"$lt": max_mass},
                  "Charge": 0, "Pos_CFM_spectra": {"$exists": 0}}
    elif type == 'neg':
        search = {"_id": {"$regex": "^C"}, "Mass": {"$lt": max_mass},
                  "Charge": 0, "Neg_CFM_spectra": {"$exists": 0}}
    elif type == 'ei':
        search = {"_id": {"$regex": "^C"}, "Mass": {"$lt": max_mass},
                  "Charge": 0,  "EI_CFM_spectra": {"$exists": 0}}
    else:
        raise ValueError('invalid spectrum type')

    for compound in db.compounds.find(search, {"SMILES": 1, "Inchikey": 1}):
        con_block = compound['Inchikey'].split('-')[0]
        if con_block not in inchi_set:
            out_file.write("%s %s\n" % (con_block, compound['SMILES']))
            i += 1
            inchi_set.add(con_block)
            if not (i % job_comp_number):
                out_file.close()
                j += 1
                out_file = open(db_name+"_"+str(j)+".txt", 'w')
    out_file.close()

    if job_template:
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
            #edit the job template data for the specific job
            job_data[1] = "#$ -e logs/cfm_%s.err\n" % (comp_file.strip('.txt'))
            job_data[2] = "#$ -o logs/cfm_%s.log\n" % (comp_file.strip('.txt'))
            job_data[11] = job_data[11].replace("{smile_file}", file_dir + "/" + comp_file.strip())

            #write the specific job file
            job_file = "cfm_%s.job" % comp_file.strip('.txt')
            outfile = open(job_file, 'w')
            outfile.writelines(job_data)
            outfile.close()
            rc = call(["qsub", job_file], shell=True)
            if not rc:
                os.remove(job_file)


def load_cfm_results(result_dir, db_list, type=True, cleanup=False):
    for spec_file in os.listdir(result_dir):
        if ('param' in spec_file) or (spec_file[-4:] != ".log"):
            continue
        with open(os.path.join(result_dir, spec_file)) as infile:
            data = []
            for x in re.split("energy\d\n", infile.read())[1:]:
                split_data = []
                for y in x.strip().split('\n'):
                    split_data.append([float(z) for z in y.split()])
                data.append(split_data)
            try:
                for db in db_list:
                    if type == 'pos':
                        db.compounds.update_many(
                            {"Inchikey": {"$regex": "^" + spec_file[:-4]}},
                            {"$set": {"Pos_CFM_spectra": {"10 V": data[0],
                                                          "20 V": data[1],
                                                          "40 V": data[2]}}})
                    elif type == 'neg':
                        db.compounds.update_many(
                            {"Inchikey": {"$regex": "^" + spec_file[:-4]}},
                            {"$set": {"Neg_CFM_spectra": {"10 V": data[0],
                                                          "20 V": data[1],
                                                          "40 V": data[2]}}})
                    elif type == 'ei':
                        db.compounds.update_many(
                            {"Inchikey": {"$regex": "^" + spec_file[:-4]}},
                            {"$set": {"EI_CFM_spectra": {"70 V": data[0]}}})
                    else:
                        raise ValueError('invalid spectrum type')
            except IndexError:
                print(spec_file)
        if cleanup:
            os.remove(os.path.join(result_dir, spec_file))


if __name__ == "__main__":
    #collect user input
    if sys.argv[1] == "calculate":
        db_name = sys.argv[2]
        file_dir = sys.argv[3]
        db = MINE(db_name)
        job_comp_number = int(sys.argv[4])
        if len(sys.argv) == 7:
            job_template = sys.argv[6]
        else:
            job_template = False
        start_cfm_jobs(type=sys.argv[5])

    if sys.argv[1] == 'load':
        dbs = [MINE(x) for x in sys.argv[4:]]
        load_cfm_results(sys.argv[2], dbs, type=sys.argv[3], cleanup=True)
