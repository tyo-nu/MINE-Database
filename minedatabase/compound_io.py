import csv
import datetime
import os
import sys

from . import utils
from .databases import MINE
from rdkit.Chem import AllChem

"""Databases.py: Functions to load MINE databases from and dump compounds into common cheminformatics formats"""


def export_sdf(mine_db, dir_path, max_compounds=None):
    """
    Exports compounds from the database as an MDL SDF file
    :param mine_db: a MINE object
    :param path: directory for files
    :param max_compounds: maximum number of compounds per file (defaults to 100000)
    :return:
    """
    if not mine_db.compounds.find_one({"Product_of": {'$exists': 1}}):
        mine_db.add_rxn_pointers()

    print("Exporting %s compounds from %s as an SDF file" %(mine_db.compounds.count(), mine_db.name))
    target = utils.prevent_overwrite(os.path.join(dir_path, mine_db.name) + "_1.sdf")
    w = AllChem.SDWriter(target)
    w.SetKekulize(True)
    n_files = 1
    for compound in mine_db.compounds.find():
        mol = AllChem.MolFromSmiles(compound['SMILES'], True, {'CoA': '*', 'R': "*"})
        if mol:
            mol.SetProp('_id', compound['_id'])
            mol.SetProp('Generation', str(compound['Generation']))
            if 'Reactant_in' in compound:
                mol.SetProp('Reactant_in', str(compound['Reactant_in']))
            if 'Product_of' in compound:
                mol.SetProp('Product_of', str(compound['Product_of']))
            w.write(mol)
            if max_compounds and (w.NumMols() >= max_compounds):
                n_files += 1
                target = utils.prevent_overwrite(os.path.join(dir_path, mine_db.name) + "_%s.sdf" % n_files)
                w = AllChem.SmilesWriter(target)
    w.close()


def export_smiles(mine_db, dir_path, max_compounds=None):
    """
    Exports compounds from the database as an SMILES file
    :param mine_db: a MINE object
    :param dir_path: directory for files
    :param max_compounds: maximum number of compounds per file (defaults to unlimited)
    :return:
    """
    header = ['SMILES', "_id", "Generation", 'Reactant_in', 'Product_of']
    if not mine_db.compounds.find_one({"Product_of": {'$exists': 1}}):
        mine_db.add_rxn_pointers()

    print("Exporting %s compounds from %s as SMILES file" % (mine_db.compounds.count(), mine_db.name))
    target = open(utils.prevent_overwrite(os.path.join(dir_path, mine_db.name) + "_1.smiles"), 'w')

    w = csv.DictWriter(target, fieldnames=header, dialect='excel-tab')
    n_files = 1
    i = 0
    for compound in mine_db.compounds.find({}, dict([(x, 1) for x in header])):
        w.writerow(compound)
        i += 1
        if max_compounds and not (i % max_compounds):
            n_files += 1
            target = open(utils.prevent_overwrite(os.path.join(dir_path, mine_db.name) + "_%s.smiles" % n_files), 'w')
            w = csv.DictWriter(target, fieldnames=header, dialect='excel-tab')


def export_mol(mine_db, target, name_field='_id'):
    """
    Exports compounds from the database as MDL molfiles
    :param mine_db: The database to export
    :type mine_db: a MINE object
    :param target: a directory in which to place the files
    :type target: string
    :param name_field: the field to provide names for the mol files. Must be unique & universal
    :type name_field: string
    :return:
    :rtype:
    """
    if not os.path.exists(target):
        os.mkdir(target)

    if mine_db.compounds.find().count() != mine_db.compounds.find({name_field: {'$exists': 1}}).count():
        raise ValueError('%s does not exist for every compound in the database' % name_field)

    for compound in mine_db.compounds.find({'_id': {'$regex': '^C'}}):
        mol = AllChem.MolFromSmiles(compound['SMILES'], True, {'CoA': '*', 'R': "*"})
        if "." in name_field:
            compound[name_field] = utils.get_dotted_field(compound, name_field)
        if isinstance(compound[name_field], list):
            compound[name_field] = ','.join(compound[name_field])
        AllChem.MolToMolFile(mol, os.path.join(target, compound[name_field]+'.mol'))


def import_sdf(mine_db, target,):
    """
    Imports a SDF file as a MINE database
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the SDF file to be loaded
    :return:
    """
    sdf_gen = AllChem.SDMolSupplier(target)
    for mol in sdf_gen:
        mine_db.insert_compound(mol, compound_dict=mol.GetPropsAsDict(), pubchem_db=None, kegg_db=None, modelseed_db=None)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "SDF Imported", "Filepath": target})


def import_smiles(mine_db, target,):
    """
    Imports a smiles file as a MINE database
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the SDF file to be loaded
    :return:
    """
    mols = AllChem.SmilesMolSupplier(target, delimiter='\t', nameColumn=0)
    for mol in mols:
        if mol:
            mine_db.insert_compound(mol, compound_dict=mol.GetPropsAsDict(), pubchem_db=None, kegg_db=None, modelseed_db=None)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "SDF Imported", "Filepath": target})


def import_mol_dir(mine_db, target,):
    """
    Imports a directory of molfiles as a MINE database
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the molfile directory to be loaded
    :return:
    """
    for file in os.listdir(target):
        if ".mol" in file:
            mol = AllChem.MolFromMolFile(target+'/'+file)
            if mol:
                mine_db.insert_compound(mol, compound_dict={"Name": file.strip('.mol'), 'Generation': 0},
                                        pubchem_db=None, kegg_db=None, modelseed_db=None)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "MolFiles Imported", "Filepath": target})

if __name__ == '__main__':
    task = sys.argv[1]
    db_name = sys.argv[2]
    path = sys.argv[3]
    database = MINE(db_name)
    if task == 'export-sdf':
        if len(sys.argv) == 5:
            export_sdf(database, path, int(sys.argv[4]))
        else:
            export_sdf(database, path)
    elif task == 'export-smi':
        if len(sys.argv) == 5:
            export_smiles(database, path, int(sys.argv[4]))
        else:
            export_smiles(database, path)
    elif task == 'export-mol':
        if len(sys.argv) == 5:
            export_mol(database, path, sys.argv[4])
        else:
            export_mol(database, path)
    elif task == 'import-sdf':
        import_sdf(database, path)
    elif task == 'import-smi':
        import_smiles(database, path)
    elif task == 'import-mol':
        import_mol_dir(database, path)
    else:
        print("ERROR: Unrecognised Task")
