#!/usr/bin/env python
"""Databases.py: Functions to load MINE databases from and dump compounds into common cheminformatics formats"""

__author__ = 'JGJeffryes'

from pymongo import DESCENDING
from rdkit.Chem import AllChem
from databases import MINE
import utils
import sys
import hashlib
import os
import datetime

def export_sdf(mine_db, target):
    """
    Exports compounds from the database as an MDL SDF file
    :param mine_db: a MINE object
    :param path: directory for files
    :param max_compounds: maximum number of compounds per file (defaults to 100000)
    :return:
    """
    if not mine_db.compounds.find_one({"Product_of": {'$exists': 1}}):
        mine_db.add_rxn_pointers()

    print("Exporting %s compounds from %s as an SDF file") %(mine_db.compounds.count(), mine_db.name)
    target = utils.prevent_overwrite(target)
    w = AllChem.SDWriter(target)
    for compound in mine_db.compounds.find():
        mol = AllChem.MolFromSmiles(compound['SMILES'], True, {'CoA': '*', 'R': "*"})
        mol.SetProp('_id', compound['_id'])
        mol.SetProp('Generation', compound['Generation'])
        if 'Reactant_in' in compound:
            mol.SetProp('Reactant_in', compound['Reactant_in'])
        if 'Reactant_in' in compound:
            mol.SetProp('Product_of', compound['Product_of'])
        w.write(mol)
    w.close()


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
        compound = {}
        for key in mol.GetPropNames():
            compound[key] = mol.GetProp(key)
        mine_db.insert_compound(mol, compound)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "SDF Imported", "Filepath": target})


if __name__ == '__main__':
    task = sys.argv[1]
    db_name = sys.argv[2]
    path = sys.argv[3]
    database = MINE(db_name)
    if task == 'export-sdf':
        export_sdf(database, path)
    if task == 'export-mol':
        if len(sys.argv) == 5:
            export_mol(database, path, sys.argv[4])
        else:
            export_mol(database, path)
    if task == 'import-sdf':
        import_sdf(database, path)
