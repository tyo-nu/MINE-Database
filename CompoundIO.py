#!/usr/bin/env python
"""Databases.py: Functions to load MINE databases from and dump compounds into common cheminformatics formats"""

__author__ = 'JGJeffryes'

from rdkit import Chem
from Databases import MINE
import sys


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
    w = Chem.SDWriter(target)
    for compound in mine_db.compounds.find():
        mol = Chem.MolFromSmiles(compound['SMILES'])
        mol.SetProp('_id', compound['_id'])
        if 'Reactant_in' in compound:
            mol.SetProp('Reactant_in', compound['Reactant_in'])
        if 'Reactant_in' in compound:
            mol.SetProp('Product_of', compound['Product_of'])
        w.write(mol)
    w.close()


if __name__ == '__main__':
    task = sys.argv[1]
    db_name = sys.argv[2]
    path = sys.argv[3]
    database = MINE(db_name)
    if task == 'export-sdf':
        export_sdf(database, path)
