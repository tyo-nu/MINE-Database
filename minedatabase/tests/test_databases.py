import pymongo
import unittest
from minedatabase import databases
from minedatabase import utils
from rdkit.Chem import AllChem
import os
from shutil import rmtree

test_db = databases.MINE('mongotest')


@unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                 "Skipping this test on Travis CI.")
def test_generate_image_files():
    img_dir = os.getcwd() + '/imgs'
    test_db.generate_image_files(img_dir)
    try:
        assert os.path.exists(
            os.path.join(img_dir,
                         './C455bc3dc93cd3bb3bef92a34767693a4716aa3fb.svg'))
        print(len(os.listdir(img_dir)))
        assert len(os.listdir(img_dir)) == 26
    finally:
        rmtree(img_dir)
    test_db.generate_image_files(img_dir, {'Generation': 1}, 3, 'png')
    try:
        assert os.path.exists(
            os.path.join(img_dir, 'C', 'c', 'f',
                         'Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f.png'))
    finally:
        rmtree(img_dir)


def test_insert_compound():
    smiles = 'CC(=O)O'
    mol = AllChem.MolFromSmiles(smiles)
    test_db.insert_compound(mol, {'Generation': 0.0})
    try:
        entry = test_db.compounds.find_one({"SMILES": smiles})
        assert entry
        assert isinstance(entry['Mass'], float)
        assert len(entry['RDKit'])
        assert len(entry['RDKit']) == entry['len_RDKit']
        assert entry["NP_likeness"]
        assert entry['logP']
    finally:
        test_db.compounds.remove({"SMILES": smiles})


def test_insert_reaction():
    rxn = {'Equation': 'A + B = C + D'}
    rxn['Reactants'], rxn['Products'] = utils.parse_text_rxn(rxn['Equation'], ' = ', ' + ')
    test_db.insert_reaction(rxn)
    entry = test_db.reactions.find_one({"_id": "4542c96f4bca04bfe2db15bc71e9eaee38bee5b87ad8a6752a5c4718ba1974c1"})
    assert entry
    assert isinstance(entry['_id'], str) and len(entry['_id'])
    assert entry['Reactants'] == [{"stoich": 1, "c_id": "A"}, {"stoich": 1, "c_id": "B"}]
    assert entry['Products'] == [{"stoich": 1, "c_id": "C"}, {"stoich": 1, "c_id": "D"}]


def test_init():
    assert isinstance(test_db.compounds, pymongo.collection.Collection)
    assert isinstance(test_db.reactions, pymongo.collection.Collection)
    assert isinstance(test_db.operators, pymongo.collection.Collection)
    assert isinstance(test_db._db, pymongo.database.Database)
    assert isinstance(test_db.id_db, pymongo.database.Database)

