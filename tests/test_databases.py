"""Tests for databases.py using pytest."""
# pylint: disable=redefined-outer-name
# pylint: disable=protected-access

import json
import os
import unittest
from shutil import rmtree

import pymongo
import pytest
from pymongo.errors import ServerSelectionTimeoutError
from rdkit.Chem import AllChem

from minedatabase import utils

from minedatabase.databases import MINE


@pytest.fixture()
def test_db():
    """Create a test MINE database. Created and torn down before and after each
    test it is used in."""
    print(os.path.dirname(__file__))
    datafile_path = os.path.join(os.path.dirname(__file__),
                                 'data/testing_db.json')
    delete_database("mongotest")
    try:
        testdb = MINE("mongotest")
        with open(datafile_path) as infile:
            jsondb = json.load(infile)
        for doc in jsondb[0]:
            if testdb.compounds.find_one({'_id': doc['_id']}):
                testdb.compounds.replace_one({'_id': doc['_id']}, doc)
            else:
                testdb.compounds.insert_one(doc)
        for doc in jsondb[1]:
            if testdb.reactions.find_one({'_id': doc['_id']}):
                testdb.reactions.replace_one({'_id': doc['_id']}, doc)
            else:
                testdb.reactions.insert_one(doc)
        for doc in jsondb[2]:
            if testdb.operators.find_one({'_id': doc['_id']}):
                testdb.operators.replace_one({'_id': doc['_id']}, doc)
            else:
                testdb.operators.insert_one(doc)

    except ServerSelectionTimeoutError:
        print('No Mongo DB server detected')

    yield testdb

def delete_database(name):
    mine = MINE(name)
    mine.client.drop_database(name)
    mine.client.close()


@unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                 "Skipping this test on Travis CI.")
def test_generate_image_files(test_db):
    """
    GIVEN a MINE database
    WHEN image files are generated from that database
    THEN make sure that they are all generated correctly
    """
    img_dir = os.path.join(os.path.dirname(__file__), 'imgs')
    test_db.generate_image_files(img_dir)
    try:
        assert os.path.exists(
            os.path.join(img_dir,
                         './C455bc3dc93cd3bb3bef92a34767693a4716aa3fb.svg'))
        assert len(os.listdir(img_dir)) == 26
    finally:
        rmtree(img_dir)
    # Use 'Generation': 1 to do this for just the starting compound
    test_db.generate_image_files(img_dir, {'Generation': 1}, 3, 'png')
    try:
        assert os.path.exists(
            os.path.join(img_dir, 'C', 'c', 'f',
                         'Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f.png'))
    finally:
        rmtree(img_dir)


def test_insert_single_mine_compound(test_db):
    """
    GIVEN a pickaxe object with compounds in it
    WHEN those compounds and their properties are added to a MINE database
    THEN check that the compound and its properties are correctly stored
    """
    
    smiles = 'CC(=O)O'
    compound_dict = {'_id': 'test_cpd', 'SMILES': smiles}
    test_db.insert_mine_compound(compound_dict, requests=None) 
    
    try:
        entry = test_db.compounds.find_one({'SMILES': smiles})
        assert entry
        assert entry['_id'] == 'test_cpd'
    finally:
        delete_database("mongotest")

def test_insert_bulk_mine_compounds(test_db):
    """
    GIVEN multiple compound informations
    WHEN those compounds and their properties are added to a MINE database in bulk
    THEN check that the compound and its properties are correctly stored
    """
    
    smiles1 = 'CC(=O)O'
    smiles2 = 'CCN'
    compound_dicts = [{'_id': 'test_cpd1', 'SMILES': smiles1},
                      {'_id': 'test_cpd2', 'SMILES': smiles2}]

    requests = []
    for compound_dict in compound_dicts:
        test_db.insert_mine_compound(compound_dict, requests=requests)

    test_db.compounds.bulk_write(requests, ordered=False)
    
    try:
        entry1 = test_db.compounds.find_one({'SMILES': smiles1})
        assert entry1
        assert entry1['_id'] == 'test_cpd1'

        entry2 = test_db.compounds.find_one({'SMILES': smiles2})
        assert entry2
        assert entry2['_id'] == 'test_cpd2'
    finally:
        delete_database("mongotest")

def test_insert_single_core_compound(test_db):
    """
    GIVEN a compound (Mol object) with associated properties
    WHEN that compound and its properties are added to a MINE database
    THEN check that the compound and its properties are correctly stored
    """
    smiles = "CC(=O)O"
    mol = AllChem.MolFromSmiles(smiles)
    test_db.insert_core_compound(mol, 'test_mine_cpd', requests=None)
    # test_db._core_db.bulk_write(insert_request, ordered=False)  
    #   
    try:
        entry = test_db.core_compounds.find_one({'_id': 'test_mine_cpd'})
        assert entry
        assert isinstance(entry['Mass'], float)
        assert entry['Inchi']
        assert entry['Inchikey']
        assert entry['MINES']
        assert entry["NP_likeness"]
        assert entry['logP']
    finally:
        test_db.core_compounds.delete_many({'_id': 'test_mine_cpd'})

def test_insert_bulk_core_compound(test_db):
    """
    GIVEN a compound (Mol object) with associated properties
    WHEN that compound and its properties are added to a MINE database
    THEN check that the compound and its properties are correctly stored
    """
    smiles1 = 'CC(=O)O'
    smiles2 = 'CCN'

    mol1 = AllChem.MolFromSmiles(smiles1)
    mol2 = AllChem.MolFromSmiles(smiles2)

    requests = []
    for i, mol in enumerate([mol1, mol2]):
        test_db.insert_core_compound(mol, f'test_mine_cpd{i}', requests=requests)
    
    test_db.core_compounds.bulk_write(requests, ordered=False)  
    #   
    try:
        for smiles in [smiles1, smiles2]:
            entry = test_db.core_compounds.find_one({'SMILES': smiles})
            assert entry
            assert isinstance(entry['Mass'], float)
            assert entry['Inchi']
            assert entry['Inchikey']
            assert entry['MINES']
            assert entry["NP_likeness"]
            assert entry['logP']
    finally:
        for i in [0, 1]:
            test_db.core_compounds.delete_many({'_id': f'test_mine_cpd{i}'})

def test_insert_reaction(test_db):
    """
    GIVEN a reaction (dict with 'Equation', 'Reactants', and 'Products' keys)
    WHEN that reaction is inserted into a MINE database
    THEN check that it is present and that the stoichiometry is correct
    """
    rxn = {'Equation': 'A + B = C + D'}
    rxn['Reactants'], rxn['Products'] = utils.parse_text_rxn(rxn['Equation'],
                                                             ' = ', ' + ')
    test_db.insert_reaction(rxn)
    entry = test_db.reactions.find_one({"_id": "R4542c96f4bca04bfe2db15bc71e9"
                                               "eaee38bee5b87ad8a6752a5c4718"
                                               "ba1974c1"})
    assert entry
    assert isinstance(entry['_id'], str) and (len(entry['_id']) > 0)
    assert entry['Reactants'] == [{"stoich": 1, "c_id": "A"},
                                  {"stoich": 1, "c_id": "B"}]
    assert entry['Products'] == [{"stoich": 1, "c_id": "C"},
                                 {"stoich": 1, "c_id": "D"}]


def test_init(test_db):
    """
    GIVEN a MINE database
    WHEN N/A
    THEN make sure that all database field types are correct
    """
    assert isinstance(test_db.compounds, pymongo.collection.Collection)
    assert isinstance(test_db.reactions, pymongo.collection.Collection)
    assert isinstance(test_db.operators, pymongo.collection.Collection)
    assert isinstance(test_db.core_compounds, pymongo.collection.Collection)
    assert isinstance(test_db._db, pymongo.database.Database)
