"""Tests for databases.py using pytest."""
# pylint: disable=redefined-outer-name
# pylint: disable=protected-access

import json
import os
import unittest
from shutil import rmtree
from copy import copy

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

@pytest.fixture()
def cpd_dict():
    cpd_dict = {'_id':'Xtest',
                'SMILES':'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OS(=O)(=O)O)[C@@H](OP(=O)(O)O)[C@H]1O',
                'Inchi':'InChI=1S/C10H15N5O13P2S/c11-8-5-9(13-2-12-8)15(3-14-5)10-6(16)7(27-29(17,18)19)4(26-10)1-25-30(20,21)28-31(22,23)24/h2-4,6-7,10,16H,1H2,(H,20,21)(H2,11,12,13)(H2,17,18,19)(H,22,23,24)/t4-,6-,7-,10-/m1/s1',
                'Type':'Coreactant',
                'Generation':0,
                'Formula':'C10H15N5O13P2S',
                'Expand':False,
                }

    yield cpd_dict

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

def test_insert_single_core_compound(test_db, cpd_dict):
    """
    GIVEN a compound (Mol object) with associated properties
    WHEN that compound and its properties are added to a MINE database
    THEN check that the compound and its properties are correctly stored
    """
    
    test_db.insert_core_compound(cpd_dict, requests=None)

    try:
        entry = test_db.core_compounds.find_one({'_id': cpd_dict['_id']})
        assert entry
        assert isinstance(entry['Mass'], float)
        assert entry['Inchi']
        assert entry['Inchikey']
        assert entry["NP_likeness"]
        assert entry['logP']
    finally:
        test_db.core_compounds.delete_many({'_id': cpd_dict['_id']})

def test_insert_bulk_core_compound(test_db, cpd_dict):
    """
    GIVEN a compound (Mol object) with associated properties
    WHEN that compound and its properties are added to a MINE database
    THEN check that the compound and its properties are correctly stored
    """

    cpd_dict2 = copy(cpd_dict)
    cpd_dict2['_id'] = 'cpd_dict2'

    requests = []
    for i, cpd in enumerate([cpd_dict, cpd_dict2]):
        test_db.insert_core_compound(cpd, requests=requests)
    
    test_db.core_compounds.bulk_write(requests, ordered=False)  
    #   
    try:
        for smiles in [cpd_dict['SMILES'], cpd_dict2['SMILES']]:
            entry = test_db.core_compounds.find_one({'SMILES': smiles})
            assert entry
            assert isinstance(entry['Mass'], float)
            assert entry['Inchi']
            assert entry['Inchikey']
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
    # TODO: update test to new format of prod/react
    rxn = {'Equation': 'A + B = C + D'}
    rxn['Reactants'], rxn['Products'] = utils.parse_text_rxn(rxn['Equation'],
                                                             ' = ', ' + ')
    rxn['_id'] = 'R_testrxn'
    test_db.insert_reaction(rxn)
    entry = test_db.reactions.find_one({"_id": "R_testrxn"})
    assert entry
    assert isinstance(entry['_id'], str) and (len(entry['_id']) > 0)
    assert entry['Reactants'] == [[1, "A"], [1, "B"]]
    assert entry['Products'] == [[1, "C"], [1, "D"]]


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
