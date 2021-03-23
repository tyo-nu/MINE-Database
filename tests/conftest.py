"""Fixtures used by tests from different modules go here."""

import json
from pathlib import Path

import pytest
from pymongo.errors import ServerSelectionTimeoutError

from minedatabase.databases import MINE
from minedatabase.pickaxe import Pickaxe


file_path = Path(__file__)
file_dir = file_path.parent


def delete_database(name):
    """Delete a mine database."""
    mine = MINE(name)
    mine.client.drop_database(name)
    mine.client.close()


@pytest.fixture(scope="module")
def test_db():
    """Create a test MINE database for testing."""
    datafile_path = file_dir / "data/testing_db.json"
    try:
        testdb = MINE("mongotest")
        with open(datafile_path) as infile:
            jsondb = json.load(infile)
        for doc in jsondb[0]:
            if testdb.compounds.find_one({"_id": doc["_id"]}):
                testdb.compounds.replace_one({"_id": doc["_id"]}, doc)
            else:
                testdb.compounds.insert_one(doc)
        for doc in jsondb[1]:
            if testdb.reactions.find_one({"_id": doc["_id"]}):
                testdb.reactions.replace_one({"_id": doc["_id"]}, doc)
            else:
                testdb.reactions.insert_one(doc)
        for doc in jsondb[2]:
            if testdb.operators.find_one({"_id": doc["_id"]}):
                testdb.operators.replace_one({"_id": doc["_id"]}, doc)
            else:
                testdb.operators.insert_one(doc)

    except ServerSelectionTimeoutError:
        print("No Mongo DB server detected")

    yield testdb
    delete_database("mongotest")


@pytest.fixture(scope="function")
def pk():
    """Create default Pickaxe object."""
    return Pickaxe(
        coreactant_list=file_dir / "data/test_coreactants.tsv",
        rule_list=file_dir / "data/test_reaction_rules.tsv",
        explicit_h=True,
        quiet=False,
    )


# test pickaxe
@pytest.fixture
def smiles_dict():
    """Store SMILES for compounds used in test cases here."""
    smiles = {
        "ATP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C"
        + "@@H](O)[C@H]1O",
        "ADP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C" + "@H]1O",
        "meh": "CCC(=O)C(=O)O",
        "l_ala": "C[C@H](N)C(=O)O",
        "d_ala": "C[C@@H](N)C(=O)O",
        "FADH": "Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc"
        + "4c(N)ncnc43)C(O)C1O)c1nc(O)nc(O)c1N2",
        "S-Adenosylmethionine": "C[S+](CC[C@H](N)C(=O)O)C[C@H]1O[C@@H](n2cnc"
        + "3c(N)ncnc32)[C@H](O)[C@@H]1O",
    }
    return smiles


@pytest.fixture
def coreactant_dict(smiles_dict):
    """Create tab-formatted coreactant entries."""
    coreactants = {
        "ATP": "ATP		" + smiles_dict["ATP"],
        "ADP": "ADP		" + smiles_dict["ADP"],
        "S-Adenosylmethionine": "S-Adenosylmethionine		"
        + smiles_dict["S-Adenosylmethionine"],
    }
    return coreactants


@pytest.fixture
def default_rule(pk):
    """Set default operator."""
    return pk.operators["2.7.1.a"]


@pytest.fixture
def pk_transformed(default_rule, smiles_dict, coreactant_dict):
    """Create Pickaxe object with a few predicted reactions."""
    pk_transformed = Pickaxe(explicit_h=True)
    pk_transformed._add_compound(
        "Start", smi=smiles_dict["FADH"], cpd_type="Starting Compound"
    )
    pk_transformed._load_coreactant(coreactant_dict["ATP"])
    pk_transformed._load_coreactant(coreactant_dict["ADP"])
    pk_transformed.operators["2.7.1.a"] = default_rule
    pk_transformed.transform_all()
    pk_transformed.assign_ids()
    return pk_transformed


# Queries
@pytest.fixture()
def test_molfile():
    """Mol file for glucose compound."""
    test_molfile = open((file_dir / "data/glucose.mol"), "r").read()
    return test_molfile


@pytest.fixture()
def glucose():
    """MongoDB document (.json) for glucose compound."""
    print(file_dir / "/data/glucose.json")
    with open((file_dir / "data/glucose.json")) as infile:
        glucose = json.load(infile)
    return glucose


@pytest.fixture
def glucose_id():
    """ID in MongoDB for glucose."""
    glucose_id = {"_id": "Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f"}
    return glucose_id


# Filters
@pytest.fixture(scope="function")
def pk_target():
    """Generate pickaxe to expand and filter."""
    new_pk = Pickaxe(
        rule_list="./tests/data/test_filters/test_filter_rules.tsv",
        coreactant_list="./tests/data/test_filters/metacyc_coreactants.tsv",
        filter_after_final_gen=False,
    )
    new_pk.load_compound_set("./tests/data/test_filters/test_filter_compounds.csv")
    new_pk.load_targets("./tests/data/test_filters/test_filter_targets.csv")

    return new_pk
