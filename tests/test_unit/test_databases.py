"""Tests for databases.py using pytest."""
# pylint: disable=redefined-outer-name
# pylint: disable=protected-access

import json
import os
import shutil
from copy import copy, deepcopy
from pathlib import Path
from shutil import rmtree

import pymongo
import pytest
from pymongo.errors import ServerSelectionTimeoutError

from minedatabase import utils
from minedatabase.databases import (
    MINE,
    write_compounds_to_mine,
    write_core_compounds,
    write_reactions_to_mine,
)


try:
    client = pymongo.MongoClient(ServerSelectionTimeoutMS=10)
    client.server_info()
    del client
except ServerSelectionTimeoutError as err:
    pytest.skip("No Mongo Server Found", allow_module_level=True)

file_path = Path(__file__)
file_dir = file_path.parent

DATA_DIR = (file_dir / "../data/").resolve()


@pytest.fixture()
def cpd_dict():
    """Create a compound dict for testing."""
    cpd_dict = {
        "_id": "Ctest",
        "SMILES": (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OS(=O)(=O)O)[C@@H]"
            "(OP(=O)(O)O)[C@H]1O"
        ),
        "Inchi": (
            "InChI=1S/C10H15N5O13P2S/c11-8-5-9(13-2-12-8)15(3-14-5)10-"
            "6(16)7(27-29(17,18)19)4(26-10)1-25-30(20,21)28-31(22,23)"
            "24/h2-4,6-7,10,16H,1H2,(H,20,21)(H2,11,12,13)(H2,17,18,19"
            ")(H,22,23,24)/t4-,6-,7-,10-/m1/s1"
        ),
        "Type": "Coreactant",
        "Generation": 0,
        "Formula": "C10H15N5O13P2S",
        "Expand": False,
    }

    yield cpd_dict


@pytest.fixture()
def rxn_dicts():
    """Create two reaction dictionaries for testing."""
    rxn_dicts = [
        {
            "_id": "RXN1",
            "Reactants": [[1, "cpd1"], [2, "cpd2"]],
            "Products": [[1, "cpd3"], [1, "cpd4"]],
            "Operators": ["op1", "op2"],
            "SMILES_rxn": "(1) A + (2) B -> (1) C + (2) D",
        },
        {
            "_id": "RXN2",
            "Reactants": [[1, "cpd1"], [2, "cpd2"]],
            "Products": [[1, "cpd3"], [1, "cpd4"]],
            "Operators": ["op1"],
            "SMILES_rxn": "(1) A + (2) B -> (1) C + (2) D",
        },
    ]
    return rxn_dicts


@pytest.mark.skipif(
    os.name == "nt", reason="MolConvert fails on Windows due to permissions errors"
)
def test_generate_image_files(test_db):
    """Test image generation."""
    img_dir = DATA_DIR / "imgs/"
    test_db.generate_image_files(img_dir)
    try:
        assert (img_dir / "C455bc3dc93cd3bb3bef92a34767693a4716aa3fb.svg").is_file()
        assert len(os.listdir(img_dir)) == 26
    finally:
        rmtree(img_dir)
    # Use 'Generation': 1 to do this for just the starting compound
    test_db.generate_image_files(img_dir, {"Generation": 1}, 3, "png")
    try:
        assert (
            img_dir / "C/c/f" / "Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f.png"
        ).is_file()
    finally:
        shutil.rmtree(img_dir)


def test_insert_single_mine_compound(test_db):
    """Test single mine insertion."""
    smiles = "CC(=O)O"
    compound_dict = [{"_id": "test_cpd", "SMILES": smiles}]
    write_compounds_to_mine(compound_dict, test_db)

    entry = test_db.compounds.find_one({"SMILES": smiles})
    assert entry
    assert entry["_id"] == "test_cpd"


def test_insert_bulk_mine_compounds(test_db):
    """Test inserting bulk compounds."""

    smiles1 = "CC(=O)O"
    smiles2 = "CCN"
    compound_dict = [
        {"_id": "test_cpd1", "SMILES": smiles1},
        {"_id": "test_cpd2", "SMILES": smiles2},
    ]

    write_compounds_to_mine(compound_dict, test_db)

    entry1 = test_db.compounds.find_one({"SMILES": smiles1})
    assert entry1
    assert entry1["_id"] == "test_cpd"

    entry2 = test_db.compounds.find_one({"SMILES": smiles2})
    assert entry2
    assert entry2["_id"] == "test_cpd2"


def test_insert_single_core_compound(test_db, cpd_dict):
    """Test inserting a single core compound."""
    write_core_compounds([cpd_dict], test_db, "test")

    try:
        entry = test_db.core_compounds.find_one({"_id": cpd_dict["_id"]})
        assert entry
        assert isinstance(entry["Mass"], float)
        assert entry["Inchi"]
        assert entry["Inchikey"]
    finally:
        test_db.core_compounds.delete_many({"_id": cpd_dict["_id"]})


def test_insert_bulk_core_compound(test_db, cpd_dict):
    """Test inserting many core compounds."""

    cpd_dict2 = copy(cpd_dict)
    cpd_dict2["_id"] = "cpd_dict2"

    write_core_compounds([cpd_dict, cpd_dict2], test_db, "test")
    #
    try:
        for smiles in [cpd_dict["SMILES"], cpd_dict2["SMILES"]]:
            entry = test_db.core_compounds.find_one({"SMILES": smiles})
            assert entry
            assert isinstance(entry["Mass"], float)
            assert entry["Inchi"]
            assert entry["Inchikey"]
    finally:
        for i in [0, 1]:
            test_db.core_compounds.delete_many({"_id": f"test_mine_cpd{i}"})


def test_write_reaction(test_db, rxn_dicts):
    """Test writing reactions."""
    write_reactions_to_mine(rxn_dicts, test_db)

    assert (deepcopy(rxn_dicts[0])) == deepcopy(
        test_db.reactions.find_one({"_id": "RXN1"})
    )
    assert (deepcopy(rxn_dicts[1])) == deepcopy(
        test_db.reactions.find_one({"_id": "RXN2"})
    )


def db_driver(test_db):
    """Test database driver."""
    assert isinstance(test_db.compounds, pymongo.collection.Collection)
    assert isinstance(test_db.reactions, pymongo.collection.Collection)
    assert isinstance(test_db.operators, pymongo.collection.Collection)
    assert isinstance(test_db.core_compounds, pymongo.collection.Collection)
    assert isinstance(test_db._db, pymongo.database.Database)
