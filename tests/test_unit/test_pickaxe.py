"""Tests for pickaxe.py using pytest."""
# pylint: disable=redefined-outer-name
# pylint: disable=protected-access

import hashlib
import os
import re
import subprocess
from filecmp import cmp
from pathlib import Path

import pymongo
import pytest
from pymongo.errors import ServerSelectionTimeoutError
from rdkit.Chem import AllChem

from minedatabase import pickaxe
from minedatabase.databases import MINE


try:
    client = pymongo.MongoClient(ServerSelectionTimeoutMS=20)
    client.server_info()
    del client
    is_mongo = True
except ServerSelectionTimeoutError as err:
    is_mongo = False
valid_db = pytest.mark.skipif(not is_mongo, reason="No MongoDB Connection")

file_path = Path(__file__)
file_dir = file_path.parent

DATA_DIR = (file_dir / "../data/").resolve()


def purge(directory, pattern):
    """Delete all files in a directory matching a regex pattern."""
    for filename in os.listdir(directory):
        if re.search(pattern, filename):
            os.remove(os.path.join(directory, filename))


def delete_database(name):
    """Delete database."""
    mine = MINE(name)
    mine.client.drop_database(name)
    mine.client.close()


def test_cofactor_loading(pk):
    """Test loading cofactors.

    GIVEN a default Pickaxe object
    WHEN cofactors are loaded into the Pickaxe object in its creation
    THEN make sure those cofactors were loaded correctly
    """
    c_id = "X73bc8ef21db580aefe4dbc0af17d4013961d9d17"

    assert c_id in pk.compounds
    assert pk.compounds[c_id]["Formula"] == "H2O"
    assert pk.compounds[c_id]["Type"] == "Coreactant"
    assert isinstance(pk.coreactants["Water"][0], AllChem.Mol)
    assert pk.coreactants["Water"][1][0] == "X"


def test_reaction_rule_loading(default_rule):
    """Test loading rules.

    GIVEN a reaction rule dict
    WHEN reaction rules are loaded during Pickaxe object initialization
    THEN make sure it is formatted correctly
    """
    assert isinstance(default_rule[0], AllChem.ChemicalReaction)
    assert isinstance(default_rule[1], dict)
    assert default_rule[1]["Reactants"] == ["ATP", "Any"]
    assert "Products" in default_rule[1]
    assert "Comments" in default_rule[1]


def test_compound_loading(pk):
    """Test loading compounds.

    GIVEN a default Pickaxe object
    WHEN compounds are loaded
    THEN check that they are loaded correctly
    """
    compound_smiles = pk.load_compound_set(
        compound_file=file_dir / "../data/test_compounds.tsv"
    )
    assert len(compound_smiles) == 14


def test_transform_all(default_rule, smiles_dict, coreactant_dict):
    """Test transform function.

    GIVEN a set of rules and starting compounds
    WHEN we run pickaxe to predict potential transformations
    THEN make sure all expected transformations are predicted
    """
    pk = pickaxe.Pickaxe(errors=False, explicit_h=True)
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound(
        smiles_dict["FADH"], smiles_dict["FADH"], cpd_type="Starting Compound"
    )
    pk.operators["2.7.1.a"] = default_rule
    pk.transform_all(generations=2)
    assert len(pk.compounds) == 31
    assert len(pk.reactions) == 49
    comp_gens = set([x["Generation"] for x in pk.compounds.values()])
    assert comp_gens == {0, 1, 2}


def test_compound_output_writing(pk_transformed):
    """Test compound output writing.

    GIVEN a Pickaxe object with predicted transformations
    WHEN all compounds (including predicted) are written to an output file
    THEN make sure they are correctly written, and that they are all present
    """

    with open(file_dir / "../data/testcompoundsout.tsv", "rb") as infile:
        expected = hashlib.sha256(infile.read()).hexdigest()
    pk_transformed.write_compound_output_file(file_dir / "../data/testcompoundsout.tsv")
    assert os.path.exists(file_dir / "../data/testcompoundsout_new.tsv")
    try:
        with open(file_dir / "../data/testcompoundsout_new.tsv", "rb") as infile:
            output_compounds = hashlib.sha256(infile.read()).hexdigest()
        assert expected == output_compounds
    finally:
        os.remove(file_dir / "../data/testcompoundsout_new.tsv")


def test_reaction_output_writing(pk_transformed):
    """Test writing reaction output.

    GIVEN a Pickaxe object with predicted transformations
    WHEN all reactions (including predicted) are written to an output file
    THEN make sure they are correctly written, and that they are all present
    """
    with open(file_dir / "../data/testreactionsout.tsv", "rb") as infile:
        expected = hashlib.sha256(infile.read()).hexdigest()
    pk_transformed.write_reaction_output_file(file_dir / "../data/testreactionsout.tsv")
    assert os.path.exists(file_dir / "../data/testreactionsout_new.tsv")
    try:
        with open(file_dir / "../data/testreactionsout_new.tsv", "rb") as infile:
            output_compounds = hashlib.sha256(infile.read()).hexdigest()
        assert expected == output_compounds
    finally:
        os.remove(file_dir / "../data/testreactionsout_new.tsv")


def test_multiprocessing(pk, smiles_dict, coreactant_dict):
    """Test multiprocessing.

    GIVEN a Pickaxe object
    WHEN we use multiprocessing to enumerate predicted reactions
    THEN make sure those predictions are correct
    """
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound("FADH", smiles_dict["FADH"], cpd_type="Starting Compound")
    pk.transform_all(generations=2, processes=2)

    assert len(pk.compounds) == 67
    assert len(pk.reactions) == 49
    comp_gens = set([x["Generation"] for x in pk.compounds.values()])
    assert comp_gens == {0, 1, 2}


def test_pruning(default_rule, smiles_dict, coreactant_dict):
    """Test pruning network to targets.

    GIVEN a Pickaxe expansion
    WHEN that expansion is pruned via Pickaxe.prune_network()
    THEN make sure that the pruned compounds no longer exist in the network
    """

    pk = pickaxe.Pickaxe(explicit_h=True)
    pk.operators["2.7.1.a"] = default_rule
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound("FADH", smiles_dict["FADH"], cpd_type="Starting Compound")
    pk.transform_all(generations=2)

    ids = [
        "C89d19c432cbe8729c117cfe50ff6ae4704a4e6c1",
        "C750e93db23dd3f796ffdf9bdefabe32b10710053",
        "C41",
    ]
    pk.prune_network(ids)
    pk.assign_ids()

    DATA_DIR = (file_dir / "../data").resolve()
    pk.write_compound_output_file(DATA_DIR / "pruned_comps.tsv")
    pk.write_reaction_output_file(DATA_DIR / "pruned_rxns.tsv")
    assert os.path.exists(DATA_DIR / "pruned_comps_new.tsv")
    assert os.path.exists(DATA_DIR / "pruned_rxns_new.tsv")
    try:
        assert cmp(DATA_DIR / "pruned_comps.tsv", DATA_DIR / "pruned_comps_new.tsv")
        assert cmp(DATA_DIR / "pruned_rxns.tsv", DATA_DIR / "pruned_rxns_new.tsv")
    finally:
        os.remove((DATA_DIR / "pruned_comps_new.tsv").resolve())
        os.remove((DATA_DIR / "pruned_rxns_new.tsv").resolve())


def test_target_generation(default_rule, smiles_dict, coreactant_dict):
    """Test generating a target from starting compounds."""
    pk = pickaxe.Pickaxe(explicit_h=True)
    pk.operators["2.7.1.a"] = default_rule
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound("FADH", smiles_dict["FADH"], cpd_type="Starting Compound")
    pk.load_targets(file_dir / "../data/test_targets.csv")
    pk.transform_all(generations=2)
    pk.prune_network_to_targets()

    assert "C11088915f64b93293e70af9c3b7822a4f131225d" in pk.compounds
    assert len(pk.reactions) == 4
    assert len(pk.compounds) == 6


@valid_db
def test_save_as_mine(default_rule, smiles_dict, coreactant_dict):
    """Test saving compounds to database.

    GIVEN a Pickaxe expansion
    WHEN that expansion is saved as a MINE DB in the MongoDB
    THEN make sure that all features are saved in the MongoDB as expected
    """
    DATA_DIR = (file_dir / "../data").resolve()
    delete_database("MINE_test")
    pk = pickaxe.Pickaxe(database="MINE_test", image_dir=DATA_DIR, explicit_h=True)
    pk.operators["2.7.1.a"] = default_rule
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound("FADH", smiles_dict["FADH"], cpd_type="Starting Compound")
    pk.transform_all(generations=2)
    pk.save_to_mine(processes=1)
    mine_db = MINE("MINE_test")

    try:

        assert mine_db.compounds.estimated_document_count() == 31
        assert mine_db.reactions.estimated_document_count() == 49
        assert mine_db.operators.estimated_document_count() == 1
        assert mine_db.operators.find_one()["Reactions_predicted"] == 49
        assert os.path.exists(
            DATA_DIR / "X9c29f84930a190d9086a46c344020283c85fb917.svg"
        )
        start_comp = mine_db.compounds.find_one({"Type": "Starting Compound"})
        assert len(start_comp["Reactant_in"]) > 0
        # Don't track sources of coreactants
        coreactant = mine_db.compounds.find_one({"Type": "Coreactant"})
        assert "Product_of" not in coreactant
        assert "Reactant_in" not in coreactant
        product = mine_db.compounds.find_one({"Generation": 2})
        assert len(product["Product_of"]) > 0
        assert product["Type"] == "Predicted"
    finally:
        delete_database("MINE_test")
        purge(DATA_DIR, r".*\.svg$")


@valid_db
def test_save_target_mine(default_rule, smiles_dict, coreactant_dict):
    """Test saving the target run to a MINE."""
    delete_database("MINE_test")
    pk = pickaxe.Pickaxe(database="MINE_test", explicit_h=True)
    pk.operators["2.7.1.a"] = default_rule
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound("FADH", smiles_dict["FADH"], cpd_type="Starting Compound")
    pk.load_targets(file_dir / "../data/test_targets.csv")
    pk.transform_all(generations=2)
    pk.prune_network_to_targets()

    pk.save_to_mine()
    mine_db = MINE("MINE_test")

    try:
        assert mine_db.compounds.estimated_document_count() == 6
        assert mine_db.reactions.estimated_document_count() == 4
        assert mine_db.operators.estimated_document_count() == 1
        assert mine_db.operators.find_one()["Reactions_predicted"] == 4
        start_comp = mine_db.target_compounds.find_one()
        assert start_comp["InChI_key"] == "RYNUDNWPSBJQQY-UHFFFAOYSA-N"
        assert all([i in start_comp.keys() for i in ["_id", "SMILES", "InChI_key"]])
    finally:
        delete_database("MINE_test")


@valid_db
def test_database_already_exists(default_rule, smiles_dict, coreactant_dict):
    """Test database collision.

    GIVEN an existing MINE
    WHEN a new pickaxe object is defined
    THEN make sure program exits with database collision
    """
    delete_database("MINE_test")
    pk = pickaxe.Pickaxe(database="MINE_test")
    pk.operators["2.7.1.a"] = default_rule
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound("FADH", smiles_dict["FADH"], cpd_type="Starting Compound")
    pk.transform_all(generations=2)
    pk.save_to_mine(processes=1)

    try:
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            pk = pickaxe.Pickaxe(database="MINE_test")
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == (
            "Exiting due to database name collision."
        )
    finally:
        delete_database("MINE_test")


def test_pickle(coreactant_dict, smiles_dict, default_rule):
    """Test pickling of pickaxe objects."""
    pickle_path = Path("test_pickle.pk")

    pk = pickaxe.Pickaxe(errors=False, explicit_h=True)
    pk._load_coreactant(coreactant_dict["ATP"])
    pk._load_coreactant(coreactant_dict["ADP"])
    pk._add_compound(
        smiles_dict["FADH"], smiles_dict["FADH"], cpd_type="Starting Compound"
    )
    pk.operators["2.7.1.a"] = default_rule

    pk.transform_all(generations=2)

    pk.pickle_pickaxe(pickle_path)
    del pk
    pk = pickaxe.Pickaxe(errors=False)
    pk.load_pickled_pickaxe(pickle_path)

    assert len(pk.compounds) == 31
    assert len(pk.reactions) == 49
    comp_gens = set([x["Generation"] for x in pk.compounds.values()])
    assert comp_gens == {0, 1, 2}

    pickle_path.unlink()


def test_local_cli():
    """Test command line interface writing locally.

    GIVEN the pickaxe CLI
    WHEN pickaxe is run from the command line
    THEN make sure it exits with exit code 0 (no errors)
    """
    os.chdir(file_dir / "../data/../..")
    rc = subprocess.call(
        f"python minedatabase/pickaxe.py -o tests/ -r tests/data/test_cd_rxn_rule.tsv",
        shell=True,
    )
    assert not rc
    purge(DATA_DIR / "..", r".*\.tsv$")


@valid_db
def test_mongo_cli():
    """Test command line interface writing to mongo."""
    mine = MINE("tests")
    os.chdir(file_dir / "../data/../..")
    rc = subprocess.call(
        "python minedatabase/pickaxe.py -d tests -r tests/data/test_cd_rxn_rule.tsv",
        shell=True,
    )
    assert not rc

    try:
        assert mine.compounds.estimated_document_count() == 51
    finally:
        mine.client.drop_database("tests")

    purge(file_dir / "..", r".*\.svg$")
