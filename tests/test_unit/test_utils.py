"""Tests for utils.py using pytest."""

import os
from collections import OrderedDict
from pathlib import Path

from rdkit.DataStructs.cDataStructs import ExplicitBitVect

import minedatabase.utils as utils


file_path = Path(__file__)
file_dir = file_path.parent

DATA_DIR = (file_dir / "../data/").resolve()


def test_get_compound_hash():
    """Test compound to hash."""
    assert utils.get_compound_hash("CCO", "Coreactant") == (
        "Xa41fe8492d86f214ba494e3d04da2f0854c0e2ea",
        "LFQSCWFLJHTTHZ",
    )
    assert utils.get_compound_hash("CCO", "Predicted") == (
        "Ca41fe8492d86f214ba494e3d04da2f0854c0e2ea",
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )
    assert utils.get_compound_hash("CCO", "Starting Compound") == (
        "Ca41fe8492d86f214ba494e3d04da2f0854c0e2ea",
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    )


def test_get_fp():
    smiles = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
    assert type(utils.get_fp(smiles)) == ExplicitBitVect


def test_get_compound_hash_two_blocks():
    smiles = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
    assert utils.get_compound_hash(smiles, "Starting Compound", inchi_blocks=2) == (
        "Cf95a3c17f908e427c3127b4e8c3d8575c286d6ce",
        "WQZGKKKJIJFFOK-DVKNGEFBSA-N",
    )

    assert utils.get_compound_hash(smiles, "Starting Compound", inchi_blocks=1) == (
        "C9ab1a08d72c90a8167d1f3a668d8f1138e534a07",
        "WQZGKKKJIJFFOK-DVKNGEFBSA-N",
    )


def test_file_to_dict_list():
    """
    GIVEN compound files of different types (.tsv, .csv, and .json)
    WHEN the file contents are converted to a list of compound dicts
    THEN check that the list of compound dicts is produced as expected
    """
    res_7 = OrderedDict(
        [
            ("id", "cpd01211"),
            ("abbreviation", "tcynt"),
            ("name", "Thiocyanate"),
            ("formula", "CNS"),
            ("mass", "58"),
            ("source", "ModelSEED"),
            ("structure", "InChI=1S/CHNS/c2-1-3/h3H"),
            ("charge", "-1"),
            ("is_core", "1"),
            ("is_obsolete", "0"),
            ("linked_compound", "null"),
            ("is_cofactor", "0"),
            ("deltag", "22.2"),
            ("deltagerr", "5.68687"),
            ("pka", "3:0.5"),
            ("pkb", ""),
            ("abstract_compound", "null"),
            ("comprised_of", "null"),
            ("aliases", "null"),
        ]
    )

    filenames = ["test_compounds.tsv", "test_compounds.csv", "test_compounds.json"]

    for file in filenames:
        res = utils.file_to_dict_list(DATA_DIR / file)
        assert len(res) == 15
        assert res[7] == res_7
