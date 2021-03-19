"""Tests for queries.py using pytest."""
# pylint: disable=redefined-outer-name

import json
from os.path import dirname

import pymongo
import pytest
from pymongo.errors import ServerSelectionTimeoutError

from minedatabase import databases, queries
from minedatabase.databases import MINE


try:
    client = pymongo.MongoClient(ServerSelectionTimeoutMS=2000)
    client.server_info()
    del client
    is_mongo = True
except ServerSelectionTimeoutError as err:
    is_mongo = False

valid_db = pytest.mark.skipif(not is_mongo, reason="No MongoDB Connection")


def test_quick_search(test_db, glucose, glucose_id):
    """
    GIVEN a quick search query (e.g. glucose identifiers)
    WHEN quick search is used to search based on that query
    THEN make sure that quick search provides the correct results
    """
    assert glucose not in queries.quick_search(test_db, "WQZGKKKJIJFFOK-UHFFFAOYSA-N")
    assert glucose in queries.quick_search(
        test_db, "InChIKey=WQZGKKKJIJFFOK-GASJEMHNSA-N"
    )
    assert glucose in queries.quick_search(
        test_db, "Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f"
    )
    assert glucose in queries.quick_search(test_db, "917030")
    assert glucose in queries.quick_search(test_db, "cpd00027")
    assert glucose in queries.quick_search(test_db, "C00031")
    assert glucose in queries.quick_search(test_db, "Glucose")
    assert glucose_id in queries.quick_search(
        test_db, "WQZGKKKJIJFFOK-GASJEMHNSA-N", {"_id": 1}
    )


def test_database_query(test_db, glucose, glucose_id):
    """
    GIVEN an andvanced search query (e.g. a MINE id)
    WHEN advanced search is used to search based on that query
    THEN make sure that advanced search provides the correct results
    """
    with pytest.raises(ValueError):
        queries.advanced_search(databases.MINE("admin"), "{'MINE_id': 19160}")
    with pytest.raises(ValueError):
        queries.advanced_search(test_db, "")
    assert queries.advanced_search(test_db, "{'MINE_id': 917030}") == [glucose]
    assert queries.advanced_search(test_db, "{'Names': 'Glucose'}") == [glucose]
    assert queries.advanced_search(test_db, "{'MINE_id': 917030}", {"_id": 1}) == [
        glucose_id
    ]


def test_similarity_search(test_db, test_molfile, glucose):
    """
    GIVEN a similary search query
    WHEN similarity search is used to search based on that query
    THEN make sure the similarity search provides the correct results
    """
    assert (
        len(
            queries.similarity_search(
                test_db,
                "Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cn"
                "c4c(N)ncnc43)C(O)C1O)c1nc(O)nc(O)c1N2",
                0.9,
                100,
            )
        )
        == 8
    )
    result = queries.similarity_search(test_db, test_molfile, 0.5, 100, fp_type="MACCS")
    assert glucose in result
    assert len(result) == 3


def test_substructure_search(test_db, glucose):
    """
    GIVEN a substructure search query
    WHEN substructure search is used to search based on that query
    THEN make sure that substructure search provides the correct results
    """
    result = queries.substructure_search(test_db, "CO", 100)
    assert glucose in result
    assert len(result) == 22
    result = queries.substructure_search(test_db, "O=P(O)(O)O", 100)
    assert len(result) == 15
    assert isinstance(result[0], dict)


def test_structure_search(test_db, test_molfile, glucose):
    """
    GIVEN a structure search query
    WHEN structure search is used to search based on that query
    THEN make sure that structure search provides the correct results
    """
    assert glucose in queries.structure_search(
        test_db, "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", True
    )
    assert glucose in queries.structure_search(test_db, test_molfile, False)
