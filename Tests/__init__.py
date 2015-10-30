__author__ = 'JGJeffryes'
from databases import MINE
from databases import establish_db_client
import json

def setup_package():
    print(__name__, '__init__.py : setup_package() ========================================')
    testdb = MINE("mongotest")
    with open('Tests/testing_db.json') as infile:
        comps = json.load(infile)
    testdb.compounds.insert(comps)


def teardown_package():
    print(__name__, '__init__.py : teardown_package() ========================================')
    client = establish_db_client()
    client.drop_database('mongotest')
