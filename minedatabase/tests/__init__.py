import json
import os
from ..databases import MINE, establish_db_client


def setup_package():
    print(__name__, '__init__.py : setup_package() ========================================')
    testdb = MINE("mongotest")
    with open(os.path.dirname(__file__)+'/data/testing_db.json') as infile:
        jsondb = json.load(infile)
    for doc in jsondb[0]:
        testdb.compounds.save(doc)
    for doc in jsondb[1]:
        testdb.reactions.save(doc)


def teardown_package():
    print(__name__, '__init__.py : teardown_package() ========================================')
    client = establish_db_client()
    client.drop_database('mongotest')
