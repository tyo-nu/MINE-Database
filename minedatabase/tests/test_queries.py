from os.path import dirname
from nose.tools import assert_raises

from .. import databases, queries

data_dir = dirname(__file__)+'/data'

test_db = databases.MINE('mongotest')
test_molfile = open(data_dir+"/glucose.mol", "r").read()
glucose = {u'SMILES': u'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', u'Inchikey': u'WQZGKKKJIJFFOK-GASJEMHNSA-N',
           u'Generation': 1, u'MINE_id': 917030, u'NP_likeness': 2.62691337083175,
           u'Sources': [{u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'C2ce917f1d3aaef7501595123894a68a7d786a9e7']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C2c3b118fb2dbd237cc8a4878b0643385860fb427',
                                                                    u'C08a914cde05039694ef0194d9ee79ff9a79dde33']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Ca67dd0794223d284a9517566c3fd4107728e0808']},
                        {u'Operators': [u'3.1.1.b'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'C71c9c2296f22b1aa1bd031d91e391fa415dc72fd']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'C1f719e8cf42a266bc93fe185816e244f5498fc1e']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'Cd51c567a382bff0cba43ecc97d5807e0ffebb5f8',
                                                                    u'C08a914cde05039694ef0194d9ee79ff9a79dde33']},
                        {u'Operators': [u'1.1.-1.h'], u'Compounds': [u'C20fa60742147ed7b7c328d7835b84d14036d7f36',
                                                                     u'C1b424c6768e3787b3453fa5d572b06106b8933d4',
                                                                     u'Cac741137c0d633e53486ed958334b206a85dcc03']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C8d3efae75c074cd20730fbc4752c0498d94ba40e',
                                                                    u'C08a914cde05039694ef0194d9ee79ff9a79dde33']},
                        {u'Operators': [u'1.1.1.g'], u'Compounds': [u'C6d0b4d05908b6e7fd3f71b56ca0e81c7a9b9027d',
                                                                    u'C8e823fcc54dec371d81eebaa8e01dc4796e4add6']},
                        {u'Operators': [u'3.1.3.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Cff41aa31cc2f8f2f13759ff3d4283726265072ee']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Cd51de949ca3a35a5b1bdd3959225766ffd2fc027']},
                        {u'Operators': [u'3.2.1.a'], u'Compounds': [u'C08a914cde05039694ef0194d9ee79ff9a79dde33',
                                                                    u'Cfdd7350d4ca8405dd44990548ae18d5d56364003']},
                        {u'Operators': [u'1.14.13.d'], u'Compounds': [u'C20fa60742147ed7b7c328d7835b84d14036d7f36',
                                                                      u'C71c9c2296f22b1aa1bd031d91e391fa415dc72fd',
                                                                      u'C1aa818910461f5961959eb05dd06b8f7a4fbbb9f',
                                                                      u'Cac741137c0d633e53486ed958334b206a85dcc03']}],
           u'Mass': 180.063388104, u'Names': [u'Glucose', u'Dextrose', u'D-Glucose', u'Grape sugar', u'D-Glucopyranose'],
           u'Formula': u'C6H12O6', u'_id': u'Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f'}
glucose_id = {'_id': 'Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f'}


def test_quick_search():
    assert glucose not in queries.quick_search(test_db, 'WQZGKKKJIJFFOK-UHFFFAOYSA-N')
    assert glucose in queries.quick_search(test_db, 'C00031')
    assert glucose in queries.quick_search(test_db, 'Glucose')
    assert glucose_id in queries.quick_search(test_db, 'WQZGKKKJIJFFOK-GASJEMHNSA-N', {'_id': 1})


def test_database_query():
    with assert_raises(ValueError):
        queries.advanced_search(databases.MINE('admin'), "{'MINE_id': 19160}")
    with assert_raises(ValueError):
        queries.advanced_search(test_db, "")
    assert queries.advanced_search(test_db, "{'MINE_id': 917030}") == [glucose]
    assert queries.advanced_search(test_db, "{'Names': 'Glucose'}") == [glucose]
    assert queries.advanced_search(test_db, "{'MINE_id': 917030}", {'_id': 1}) == [glucose_id]


def test_similarity_search():
    assert len(queries.similarity_search(test_db, 'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc'
                                                  '43)C(O)C1O)c1nc(O)nc(O)c1N2', 0.9, "RDKit", 100)) == 8
    foo = queries.similarity_search(test_db, test_molfile, 0.5, 'MACCS', 100)
    assert glucose in foo
    assert len(foo) == 3


def test_substructure_search():
    foo = queries.substructure_search(test_db, "CO", 100)
    assert glucose in foo
    assert len(foo) == 22
    bar = queries.substructure_search(test_db, 'O=P(O)(O)O', 100)
    assert len(bar) == 15
    assert isinstance(bar[0], dict)


def test_structure_search():
    assert glucose in queries.structure_search(test_db, 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', True)
    assert glucose in queries.structure_search(test_db, test_molfile, False)
