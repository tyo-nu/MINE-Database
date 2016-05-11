__author__ = 'JGJeffryes'
import queries
import databases

test_db = databases.MINE('mongotest')
test_molfile = open("./Tests/glucose.mol", "r").read()
glucose = {'Names': ['Hexose', 'D-Idose', 'Glucose', 'Mannose', 'D-Gulose', 'D-Allose', 'D-Hexose', 'Dextrose',
                     'Seminose', 'L-Gulose', 'D-Talose', 'D-Aldose', 'D-Mannose', 'D-Aldose2', 'D-Aldose1', 'D-Glucose',
                     'D-Altrose', 'Carubinose', 'Grape sugar', 'L-Galactose', 'D-Galactose', 'D-ido-Hexose',
                     'D-gulo-Hexose', 'D-talo-Hexose', 'beta-D-Mannose', 'beta-D-Glucose', 'D-altro-Hexose',
                     'alpha-D-Glucose', 'alpha-D-Mannose', 'D-glucopyranose', 'beta-D-Galactose', 'alpha-D-Galactose',
                     'D-galactopyranose', '1,4-beta-D-Mannooligosaccharide'], 'SMILES': 'OCC1OC(O)C(C(C1O)O)O',
           'Generation': 0.0, 'Mass': 180.063388104, '_id': 'Cb5b3273ab083d77ed29fbef8f7e464929af29c13',
           'MINE_id': 19160, 'NP_likeness': 0, 'Formula': 'C6H12O6', 'Inchikey': 'WQZGKKKJIJFFOK-UHFFFAOYSA-N'}
glucose_id = {'_id': 'Cb5b3273ab083d77ed29fbef8f7e464929af29c13'}

def test_quick_search():
    assert glucose in queries.quick_search(test_db, 'WQZGKKKJIJFFOK-UHFFFAOYSA-N')
    assert glucose in queries.quick_search(test_db, 'C00031')
    assert glucose in queries.quick_search(test_db, 'Glucose')
    assert glucose_id in queries.quick_search(test_db, 'WQZGKKKJIJFFOK-UHFFFAOYSA-N', {'_id': 1})


def test_database_query():
    assert queries.advanced_search('admin', '') == ['Illegal query']
    assert queries.advanced_search(test_db, "{'MINE_id': 19160}") == [glucose]
    assert queries.advanced_search(test_db, "{'Names': 'Glucose'}") == [glucose]
    assert queries.advanced_search(test_db, "{'MINE_id': 19160}", {'_id': 1}) == [glucose_id]

def test_similarity_search():
    assert len(queries.similarity_search(test_db, 'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc'
                                                  '43)C(O)C1O)c1nc(O)nc(O)c1N2', 0.9, "RDKit", 100)) == 8
    assert len(queries.similarity_search(test_db, test_molfile, 0.5, 'MACCS', 100)) == 2

def test_substructure_search():
    meh = queries.substructure_search(test_db, 'O=P(O)(O)O', 100)
    assert len(meh) == 15
    assert isinstance(meh[0], dict)

def test_structure_search():
    assert glucose in queries.structure_search(test_db, 'OCC1OC(O)C(C(C1O)O)O')
    assert glucose in queries.structure_search(test_db, test_molfile)
