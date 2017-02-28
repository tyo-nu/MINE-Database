import pymongo
from .. import databases
from .. import utils
from rdkit.Chem import AllChem

test_db = databases.MINE('mongotest')


def test_insert_compound():
    smiles = 'CC(=O)O'
    mol = AllChem.MolFromSmiles(smiles)
    test_db.insert_compound(mol, {'Generation': 0.0})
    try:
        entry = test_db.compounds.find_one({"SMILES": smiles})
        assert entry
        assert isinstance(entry['Mass'], float)
        assert len(entry['RDKit'])
        assert len(entry['RDKit']) == entry['len_RDKit']
        assert len(entry['Names'])
        assert entry["NP_likeness"]
        assert entry['logP']
    finally:
        test_db.compounds.remove({"SMILES": smiles})


def test_insert_reaction():
    rxn = {'Equation': 'A + B = C + D'}
    rxn['Reactants'], rxn['Products'] = utils.parse_text_rxn(rxn['Equation'], ' = ', ' + ')
    test_db.insert_reaction(rxn)
    entry = test_db.reactions.find_one({"_id": "4542c96f4bca04bfe2db15bc71e9eaee38bee5b87ad8a6752a5c4718ba1974c1"})
    assert entry
    assert isinstance(entry['_id'], str) and len(entry['_id'])
    assert entry['Reactants'] == [{"stoich": 1, "c_id": "A"}, {"stoich": 1, "c_id": "B"}]
    assert entry['Products'] == [{"stoich": 1, "c_id": "C"}, {"stoich": 1, "c_id": "D"}]


def test_init():
    assert isinstance(test_db.compounds, pymongo.collection.Collection)
    assert isinstance(test_db.reactions, pymongo.collection.Collection)
    assert isinstance(test_db.operators, pymongo.collection.Collection)
    assert isinstance(test_db._db, pymongo.database.Database)
    assert isinstance(test_db.id_db, pymongo.database.Database)

