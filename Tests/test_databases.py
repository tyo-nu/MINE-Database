import databases
from rdkit.Chem import AllChem
import os

test_db = databases.MINE('mongotest')


def test_insert_compound():
    smiles = 'CC(=O)O'
    mol = AllChem.MolFromSmiles(smiles)
    test_db.insert_compound(mol, {'Generation': 0.0})
    try:
        entry = test_db.compounds.find_one({"SMILES": smiles})
        assert entry
        assert 'C00033' in entry['DB_links']['KEGG']
        assert 'cpd00029' in entry['DB_links']['Model_SEED']
        assert '3335' in entry['DB_links']['PubChem']
        assert isinstance(entry['Mass'], float)
        assert len(entry['RDKit'])
        assert len(entry['RDKit']) == entry['len_RDKit']
        assert len(entry['Names'])
        assert entry["NP_likeness"]
        assert entry['logP']
    finally:
        test_db.compounds.remove({"SMILES": smiles})


def test_add_rxn_pointers():
    raise NotImplementedError


def test_add_compound_sources():
    raise NotImplementedError
