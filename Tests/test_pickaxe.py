__author__ = 'JGJeffryes'

import pickaxe
import rdkit
import filecmp
import os
import re
from databases import MINE

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))

pk = pickaxe.Pickaxe(image_dir="Tests/", database='MINE_test')
meh = 'CCC(=O)C(=O)O'
l_ala = 'C[C@H](N)C(=O)O'
d_ala = 'C[C@@H](N)C(=O)O'
fadh = 'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1nc(O)nc(O)c1N2'

def test_cofactor_loading():
    pk2 = pickaxe.Pickaxe(cofactor_list='Tests/Cofactor_SMILES.tsv')
    assert "O=C=O" in pk2._raw_compounds
    assert isinstance(pk2.cofactors['ATP'], rdkit.Chem.rdchem.Mol)

def test_reaction_rule_loading():
    pk2 = pickaxe.Pickaxe(rule_list='Tests/test_operators.tsv')
    rule = pk2.rxn_rules['2.7.1.a']
    assert isinstance(rule[0], rdkit.Chem.rdChemReactions.ChemicalReaction)
    assert isinstance(rule[1], dict)
    assert rule[1]['Reactants'] == ['Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O', 'Any']

def test_compound_loading():
    compound_smiles = pk.load_compound_set(compound_file='Tests/test_compounds.tsv')
    assert len(compound_smiles) == 15
    pk2 = pickaxe.Pickaxe(database='mongotest')
    compound_smiles = pk2.load_compound_set()
    assert len(compound_smiles) == 26

def test_transform_compounds():
    pk._load_cofactor('ATP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk._load_cofactor('ADP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk.load_rxn_rule('2.7.1.a	ATP;Any	[#6;H2D4:8][#8;H0D2:7][#15;H0D4:6][#8;H0D2:5][#15;H0D4:4][#8;H0D2:3]'
                     '[#15;H0D4:2][#8;H1D2R0:1].[#1;D1R0:11][#8;H1D2R0:10][#6:9]>>[*:1]-[*:2]-[*:10]-[*:9].[*:8]-[*:7]'
                     '-[*:6]-[*:5]-[*:4]-[*:3]-[*:11]')
    pk.transform_compound(fadh)
    pk._assign_ids()

def test_hashing():
    pk2 = pickaxe.Pickaxe(explicit_h=False, kekulize=False)
    pk2._load_cofactor('S-Adenosylmethionine	C[S+](CC[C@H](N)C(=O)O)C[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O')
    #pk2.load_rxn_rule('Methylation_1	Any;S-Adenosylmethionine	[#6;R0:2]-[#7;R0;H1,H2:1].[#6H3:3]-[#16;D3:4]>>'
    #                  '[#6:2]-[#7:1]-[#6:3].[#16:4]	Any;S-Adenosylhomocystine')
    pk2.transform_compound(l_ala)
    len_rxns = len(pk2.reactions)
    pk2.transform_compound(d_ala)
    assert len(pk2.reactions) == 2 * len_rxns

def test_product_racimization():
    pk2 = pickaxe.Pickaxe(raceimze=False, rule_list='Tests/test_operators.tsv')
    comps, rxns = pk2.transform_compound(meh, rules=['2.6.1.a'])
    assert len(comps) == 2
    assert len(rxns) == 1
    pk2 = pickaxe.Pickaxe(raceimze=True, rule_list='Tests/test_operators.tsv')
    rcomps, rrxns = pk2.transform_compound(meh, rules=['2.6.1.a'])
    assert len(rcomps) == 3
    assert len(rrxns) == 2

def test_compound_output_writing():
    pk.write_compound_output_file('Tests/testcompoundsout')
    assert os.path.exists('Tests/testcompoundsout_new')
    try:
        assert filecmp.cmp('Tests/testcompoundsout', 'Tests/testcompoundsout_new')
    finally:
        os.remove('Tests/testcompoundsout_new')

def test_reaction_output_writing():
    pk.write_reaction_output_file('Tests/testreactionsout')
    assert os.path.exists('Tests/testreactionsout_new')
    try:
        assert filecmp.cmp('Tests/testreactionsout', 'Tests/testreactionsout_new')
    finally:
        os.remove('Tests/testreactionsout_new')

def test_transform_all():
    pk3 = pickaxe.Pickaxe(errors=False)
    pk3.compounds[meh] = {'ID': None, '_id': fadh, 'Inchikey': '', 'SMILES': fadh, 'Generation': 0}
    pk3._load_cofactor('ATP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk3._load_cofactor('ADP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk3.generation = 0
    pk3.load_rxn_rule('2.7.1.a	ATP;Any	[#6;H2D4:8][#8;H0D2:7][#15;H0D4:6][#8;H0D2:5][#15;H0D4:4][#8;H0D2:3]'
                     '[#15;H0D4:2][#8;H1D2R0:1].[#1;D1R0:11][#8;H1D2R0:10][#6:9]>>[*:1]-[*:2]-[*:10]-[*:9].[*:8]-[*:7]'
                     '-[*:6]-[*:5]-[*:4]-[*:3]-[*:11]')
    pk3.transform_all(max_generations=2)
    assert len(pk3.compounds) == 32
    assert len(pk3.reactions) == 49


def test_multiprocessing():
    pk3 = pickaxe.Pickaxe(errors=False)
    pk3.compounds[meh] = {'ID': None, '_id': fadh, 'Inchikey': '', 'SMILES': fadh, 'Generation': 0}
    pk3._load_cofactor('ATP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk3._load_cofactor('ADP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk3.generation = 0
    pk3.load_rxn_rule('2.7.1.a	ATP;Any	[#6;H2D4:8][#8;H0D2:7][#15;H0D4:6][#8;H0D2:5][#15;H0D4:4][#8;H0D2:3]'
                     '[#15;H0D4:2][#8;H1D2R0:1].[#1;D1R0:11][#8;H1D2R0:10][#6:9]>>[*:1]-[*:2]-[*:10]-[*:9].[*:8]-[*:7]'
                     '-[*:6]-[*:5]-[*:4]-[*:3]-[*:11]')
    pk3.transform_all(max_generations=2, num_workers=2)
    assert len(pk3.compounds) == 32
    assert len(pk3.reactions) == 49


def test_save_as_MINE():
    pk.save_to_MINE("MINE_test")
    mine_db = MINE('MINE_test')
    try:
        assert mine_db.compounds.count() == 25
        assert mine_db.reactions.count() == 7
        assert mine_db.operators.count() == 1
        assert os.path.exists("Tests/C9c69cbeb40f083118c1913599c12c7f4e5e68d03.svg")
    finally:
        mine_db.compounds.drop()
        mine_db.reactions.drop()
        mine_db.operators.drop()
        purge('Tests', ".*\.svg$")