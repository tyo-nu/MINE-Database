__author__ = 'JGJeffryes'

import pickaxe
import rdkit
import filecmp
import os
from databases import MINE

pk = pickaxe.Pickaxe()
meh = 'CCC(=O)C(=O)O'
fadh = 'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1nc(O)nc(O)c1N2'

def test_cofactor_loading():
    pk2 = pickaxe.Pickaxe(cofactor_list='Tests/Cofactor_SMILES.tsv')
    assert "O=C=O" in pk2._raw_compounds
    assert isinstance(pk2.cofactors['ATP'], rdkit.Chem.rdchem.Mol)

def test_reaction_rule_loading():
    pk2 = pickaxe.Pickaxe(rule_list='Tests/test_operators.tsv')
    rule = pk2.rxn_rules['2.7.1.a']
    assert isinstance(rule[1], rdkit.Chem.rdChemReactions.ChemicalReaction)
    assert "Any" in rule[0]

def test_compound_loading():
    raise NotImplementedError

def test_transform_compounds():
    pk._load_cofactor('ATP	Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O')
    pk.load_rxn_rule('2.7.1.a	ATP;Any	[#6;H2D4:8][#8;H0D2:7][#15;H0D4:6][#8;H0D2:5][#15;H0D4:4][#8;H0D2:3]'
                     '[#15;H0D4:2][#8;H1D2R0:1].[#1;D1R0:11][#8;H1D2R0:10][#6:9]>>[*:1]-[*:2]-[*:10]-[*:9].[*:8]-[*:7]'
                     '-[*:6]-[*:5]-[*:4]-[*:3]-[*:11]')
    pk.transform_compound(fadh)

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
    pk.write_compound_output_file('Tests/testcompounds')
    try:
        assert os.path.exists('Tests/testcompounds_new')
        assert filecmp.cmp('Tests/testcompounds', 'Tests/testcompounds_new')
    finally:
        os.remove('Tests/testcompounds_new')

def test_reaction_output_writing():
    pk.write_reaction_output_file('Tests/testreactions')
    assert os.path.exists('Tests/testreactions_new')
    try:
        assert filecmp.cmp('Tests/testreactions', 'Tests/testreactions_new')
    finally:
        os.remove('Tests/testreactions_new')

def test_save_as_MINE():
    pk.save_to_MINE("MINE_test")
    mine_db = MINE('MINE_test')
    try:
        assert mine_db.compounds.count() == 10
        assert mine_db.reactions.count() == 7
    finally:
        mine_db.compounds.drop()
        mine_db.reactions.drop()
