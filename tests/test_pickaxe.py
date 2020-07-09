"""Tests for pickaxe.py using pytest."""
# pylint: disable=redefined-outer-name
# pylint: disable=protected-access

import os
import re
import hashlib
import subprocess
from filecmp import cmp

import pytest

from rdkit.Chem import AllChem
from minedatabase import pickaxe
from minedatabase.databases import MINE

DATA_DIR = os.path.dirname(__file__) + '/data'


@pytest.fixture
def smiles_dict():
    """Store SMILES for compounds used in test cases here."""
    smiles = {
        'ATP': 'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C'
               + '@@H](O)[C@H]1O',
        'ADP': 'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C'
               + '@H]1O',
        'meh': 'CCC(=O)C(=O)O',
        'l_ala': 'C[C@H](N)C(=O)O',
        'd_ala': 'C[C@@H](N)C(=O)O',
        'FADH': 'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc'
                + '4c(N)ncnc43)C(O)C1O)c1nc(O)nc(O)c1N2',
        'S-Adenosylmethionine': 'C[S+](CC[C@H](N)C(=O)O)C[C@H]1O[C@@H](n2cnc'
                                + '3c(N)ncnc32)[C@H](O)[C@@H]1O'
    }
    return smiles


@pytest.fixture
def coreactant_dict(smiles_dict):
    """Create tab-formatted coreactant entries."""
    coreactants = {
        'ATP': 'ATP		' + smiles_dict['ATP'],
        'ADP': 'ADP		' + smiles_dict['ADP'],
        'S-Adenosylmethionine': 'S-Adenosylmethionine		'
                                + smiles_dict['S-Adenosylmethionine']
    }
    return coreactants


@pytest.fixture
def pk():
    """Create default Pickaxe object."""
    return pickaxe.Pickaxe(coreactant_list=DATA_DIR + '/test_coreactants.tsv',
                           rule_list=DATA_DIR + '/test_reaction_rules.tsv')


@pytest.fixture
def default_rule(pk):
    """Default operator to use for testing is set here."""
    return pk.operators['2.7.1.a']


@pytest.fixture
def pk_transformed(default_rule, smiles_dict, coreactant_dict):
    """Create Pickaxe object with a few predicted reactions."""
    pk_transformed = pickaxe.Pickaxe()
    pk_transformed._add_compound("Start", smi=smiles_dict['FADH'])
    pk_transformed._load_coreactant(coreactant_dict['ATP'])
    pk_transformed._load_coreactant(coreactant_dict['ADP'])
    pk_transformed.operators['2.7.1.a'] = default_rule
    pk_transformed.transform_compound(smiles_dict['FADH'])
    pk_transformed.assign_ids()
    return pk_transformed


def purge(directory, pattern):
    """Delete all files in a directory matching a regex pattern."""
    for filename in os.listdir(directory):
        if re.search(pattern, filename):
            os.remove(os.path.join(directory, filename))


def multiprocess(pk, smiles_dict, coreactant_dict):
    """Convert FADH into other compounds in parallel."""
    pk._load_coreactant(coreactant_dict['ATP'])
    pk._load_coreactant(coreactant_dict['ADP'])
    pk._add_compound(smiles_dict['FADH'], smiles_dict['FADH'],
                     cpd_type='Starting Compound')
    pk.transform_all(max_generations=2, num_workers=2)
    return pk
    
def delete_database(name):
    mine = MINE(name)
    mine.client.drop_database(name)
    mine.client.close()
    

def test_cofactor_loading(pk):
    """
    GIVEN a default Pickaxe object
    WHEN cofactors are loaded into the Pickaxe object in its creation
    THEN make sure those cofactors were loaded correctly
    """
    assert "O=C=O" in pk._raw_compounds
    c_id = pk._raw_compounds['O=C=O']
    assert c_id in pk.compounds
    assert pk.compounds[c_id]['Type'] == 'Coreactant'
    assert isinstance(pk.coreactants['ATP'][0], AllChem.Mol)
    assert pk.coreactants['ATP'][1][0] == "X"


def test_reaction_rule_loading(default_rule):
    """
    GIVEN a reaction rule dict
    WHEN reaction rules are loaded during Pickaxe object initialization
    THEN make sure it is formatted correctly
    """
    assert isinstance(default_rule[0], AllChem.ChemicalReaction)
    assert isinstance(default_rule[1], dict)
    assert default_rule[1]['Reactants'] == ['ATP', 'Any']
    assert "Products" in default_rule[1]
    assert "Comments" in default_rule[1]


def test_compound_loading(pk):
    """
    GIVEN a default Pickaxe object
    WHEN compounds are loaded
    THEN check that they are loaded correctly
    """
    compound_smiles = pk.load_compound_set(
        compound_file=DATA_DIR + '/test_compounds.tsv')
    assert len(compound_smiles) == 14


def test_transform_compounds_implicit(smiles_dict):
    """
    GIVEN a compound (meh in this case)
    WHEN that compound is transformed into another compound via pickaxe
    THEN make sure that the transformation is successful and recorded
    """
    pk = pickaxe.Pickaxe(explicit_h=False, kekulize=False,
                         coreactant_list=DATA_DIR + '/test_coreactants.tsv',
                         rule_list=DATA_DIR + '/test_cd_rxn_rule.tsv')
    pk._add_compound("Start", smi=smiles_dict['meh'])
    pk.transform_compound(smiles_dict['meh'])
    assert len(pk.compounds) == 38
    assert len(pk.reactions) == 1


def test_hashing(pk, smiles_dict, coreactant_dict):
    """
    GIVEN a default Pickaxe object
    WHEN compounds with very similar structures are both transformed
    THEN make sure they result in different products (hashed differently)
    """
    pk._load_coreactant(coreactant_dict['S-Adenosylmethionine'])
    pk.transform_compound(smiles_dict['l_ala'])
    len_rxns = len(pk.reactions)
    assert len_rxns > 0
    pk.transform_compound(smiles_dict['d_ala'])
    # If hashed differently, should now have an extra set of rxns for d_ala
    assert len(pk.reactions) == 2 * len_rxns


def test_product_racimization(smiles_dict):
    """
    GIVEN molecules of undefined steriochemistry
    WHEN transforming them via pickaxe with or without racemization
    THEN make sure that when racemized, we observe extra reaction(s)
    """
    pk = pickaxe.Pickaxe(racemize=False,
                         coreactant_list=DATA_DIR + '/test_coreactants.tsv',
                         rule_list=DATA_DIR + '/test_reaction_rules.tsv')
    comps, rxns = pk.transform_compound(smiles_dict['meh'], rules=['2.6.1.a'])
    assert len(comps) == 38
    assert len(rxns) == 1

    pk = pickaxe.Pickaxe(racemize=True,
                         coreactant_list=DATA_DIR + '/test_coreactants.tsv',
                         rule_list=DATA_DIR + '/test_reaction_rules.tsv')
    rcomps, rrxns = pk.transform_compound(smiles_dict['meh'],
                                          rules=['2.6.1.a'])
    # Extra reaction occurs due to undefined chiral center giving rise to
    # two possible starting structures in this second case
    assert len(rcomps) == 39
    assert len(rrxns) == 2


def test_compound_output_writing(pk_transformed):
    """
    GIVEN a Pickaxe object with predicted transformations
    WHEN all compounds (including predicted) are written to an output file
    THEN make sure they are correctly written, and that they are all present
    """
    with open(DATA_DIR + '/testcompoundsout.tsv', 'rb') as infile:
        expected = hashlib.sha256(infile.read()).hexdigest()
    pk_transformed.write_compound_output_file(DATA_DIR
                                              + '/testcompoundsout.tsv')
    assert os.path.exists(DATA_DIR + '/testcompoundsout_new.tsv')
    try:
        with open(DATA_DIR + '/testcompoundsout_new.tsv', 'rb') as infile:
            output_compounds = hashlib.sha256(infile.read()).hexdigest()
        assert expected == output_compounds
    finally:
        os.remove(DATA_DIR + '/testcompoundsout_new.tsv')


def test_reaction_output_writing(pk_transformed):
    """
    GIVEN a Pickaxe object with predicted transformations
    WHEN all reactions (including predicted) are written to an output file
    THEN make sure they are correctly written, and that they are all present
    """
    with open(DATA_DIR + '/testreactionsout.tsv', 'rb') as infile:
        expected = hashlib.sha256(infile.read()).hexdigest()
    pk_transformed.write_reaction_output_file(DATA_DIR
                                              + '/testreactionsout.tsv')
    assert os.path.exists(DATA_DIR + '/testreactionsout_new.tsv')
    try:
        with open(DATA_DIR + '/testreactionsout_new.tsv', 'rb') as infile:
            output_compounds = hashlib.sha256(infile.read()).hexdigest()
        assert expected == output_compounds
    finally:
        os.remove(DATA_DIR + '/testreactionsout_new.tsv')


def test_transform_all(default_rule, smiles_dict, coreactant_dict):
    """
    GIVEN a set of rules and starting compounds
    WHEN we run pickaxe to predict potential transformations
    THEN make sure all expected transformations are predicted
    """
    pk = pickaxe.Pickaxe(errors=False)
    pk._load_coreactant(coreactant_dict['ATP'])
    pk._load_coreactant(coreactant_dict['ADP'])
    pk._add_compound(smiles_dict['FADH'], smiles_dict['FADH'],
                     cpd_type='Starting Compound')
    pk.operators['2.7.1.a'] = default_rule
    pk.transform_all(max_generations=2)
    assert len(pk.compounds) == 31
    assert len(pk.reactions) == 49
    comp_gens = set([x['Generation'] for x in pk.compounds.values()])
    assert comp_gens == {0, 1, 2}


def test_multiprocessing(pk, smiles_dict, coreactant_dict):
    """
    GIVEN a Pickaxe object
    WHEN we use multiprocessing to enumerate predicted reactions
    THEN make sure those predictions are correct
    """
    pk = multiprocess(pk, smiles_dict, coreactant_dict)
    assert len(pk.compounds) == 104
    assert len(pk.reactions) == 155
    comp_gens = set([x['Generation'] for x in pk.compounds.values()])
    assert comp_gens == {0, 1, 2, 3}


def test_cli():
    """
    GIVEN the pickaxe CLI
    WHEN pickaxe is run from the command line
    THEN make sure it exits with exit code 0 (no errors)
    """
    os.chdir(DATA_DIR + "/../..")
    rc = subprocess.call(
        'python minedatabase/pickaxe.py -o tests -r '
        'tests/data/test_cd_rxn_rule.tsv',
        shell=True)
    assert not rc
    purge('tests/', r".*\.tsv$")


def test_pruning(default_rule, smiles_dict, coreactant_dict):
    """
    GIVEN a Pickaxe expansion
    WHEN that expansion is pruned via Pickaxe.prune_network()
    THEN make sure that the pruned compounds no longer exist in the network
    """

    pk = pickaxe.Pickaxe(database=None, image_dir=None)
    pk.operators['2.7.1.a'] = default_rule
    pk = multiprocess(pk, smiles_dict, coreactant_dict)
    ids = ['C9437bf42d165907392564b23c0ca132b8bd51625',
           'C70016088b2c54458296232054a1ec58d54035560', 'C41']
    pk.prune_network(ids)
    pk.assign_ids()
    pk.write_compound_output_file(DATA_DIR + '/pruned_comps')
    pk.write_reaction_output_file(DATA_DIR + '/pruned_rxns')
    assert os.path.exists(DATA_DIR + '/pruned_comps_new')
    assert os.path.exists(DATA_DIR + '/pruned_rxns_new')
    try:
        assert cmp(DATA_DIR + '/pruned_comps', DATA_DIR + '/pruned_comps_new')
    finally:
        os.remove(DATA_DIR + '/pruned_comps_new')
    try:
        assert cmp(DATA_DIR + '/pruned_rxns', DATA_DIR + '/pruned_rxns_new')
    finally:
        os.remove(DATA_DIR + '/pruned_rxns_new')


def test_database_already_exists(default_rule, smiles_dict, coreactant_dict):
    """
    GIVEN an existing MINE
    WHEN a new pickaxe object is defined
    THEN make sure program exits with database collision
    """
    delete_database('MINE_test')
    pk = pickaxe.Pickaxe(database='MINE_test')
    pk.operators['2.7.1.a'] = default_rule
    pk = multiprocess(pk, smiles_dict, coreactant_dict)
    pk.save_to_mine()

    try:     
        with pytest.raises(SystemExit) as pytest_wrapped_e:
                pk = pickaxe.Pickaxe(database='MINE_test')
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 'Exiting due to database name collision.'
    finally:
        delete_database('MINE_test')


def test_save_as_mine(default_rule, smiles_dict, coreactant_dict):
    """
    GIVEN a Pickaxe expansion
    WHEN that expansion is saved as a MINE DB in the MongoDB
    THEN make sure that all features are saved in the MongoDB as expected
    """
    delete_database('MINE_test')
    pk = pickaxe.Pickaxe(database='MINE_test', image_dir=DATA_DIR)
    pk.operators['2.7.1.a'] = default_rule
    pk = multiprocess(pk, smiles_dict, coreactant_dict)
    pk.save_to_mine()
    mine_db = MINE('MINE_test')
    try:
        assert mine_db.compounds.estimated_document_count() == 31
        assert mine_db.reactions.estimated_document_count() == 49
        assert mine_db.operators.estimated_document_count() == 1
        assert mine_db.operators.find_one()["Reactions_predicted"] == 49
        assert os.path.exists(
            DATA_DIR + '/X02fba734e41145959768095d202750b2c777b274.svg')
        start_comp = mine_db.compounds.find_one({'Type': 'Starting Compound'})
        assert len(start_comp['Reactant_in']) > 0
        # Don't track sources of coreactants
        coreactant = mine_db.compounds.find_one({'Type': 'Coreactant'})
        assert 'Product_of' not in coreactant
        assert 'Reactant_in' not in coreactant
        product = mine_db.compounds.find_one({'Generation': 2})
        assert len(product['Product_of']) > 0
        assert product['Type'] == 'Predicted'
    finally:
        delete_database('MINE_test')
        purge(DATA_DIR, r".*\.svg$")


def test_load_compounds_from_mine(default_rule, smiles_dict, coreactant_dict):
    """
    GIVEN a Pickaxe object loaded from the 'MINE_test' MINE DB
    WHEN Pickaxe.load_compound_set() is used to load all compounds from the DB
    THEN make sure that all compounds are loaded
    """
    # Generate database
    delete_database('MINE_test')
    pk = pickaxe.Pickaxe(database='MINE_test')
    pk.operators['2.7.1.a'] = default_rule
    pk = multiprocess(pk, smiles_dict, coreactant_dict)
    pk.save_to_mine()
    del(pk)
    try:
        pk = pickaxe.Pickaxe()
        compound_smiles = pk.load_compound_set(database='MINE_test')
        assert len(compound_smiles) == 1
    finally:
        delete_database('MINE_test')
        

def test_save_no_rxn_mine():
    """
    GIVEN a Pickaxe object with no expansion
    WHEN that Pickaxe object is saved into a MINE DB in the MongoDB
    THEN check that starting compounds are present and that no reactions exist
    """
    pk = pickaxe.Pickaxe(database='MINE_test')
    pk.load_compound_set(compound_file=DATA_DIR + '/test_compounds.tsv')
    pk.save_to_mine()
    mine_db = MINE('MINE_test')
    try:
        assert mine_db.compounds.estimated_document_count() == 14
        assert mine_db.reactions.estimated_document_count() == 0
    finally:
       delete_database('MINE_test')
