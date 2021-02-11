"""reaction.py: This module contains methods to carry out reactions by specifying
reaction rules and compounds. """
import collections
import copy
import multiprocessing
from functools import partial

from rdkit.RDLogger import logger
from rdkit.Chem.AllChem import (AddHs, Kekulize, MolFromSmiles, MolToSmiles,
                                RemoveHs, SanitizeMol, CalcMolFormula) 

from minedatabase import utils

from typing import Tuple
lg = logger()
lg.setLevel(0)

###############################################################################
# Functions to run transformations
# There are two distinct functions to transform two flavors of operators.
#   1. Full Operators are the operators as loaded directly from the list
#       of operatores. These operators use a single supplied molecule for all
#       "Any" instances in the rule.
#   2. Partial operators are operators for reactions with more than one
#       "Any" in the reactants. These rules are derived from individually
#       mapped reactions and are called partial because only one "Any" is
#       allowed to be novel, the other "Any"s are determined by the
#       mapped reactions.

# Both operators are preprocessed slightly differently, but yield the same
# output format back to the pickaxe object.

# there are way too many variables passed... switch to **kwargs?
# Generic reaction implementation
def _run_reaction(rule_name: str, rule: tuple, reactant_mols: dict, coreactant_mols: dict,
                  coreactant_dict: dict, local_cpds: dict, local_rxns: dict, generation: int,
                  explicit_h: bool) -> None:
    """Transform list of mols and a reaction rule into a half reaction
    describing either the reactants or products.
    """
    def _make_half_rxn(mol_list, rules):
        cpds = {}
        cpd_counter = collections.Counter()

        for mol, rule in zip(mol_list, rules):
            if rule == 'Any':
                cpd_dict = _gen_compound(mol)
                # failed compound
                if cpd_dict == None:
                    return None, None
            else:
                cpd_id = coreactant_mols[rule][1]
                cpd_dict = coreactant_dict[cpd_id]

            cpds[cpd_dict['_id']] = cpd_dict
            cpd_counter.update({cpd_dict['_id']:1})

        atom_counts = collections.Counter()
        for cpd_id, cpd_dict in cpds.items():
            for atom_id, atom_count in cpd_dict['atom_count'].items():
                atom_counts[atom_id] += atom_count*cpd_counter[cpd_id]

        cpd_returns = [(stoich, cpds[cpd_id]) for
                       cpd_id, stoich in cpd_counter.items()]

        return cpd_returns, atom_counts

    def _gen_compound(mol):
        try:
            if explicit_h:
                mol = RemoveHs(mol)
            SanitizeMol(mol)
        except:
            return None

        mol_smiles = MolToSmiles(mol, True)
        if '.' in mol_smiles:
            return None

        cpd_id = utils.compound_hash(mol_smiles, 'Predicted')
        if cpd_id:
            if cpd_id not in local_cpds:
                cpd_dict = {'ID': None, '_id': cpd_id, 'SMILES': mol_smiles,
                                'Type': 'Predicted',
                                'Generation': generation,
                                'atom_count': utils._getatom_count(mol),
                                'Reactant_in': [], 'Product_of': [],
                                'Expand': True,
                                'Formula': CalcMolFormula(mol),
                                'last_tani': 0}
            else:
                cpd_dict = local_cpds[cpd_id]

            return cpd_dict
        else:
            return None

    try:
        product_sets = rule[0].RunReactants(reactant_mols)
        reactants, reactant_atoms = _make_half_rxn(reactant_mols,
                                                   rule[1]['Reactants'])
    except:
        reactants = None

    if not reactants:
        return local_cpds, local_rxns

    for product_mols in product_sets:
        try:
            products, product_atoms = _make_half_rxn(product_mols,
                                                     rule[1]['Products'])
            if not products:
                continue

            if (
                reactant_atoms - product_atoms or
                product_atoms - reactant_atoms
               ):
                is_atom_balanced = False
            else:
                is_atom_balanced = True

            if is_atom_balanced:
                for _, cpd_dict in products:
                    if cpd_dict['_id'].startswith('C'):
                        local_cpds.update({cpd_dict['_id']: cpd_dict})

                rhash, rxn_text = utils.rxn2hash(reactants, products)
                if rhash not in local_rxns:
                    local_rxns[rhash] = {
                            '_id': rhash,
                            # give stoich and id of reactants/products
                            'Reactants': [(s, r['_id']) for s, r in reactants],
                            'Products': [(s, p['_id']) for s, p in products],
                            'Operators': {rule_name},
                            'SMILES_rxn': rxn_text
                            }
                else:
                    local_rxns[rhash]['Operators'].add(rule_name)

        except (ValueError, MemoryError) as e:
            continue
    # return compounds and reactions to be added into the local
    return local_cpds, local_rxns


# Full Operators
def _transform_ind_compound_with_full(coreactant_mols: dict, coreactant_dict: dict,
                                      operators: list, generation: int, explicit_h: bool,
                                      compound_smiles: list):
    local_cpds = dict()
    local_rxns = dict()

    mol = MolFromSmiles(compound_smiles)
    mol = RemoveHs(mol)
    if not mol:
        print(f"Unable to parse: {compound_smiles}")
        return None
    Kekulize(mol, clearAromaticFlags=True)
    if explicit_h:
        mol = AddHs(mol)
    # Apply reaction rules to prepared compound

    # run through the single compound operatores
    for rule_name, rule in operators.items():
        # Get RDKit Mol objects for reactants
        reactant_mols = tuple([mol if x == 'Any'
                               else coreactant_mols[x][0]
                               for x in rule[1]['Reactants']])
        # Perform chemical reaction on reactants for each rule
        # try:
        generated_cpds, generated_rxns = _run_reaction(
            rule_name, rule, reactant_mols, coreactant_mols,
            coreactant_dict, local_cpds, local_rxns, generation, explicit_h
        )

        local_cpds.update(generated_cpds)
        for rxn, vals in generated_rxns.items():
            if rxn in local_rxns:
                local_rxns[rxn]['Operators'].union(vals['Operators'])

    return local_cpds, local_rxns


def transform_all_compounds_with_full(compound_smiles: list, coreactants: dict,
                                       coreactant_dict: dict, operators: dict, generation: int,
                                       explicit_h: bool, processes: int) -> Tuple[dict, dict]:
    """Transform compounds given a list of rules

    Carry out the transformation of a list of compounds given operators. 
    Generates new products and returns them to be processed by pickaxe.

    Parameters
    ----------
    compound_smiles : list
        List of SMILES to react
    coreactants : dict
        Dictionary of correactants RDKit Mols defined in rules
    coreactant_dict : dict
        Dictionary of correactant compoudnds defined in rules
    operators : dict
        Dictionary of reaction rules
    generation : int
        Value of generation to expand
    explicit_h : bool
        Whether or not to have explicit Hs in reactions
    processes : int
        Number of processors being used

    Returns
    -------
    Tuple[dict, dict]
        Returns a tuple of New Compounds and New Reactants
    """
    def print_progress(done, total):
        # Use print_on to print % completion roughly every 2.5 percent
        # Include max to print no more than once per compound (e.g. if
        # less than 20 compounds)
        print_on = max(round(0.1 * total), 1)
        if not done % print_on:
            print(f"Generation {generation}: {round(done / total * 100)}"
                  " percent complete")

    # First transform
    new_cpds_master = {}
    new_rxns_master = {}

    transform_compound_partial = partial(
        _transform_ind_compound_with_full,
        coreactants,
        coreactant_dict,
        operators,
        generation,
        explicit_h
    )
    # par loop
    if processes > 1:
        # TODO chunk size?
        # chunk_size = max(
        #         [round(len(compound_smiles) / (processes)), 1])
        chunk_size = 1
        # print(f'Chunk size = {chunk_size}')
        pool = multiprocessing.Pool(processes=processes)
        for i, res in enumerate(pool.imap_unordered(
                            transform_compound_partial,
                            compound_smiles,
                            chunk_size)):
            new_cpds, new_rxns = res
            new_cpds_master.update(new_cpds)

            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    ops_set = rxn_dict["Operators"]
                    new_rxns_master[rxn]['Operators'].union(ops_set)
                else:
                    new_rxns_master.update({rxn: rxn_dict})
            print_progress(i, len(compound_smiles))

    else:
        for i, smiles in enumerate(compound_smiles):
            new_cpds, new_rxns = transform_compound_partial(smiles)
            # new_cpds as cpd_id:cpd_dict
            # new_rxns as rxn_id:rxn_dict
            new_cpds_master.update(new_cpds)
            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    ops_set = rxn_dict["Operators"]
                    new_rxns_master[rxn]['Operators'].union(ops_set)
                else:
                    new_rxns_master.update({rxn: rxn_dict})
            print_progress(i, len(compound_smiles))

    return new_cpds_master, new_rxns_master

# Partial Operators
def _transform_ind_compound_with_partial(coreactant_mols: dict, coreactant_dict: dict,
                                         operators: dict, generation: int, explicit_h: bool,
                                         partial_rules: dict, compound_smiles: list):
    # 1. See if rule matches the compound passed
    #   (rule from partial_rules dict keys)
    # 2. If match apply transform_ind_compound_with_full to each
    def generate_partial_mols(partial_rule):
        def gen_mol(smi):
            mol = MolFromSmiles(smi)
            mol = RemoveHs(mol)
            Kekulize(mol, clearAromaticFlags=True)
            if explicit_h:
                mol = AddHs(mol)
            return mol

        rule_reactants = operators[partial_rule['rule']][1]['Reactants']
        cofactor = [False if r == 'Any' else True for r in rule_reactants]
        reactant_mols = []
        for is_cofactor, smi in zip(cofactor, partial_rule['reactants']):
                if is_cofactor:
                    reactant_mols.append(coreactant_mols[smi][0])
                elif smi == 'SMARTS_match':
                    reactant_mols.append(gen_mol(compound_smiles))
                else:
                    # These reactions already happen with any;any
                    if (utils.compound_hash(smi) !=
                            utils.compound_hash(compound_smiles)):
                        reactant_mols.append(gen_mol(smi))
                    else:
                        return None
        return reactant_mols

    local_cpds = dict()
    local_rxns = dict()

    mol = MolFromSmiles(compound_smiles)
    mol = RemoveHs(mol)
    if not mol:
        print(f"Unable to parse: {compound_smiles}")
        return None
    Kekulize(mol, clearAromaticFlags=True)
    if explicit_h:
        mol = AddHs(mol)
    # Apply reaction rules to prepared compound

    # run through the single compound operatores
    for ind_SMARTS, rules in partial_rules.items():
        # does mol match vs smiles match change things?
        if QuickSmartsMatch(compound_smiles, ind_SMARTS):
            for partial_rule in rules:
                # Perform chemical reaction on reactants for each rule
                # try:
                rule_name = partial_rule['rule_reaction'].split('_')[0]
                rule = operators[partial_rule['rule']]
                reactant_mols = generate_partial_mols(partial_rule)
                if reactant_mols:
                    generated_cpds, generated_rxns = _run_reaction(rule_name, rule,
                        reactant_mols, coreactant_mols, coreactant_dict, local_cpds, local_rxns, generation, explicit_h)


                    local_cpds.update(generated_cpds)
                    for rxn, vals in generated_rxns.items():
                        if rxn in local_rxns:
                            if 'Partial Operators' in local_rxns[rxn]:
                                local_rxns[rxn]['Partial Operators'].update([partial_rule['rule_reaction']])
                            else:
                                local_rxns[rxn]['Partial Operators'] = set([partial_rule['rule_reaction']])
    return local_cpds,local_rxns


def _transform_all_compounds_with_partial(compound_smiles, coreactants,
                                          coreactant_dict, operators,
                                          generation, explicit_h, processes,
                                          partial_rules):
    """
    This function is made to reduce the memory load of parallelization.
    This function accepts in a list of cpds (cpd_list) and runs the transformation in parallel of these.
    """
    def print_progress(done, total):
            # Use print_on to print % completion roughly every 2.5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(0.1 * total), 1)
            if not done % print_on:
                print(f"Generation {generation}: {round(done / total * 100)} percent complete")

    # First transform
    new_cpds_master = {}
    new_rxns_master = {}

    transform_compound_partial = partial(_transform_ind_compound_with_partial, coreactants,
                                            coreactant_dict, operators, generation, explicit_h, partial_rules)
    # par loop
    if processes > 1:
        # TODO chunk size?
        # chunk_size = max(
        #         [round(len(compound_smiles) / (processes)), 1])
        chunk_size = 1
        # print(f'Chunk size = {chunk_size}')
        pool = multiprocessing.Pool(processes=processes)
        for i, res in enumerate(pool.imap_unordered(
                            transform_compound_partial, compound_smiles, chunk_size)):
            new_cpds, new_rxns = res
            new_cpds_master.update(new_cpds)

            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    new_rxns_master[rxn]['Operators'].union(rxn_dict['Operators'])
                else:
                    new_rxns_master.update({rxn:rxn_dict})
            print_progress(i, len(compound_smiles))

    else:
        for i, smiles in enumerate(compound_smiles):
            new_cpds, new_rxns = transform_compound_partial(smiles)
            # new_cpds as cpd_id:cpd_dict
            # new_rxns as rxn_id:rxn_dict
            new_cpds_master.update(new_cpds)
            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    new_rxns_master[rxn]['Partial Operators'] = new_rxns_master[rxn]['Partial Operators'].union(rxn_dict['Partial Operators'])
                else:
                    new_rxns_master.update({rxn:rxn_dict})
            print_progress(i, len(compound_smiles))

    return new_cpds_master, new_rxns_master

