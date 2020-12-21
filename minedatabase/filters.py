from functools import partial
import multiprocessing

from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem import rdFMCS as mcs
from rdkit.Chem import AllChem

import pandas as pd
import numpy as np
from scipy.stats import rv_discrete

from minedatabase.utils import get_fp

###############################################################################
# Tanimoto Sampling Filter


def _filter_by_sample(self, weight=None, num_workers=1):
    """
    Samples N compounds to expand based on the weighted Tanimoto
    distribution.
    """
    print(f'Filtering Generation {self.generation} via Tanimoto Sampling.')

    # Get compounds eligible for expansion in the current generation
    compounds_to_check = []

    for cpd in self.compounds.values():
        # Compounds are in generation and correct type
        if (cpd['Generation'] == self.generation
                and cpd['Type'] not in ['Coreactant', 'Target Compound']):

            # Check for targets and only react if terminal
            if self.react_targets:
                compounds_to_check.append(cpd)
            else:
                for t_id in self.targets:
                    if 'C' + t_id[1:] != cpd['_id']:
                        compounds_to_check.append(cpd)
                    else:
                        self.compounds[cpd['_id']]['Expand'] = False

    # Get compounds to keep
    cpd_info = [(cpd['_id'], cpd['SMILES']) for cpd in compounds_to_check]
    sampled_ids = sample_by_tanimoto(cpd_info, self.target_fps,
                                     self.sample_size, min_T=0.05,
                                     weighting=self.sample_weight,
                                     max_iter=10**7, n_cores=num_workers)
    # Get compounds to remove
    ids = set(i[0] for i in cpd_info)
    cpds_remove_set = ids - sampled_ids

    for c_id in cpds_remove_set:
        self.compounds[c_id]['Expand'] = False

    self._apply_filter_results(compounds_to_check)

    return None


def sample_by_tanimoto(mol_info, t_fp, n_cpds=None, min_T=0.05,
                       weighting=None, max_iter=None,
                       n_cores=1):
    """
    Given a list of ids and compounds, this function uses inverse-CDF
    sampling from a PMF generated from weighted tanimoto similarity
    to select compounds to expand.

    :param mol_info:  A list consisting of (cpd_id, SMILES)
    :type mol_info: list(tuple)
    :param t_fp: A list consistent of the target fingerprints
    :type t_fp: list(rdkit.DataStructs.cDataStructs.ExplicitBitVect)

    #TODO

    :return: A list of cpd_id to be expanded
    :rtype: list(str)
    """
    if len(mol_info) <= n_cpds:
        ids = set(x[0] for x in mol_info)
        return ids

    # Get random, discrete distribution and associated ids
    rv, ids = _gen_rv_from_tanimoto(mol_info, t_fp, weighting, min_T, n_cores)

    # Sample intervals from rv and get c_id from id
    if max_iter is None:
        max_iter = n_cpds*5

    chosen = set()
    i = 0

    while len(chosen) < n_cpds and i < max_iter:
        chosen.add(rv.rvs(size=1)[0])
        i += 1

    # Map sampled intervals into pickaxe c_ids
    chosen_ids = set(ids.iloc[chosen_id] for chosen_id in chosen)

    return chosen_ids


def _gen_rv_from_tanimoto(mol_info, t_fp, weighting=None, min_T=0.05,
                          n_cores=1):
    """Generate a CDF in a pandas dataframe from a list of IDs and smiles"""
    if weighting is None:
        def weight_f(T):
            return T**4

    # Construct target df
    t_df = pd.DataFrame(t_fp, columns=['fp'])

    # Construct targets to sample df
    df = pd.DataFrame(mol_info, columns=['_id', 'SMILES'])

    # Calculate Tanimoto for each compound and drop T < min_T
    partial_T_calc = partial(_calc_max_T, t_df, min_T)
    df = _parallelize_dataframe(df, partial_T_calc, n_cores)

    # Generate CDF for sampling
    df['T_trans'] = df['T'].map(weight_f)
    df['T_pdf'] = df['T_trans']/sum(df['T_trans'])

    # Generate CDF
    df.reset_index(inplace=True, drop=True)
    xk = df.index
    pk = df['T_pdf']
    rv = rv_discrete(values=(xk, pk))

    return rv, df['_id']


def _parallelize_dataframe(df, func, n_cores=1):
    """
    Applies a function to a dataframe in parallel by chunking it up over
    the specified number of cores.
    """
    df_split = np.array_split(df, n_cores)
    pool = multiprocessing.Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def _calc_max_T(t_df, min_T, df):
    """
    Generate the tanimoto to use to generate the PMF to sample from.
    For each compound a list of tanimoito values are obtained by a generated
    compound to every target compound and the max is taken.
    """
    df['fp'] = df['SMILES'].map(get_fp)

    df['T'] = None
    for i in range(len(df)):
        fp = df['fp'].iloc[i]
        df['T'].iloc[i] = max(t_df['fp'].map(
            lambda x: FingerprintSimilarity(x, fp))
        )
    # Filter out low Tanimoto
    df = df[df['T'] > min_T]

    return df

# End Tanimoto Sampling Filter
###############################################################################

###############################################################################
# Hard cutoff filters -- e.g. Tanimoto and MCS metric


def _filter_by_tani(self, num_workers=1):
    """
    Compares the current generation to the target compound fingerprints
    marking compounds, who have a tanimoto similarity score to a target
    compound greater than or equal to the crit_tani, for expansion.
    """

    # Set up variables required for filtering
    # Tanimoto Threshold
    if type(self.crit_tani) == list:
        crit_tani = self.crit_tani[self.generation]
    else:
        crit_tani = self.crit_tani

    # Get compounds eligible for expansion in the current generation
    compounds_to_check = []

    for cpd in self.compounds.values():
        # Compounds are in generation and correct type
        if (cpd['Generation'] == self.generation
                and cpd['Type'] not in ['Coreactant', 'Target Compound']):

            # Check for targets and only react if terminal
            if self.react_targets:
                compounds_to_check.append(cpd)
            else:
                for t_id in self.targets:
                    if 'C' + t_id[1:] != cpd['_id']:
                        compounds_to_check.append(cpd)
                    else:
                        self.compounds[cpd['_id']]['Expand'] = False

    # Run the filtering code to get a list of compounds to ignore
    print(f"Filtering Generation {self.generation} "
          f"with Tanimoto > {crit_tani}.")
    # Get input to filter code, c_id and smiles (to be turned into fingerprint)
    cpd_info = [(cpd['_id'], cpd['SMILES']) for cpd in compounds_to_check]
    cpd_filters = _filter_by_tani_helper(
                        cpd_info,
                        self.target_fps,
                        crit_tani, num_workers
                    )

    # Process filtering results
    for c_id, current_tani in cpd_filters:
        if current_tani == -1:
            self.compounds[c_id]['Expand'] = False
        else:
            # Check if tani is increasing
            if self.increasing_tani:
                if current_tani >= self.compounds[c_id]['last_tani']:
                    self.compounds[c_id]['last_tani'] = current_tani
                else:
                    # tanimoto isn't increasing
                    self.compounds[c_id]['Expand'] = False

    self._apply_filter_results(compounds_to_check)


def _filter_by_tani_helper(compounds_info, target_fps, crit_tani,
                           num_workers, retro=False):
    def print_progress(done, total, section):
        # Use print_on to print % completion roughly every 5 percent
        # Include max to print no more than once per compound (e.g. if
        # less than 20 compounds)
        print_on = max(round(.05 * total), 1)
        if not (done % print_on):
            print(f"{section} {round(done / total * 100)} percent complete")

    # compound_info = [(smiles, id)]
    cpds_to_filter = list()
    compare_target_fps_partial = partial(_compare_target_fps, target_fps,
                                         crit_tani)

    if num_workers > 1:
        # Set up parallel computing of compounds to
        chunk_size = max(
                    [round(len(compounds_info) / (num_workers * 4)), 1])
        pool = multiprocessing.Pool(num_workers)
        for i, res in enumerate(pool.imap_unordered(
                compare_target_fps_partial, compounds_info, chunk_size)):
            # If the result of comparison is false, compound is not expanded
            # Default value for a compound is True, so no need to
            # specify expansion
            if res:
                cpds_to_filter.append(res)
            print_progress(i, len(compounds_info), 'Tanimoto filter progress:')

    else:
        for i, cpd in enumerate(compounds_info):
            res = compare_target_fps_partial(cpd)
            if res:
                cpds_to_filter.append(res)
            print_progress(i, len(compounds_info), 'Tanimoto filter progress:')
    print("Tanimoto filter progress: 100 percent complete")
    return cpds_to_filter


def _compare_target_fps(target_fps, crit_tani, compound_info):
    # do finger print loop here
    """
    Helper function to allow parallel computation of tanimoto filtering.
    Works with _filter_by_tani_helper

    Returns cpd_id if a the compound is similar enough to a target.

    """
    # Generate the fingerprint of a compound and compare to the fingerprints
    # of the targets
    try:
        fp1 = get_fp(compound_info[1])
        for fp2 in target_fps:
            tani = AllChem.DataStructs.FingerprintSimilarity(fp1, fp2)
            if tani >= crit_tani:
                return (compound_info[0], tani)
    except:
        return (compound_info[0], -1)


def _filter_by_MCS(self, num_workers=1):
    """
    Compares the current generation to the target compound fingerprints
    marking compounds, who have a tanimoto similarity score to a target
    compound greater than or equal to the crit_tani, for expansion.
    """

    # Set up variables required for filtering
    # Tanimoto Threshold
    if type(self.crit_mcs) == list:
        crit_mcs = self.crit_mcs[self.generation]
    else:
        crit_mcs = self.crit_mcs

    # Get compounds eligible for expansion in the current generation
    compounds_to_check = []

    for cpd in self.compounds.values():
        # Compounds are in generation and correct type
        if (cpd['Generation'] == self.generation
                and cpd['Type'] not in ['Coreactant', 'Target Compound']):

            # Check for targets and only react if terminal
            if self.react_targets:
                compounds_to_check.append(cpd)
            else:
                for t_id in self.targets:
                    if 'C' + t_id[1:] != cpd['_id']:
                        compounds_to_check.append(cpd)
                    else:
                        self.compounds[cpd['_id']]['Expand'] = False
                
    # Run the filtering code to get a list of compounds to ignore
    print(f"Filtering Generation {self.generation} "
          f"with MCS > {crit_mcs}.")
    # Get input to filter code, c_id and smiles (to be turned into fingerprint)
    cpd_info = [(cpd['_id'], cpd['SMILES']) for cpd in compounds_to_check]
    cpd_filters = _filter_by_mcs_helper(
                        cpd_info,
                        self.target_smiles,
                        crit_mcs, num_workers
                    )

    # Process filtering results
    for c_id, current_tani in cpd_filters:
        if current_tani == -1:
            self.compounds[c_id]['Expand'] = False
        else:
            # Check if tani is increasing
            if self.increasing_tani:
                if current_tani >= self.compounds[c_id]['last_tani']:
                    self.compounds[c_id]['last_tani'] = current_tani
                else:
                    # tanimoto isn't increasing
                    self.compounds[c_id]['Expand'] = False

    self._apply_filter_results(compounds_to_check)


def _filter_by_mcs_helper(compounds_info, target_smiles, crit_mcs, num_workers,
                          retro=False):
    def print_progress(done, total, section):
        # Use print_on to print % completion roughly every 5 percent
        # Include max to print no more than once per compound (e.g. if
        # less than 20 compounds)
        print_on = max(round(.05 * total), 1)
        if not (done % print_on):
            print(f"{section} {round(done / total * 100)} percent complete")

    # compound_info = [(smiles, id)]
    cpds_to_filter = list()
    compare_target_mcs_partial = partial(_compare_target_mcs,
                                    target_smiles, crit_mcs, retro)

    if num_workers > 1:
        # Set up parallel computing of compounds to
        chunk_size = max(
                    [round(len(compounds_info) / (num_workers * 4)), 1])
        pool = multiprocessing.Pool(num_workers)
        for i, res in enumerate(pool.imap_unordered(
                compare_target_mcs_partial, compounds_info, chunk_size)):

            if res:
                cpds_to_filter.append(res)
            print_progress(i, len(compounds_info),
                            'Maximum Common Substructure filter progress:')

    else:
        for i, cpd in enumerate(compounds_info):
            res = compare_target_mcs_partial(cpd)
            if res:
                cpds_to_filter.append(res)
            print_progress(i, len(compounds_info),
                           'Maximum Common Substructure filter progress:')

    print("Maximum Common Substructure filter progress: 100 percent complete")
    return cpds_to_filter


def _compare_target_mcs(target_smiles, crit_mcs, retro, compound_info):
    """
    Helper function to allow parallel computation of MCS filtering.
    Works with _filter_by_tani_helper

    Returns cpd_id if a the compound is similar enough to a target.

    """
    def get_mcs_overlap(mol, target_mol):
        mcs_out = mcs.FindMCS([mol, target_mol],
                              matchValences=False,
                              ringMatchesRingOnly=False)

        if not mcs_out.canceled:
            ss_atoms = mcs_out.numAtoms
            ss_bonds = mcs_out.numBonds
            t_atoms = target_mol.GetNumAtoms()
            t_bonds = target_mol.GetNumBonds()

            mcs_overlap = ((ss_atoms + ss_bonds) / (t_bonds + t_atoms))
            return mcs_overlap

        else:
            return 0
    # compare MCS for filter
    try:
        mol = AllChem.MolFromSmiles(compound_info[1])

        for t_smi in target_smiles:
            t_mol = AllChem.MolFromSmiles(t_smi)
            if not retro:
                mcs_overlap = get_mcs_overlap(mol, t_mol)
            else:
                mcs_overlap = get_mcs_overlap(t_mol, mol)

            if mcs_overlap > 1:
                print("pause")
            if mcs_overlap >= crit_mcs:
                return (compound_info[0], mcs_overlap)
    except:
        return (compound_info[0], -1)
#
#
# End Hard cutoff filters
###############################################################################

###############################################################################
# Agnostic Filtering Functions


def _apply_filter_results(self, compounds_to_check):
    """
    Remove compounds and reactions that can be removed
    For a compound to be removed it must:
        1. Not be flagged for expansion
        2. Not have a coproduct in a reaction marked for expansion
        3. Start with 'C'
    """
    def should_delete_reaction(rxn_id):
        products = self.reactions[rxn_id]['Products']
        for _, c_id in products:
            if c_id.startswith('C') and c_id not in cpds_to_remove:
                return False
        # Every compound isn't in cpds_to_remove
        return True

    cpds_to_remove = set()
    rxns_to_check = set()
    for cpd_dict in compounds_to_check:
        cpd_id = cpd_dict['_id']
        if not cpd_dict['Expand'] and cpd_id.startswith('C'):
            cpds_to_remove.add(cpd_id)
            # Generate set of reactions to remove
            rxn_ids = set(self.compounds[cpd_id]['Product_of']
                          + self.compounds[cpd_id]['Reactant_in'])

            rxns_to_check = rxns_to_check.union(rxn_ids)

    # Function to check to see if should delete reaction
    # If reaction has compound that won't be deleted keep it
    # Check reactions for deletion
    for rxn_id in rxns_to_check:
        if should_delete_reaction(rxn_id):
            for _, c_id in self.reactions[rxn_id]['Products']:
                if c_id.startswith('C'):
                    if rxn_id in self.compounds[c_id]['Product_of']:
                        self.compounds[c_id]['Product_of'].remove(rxn_id)

            for _, c_id in self.reactions[rxn_id]['Reactants']:
                if c_id.startswith('C'):
                    if rxn_id in self.compounds[c_id]['Reactant_in']:
                        self.compounds[c_id]['Reactant_in'].remove(rxn_id)

            del(self.reactions[rxn_id])
        else:
            # Reaction is dependent on compound that is flagged to be
            # removed. Don't remove compound
            products = self.reactions[rxn_id]['Products']
            cpds_to_remove -= set(i[1] for i in products)

            # for _, c_id in products:
            #     if c_id in cpds_to_remove:
            #         cpds_to_remove -= {c_id}

    # Remove compounds and reactions if any found
    for cpd_id in cpds_to_remove:
        del(self.compounds[cpd_id])

# End of Agnostic Filtering Functions
###############################################################################
