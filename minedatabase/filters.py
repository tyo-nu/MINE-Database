"""Definitions of filters for pickaxe.

Use this module to define your own filter classes. All filter classes must
subclass Filter. See Filter docstring for more information.

To use any filter, import it in pickaxe_run.py, initialize it, and append it
to the .filters property of your pickaxe object.
"""

import abc
import copy
import multiprocessing
import time
from functools import partial

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem, CanonSmiles
from rdkit.Chem import rdFMCS as mcs
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.DataStructs import FingerprintSimilarity
from scipy.stats import rv_discrete

from minedatabase import utils
from minedatabase.metabolomics import MetabolomicsDataset, Peak
from minedatabase.utils import get_fp


###############################################################################
# ABC for all Filter Subclasses


class Filter(metaclass=abc.ABCMeta):
    """ABC for all Filter subclasses.

    All subclasses must implement properties and methods decorated with
    @abc.abstractmethod. Feel free to override other non-private methods as
    well, such as pre_print() and post_print().
    """

    @property
    @abc.abstractmethod
    def filter_name(self):
        """Obtain name of filter."""
        pass

    @abc.abstractmethod
    def _choose_cpds_to_filter(self, pickaxe, num_workers):
        """Return list of compounds to remove from pickaxe object."""
        pass

    def apply_filter(self, pickaxe, num_workers=1, print_on=True):
        """Apply filter from Pickaxe object."""
        time_sample = time.time()

        if print_on:
            n_total = self._get_n(pickaxe, "total")
            self.pre_print_header(pickaxe)
            self.pre_print()

        compound_ids_to_check = self._choose_cpds_to_filter(pickaxe, num_workers)

        if compound_ids_to_check:
            self._apply_filter_results(pickaxe, compound_ids_to_check)

        if print_on:
            n_filtered = self._get_n(pickaxe, "filtered")
            self.post_print(pickaxe, n_total, n_filtered, time_sample)
            self.post_print_footer(pickaxe)

    def pre_print_header(self, pickaxe):
        """Print header before filtering."""
        print("----------------------------------------")
        print(f"Filtering Generation {pickaxe.generation}\n")

    def pre_print(self):
        """Print filter being applied."""
        print(f"Applying filter: {self.filter_name}")

    def post_print(self, pickaxe, n_total, n_filtered, time_sample):
        """Print results of filtering."""
        print(
            f"{n_filtered} of {n_total} compounds remain after applying "
            f"filter: {self.filter_name}"
            f"--took {round(time.time() - time_sample, 2)}s.\n"
        )

    def post_print_footer(self, pickaxe):
        """Print end of filtering."""
        print(f"Done filtering Generation {pickaxe.generation}")
        print("----------------------------------------\n")

    def _get_n(self, pickaxe, n_type):
        """Get current number of compounds to be filtered."""
        n = 0
        for cpd_dict in pickaxe.compounds.values():
            is_in_current_gen = cpd_dict["Generation"] == pickaxe.generation
            is_predicted_compound = cpd_dict["_id"].startswith("C")
            if is_in_current_gen and is_predicted_compound:
                if n_type == "total":
                    n += 1
                elif n_type == "filtered" and cpd_dict["Expand"]:
                    n += 1
        return n

    def _apply_filter_results(self, pickaxe, compound_ids_to_check):
        """Apply filter results to Pickaxe object.

        Remove compounds and reactions that can be removed
        For a compound to be removed it must:
            1. Not be flagged for expansion
            2. Not have a coproduct in a reaction marked for expansion
            3. Start with 'C'
        """

        def should_delete_reaction(rxn_id):
            products = pickaxe.reactions[rxn_id]["Products"]
            for _, c_id in products:
                if c_id.startswith("C") and c_id not in cpds_to_remove:
                    return False
            # Every compound isn't in cpds_to_remove
            return True

        def get_compounds_to_check_from_ids(pickaxe, cpd_ids_to_check):
            cpds_to_check = []
            for cpd in pickaxe.compounds.values():
                if cpd["_id"] in cpd_ids_to_check:
                    cpds_to_check.append(cpd)
            return cpds_to_check

        compounds_to_check = get_compounds_to_check_from_ids(
            pickaxe, compound_ids_to_check
        )

        cpds_to_remove = set()
        rxns_to_check = set()
        for cpd_dict in compounds_to_check:
            cpd_id = cpd_dict["_id"]
            if not cpd_dict["Expand"] and cpd_id.startswith("C"):
                cpds_to_remove.add(cpd_id)
                # Generate set of reactions to remove
                rxn_ids = set(
                    pickaxe.compounds[cpd_id]["Product_of"]
                    + pickaxe.compounds[cpd_id]["Reactant_in"]
                )

                rxns_to_check = rxns_to_check.union(rxn_ids)

        # Function to check to see if should delete reaction
        # If reaction has compound that won't be deleted keep it
        # Check reactions for deletion
        for rxn_id in rxns_to_check:
            if should_delete_reaction(rxn_id):
                for _, c_id in pickaxe.reactions[rxn_id]["Products"]:
                    if c_id.startswith("C"):
                        if rxn_id in pickaxe.compounds[c_id]["Product_of"]:
                            pickaxe.compounds[c_id]["Product_of"].remove(rxn_id)

                for _, c_id in pickaxe.reactions[rxn_id]["Reactants"]:
                    if c_id.startswith("C"):
                        if rxn_id in pickaxe.compounds[c_id]["Reactant_in"]:
                            pickaxe.compounds[c_id]["Reactant_in"].remove(rxn_id)

                del pickaxe.reactions[rxn_id]
            else:
                # Reaction is dependent on compound that is flagged to be
                # removed. Don't remove compound
                products = pickaxe.reactions[rxn_id]["Products"]
                cpds_to_remove -= set(i[1] for i in products)

                # for _, c_id in products:
                #     if c_id in cpds_to_remove:
                #         cpds_to_remove -= {c_id}

        # Remove compounds and reactions if any found
        for cpd_id in cpds_to_remove:
            del pickaxe.compounds[cpd_id]


###############################################################################
# Tanimoto Sampling Filter


class TanimotoSamplingFilter(Filter):
    """Filter that samples randomly from weighted tanimoto."""

    def __init__(self, filter_name: str, sample_size: int, weight=None):
        """Initialize filtering class.

        Parameters
        ----------
        filter_name : [type]
            [description]
        sample_size : [type]
            [description]
        weight : [type]
            [description]
        """
        self._filter_name = filter_name
        self.sample_size = sample_size
        self.sample_weight = weight

    @property
    def filter_name(self):
        return self._filter_name

    def pre_print(self):
        print(
            (
                f"Sampling {self.sample_size} Compounds Based on a "
                f"Weighted Tanimoto Distribution"
            )
        )

    def post_print(self, pickaxe, n_total, n_filtered, time_sample):
        print(
            (
                f"{n_filtered} of {n_total} "
                "compounds selected after "
                f"Tanimoto Sampling of generation {pickaxe.generation}"
                f"--took {time.time() - time_sample}s.\n"
            )
        )

    def _choose_cpds_to_filter(self, pickaxe, num_workers):
        """
        Samples N compounds to expand based on the weighted Tanimoto
        distribution.
        """
        print(f"Filtering Generation {pickaxe.generation}" " via Tanimoto Sampling.")

        if not pickaxe.target_fps:
            print("No targets to filter for. Can't expand.")
            return None

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = []
        set_unreactive = True

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:

                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    compounds_to_check.append(cpd)
                else:
                    for t_id in pickaxe.targets:
                        if "C" + t_id[1:] != cpd["_id"]:
                            compounds_to_check.append(cpd)
                            set_unreactive = False
                            break

                    if set_unreactive:
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                    else:
                        set_unreactive = True

        # Get compounds to keep
        cpd_info = [(cpd["_id"], cpd["SMILES"]) for cpd in compounds_to_check]

        sampled_ids = self.sample_by_tanimoto(
            cpd_info,
            pickaxe.target_fps,
            self.sample_size,
            min_T=0.15,
            weighting=self.sample_weight,
            max_iter=None,
            n_cores=num_workers,
        )

        # Get compounds to remove
        ids = set(i[0] for i in cpd_info)
        cpds_remove_set = ids - sampled_ids

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set

    def sample_by_tanimoto(
        self,
        mol_info,
        t_fp,
        n_cpds=None,
        min_T=0.05,
        weighting=None,
        max_iter=None,
        n_cores=1,
    ):
        """
        Given a list of ids and compounds, this function uses inverse-CDF
        sampling from a PMF generated from weighted tanimoto similarity
        to select compounds to expand.

        :param mol_info:  A list consisting of (cpd_id, SMILES)
        :type mol_info: list(tuple)
        :param t_fp: A list of the target fingerprints
        :type t_fp: list(rdkit.DataStructs.cDataStructs.ExplicitBitVect)
        :param n_cpds: Number of compounds to sample
        :type n_cpds: int
        :param min_T: Minimum Tanimoto for consideration
        :type min_T: float
        :param weighting: Weighting function that takes Tanimoto as input
        :type weighting: func
        :param max_iter: Maximum iterations before new CDF is calculated
        :type max_iter: int
        :param n_cores: Number of CPU cores to use
        :type n_cores: int

        :return: A set of cpd_ids to be expanded
        :rtype: set(str)
        """

        # Return input if less than desired number of compounds
        if len(mol_info) <= n_cpds:
            ids = set(x[0] for x in mol_info)
            print(
                "-- Number to sample is less than number of compounds. "
                "Returning all compounds."
            )
            return ids

        # Get pandas df and ids
        df = self._gen_df_from_tanimoto(mol_info, t_fp, min_T=min_T, n_cores=n_cores)

        if len(df) <= n_cpds:
            ids = set(df["_id"])
            print(
                f"-- After filtering by minimum tanimoto ({min_T}) "
                "number to sample is less than number of compounds. "
                "Returning all compounds."
            )
            return ids

        print("-- Sampling compounds to expand.")
        then = time.time()
        # Get discrete distribution to sample randomly from
        rv, ids = self._gen_rv_from_df(df, weighting=weighting)

        # Sample intervals from rv and get c_id from id
        if max_iter is None:
            max_iter = n_cpds / 10 if n_cpds > 1000 else n_cpds / 2

        chosen_ids = set()
        i = 0
        nCDF = 0

        while len(chosen_ids) != n_cpds:
            # if current iteration if greater than max then
            # recalc distribution to exclude chosen
            if i > max_iter:
                i = 0
                nCDF += 1
                rv, ids = self._gen_rv_from_df(
                    df, chosen=chosen_ids, weighting=weighting
                )

            chosen_ids.add(ids.iloc[rv.rvs(size=1)[0]])
            i += 1

        print(
            f"-- Finished sampling in {time.time() - then} s."
            f" Recalculated CDF {nCDF} times."
        )

        return chosen_ids

    def _gen_rv_from_df(self, df, chosen=[], weighting=None):
        """Genderate a scipy.rv object to sample from."""
        if weighting is None:

            def weighting(T):
                return T ** 4

        # TODO Make more memory efficient... maybe use np directly instead?
        # Could be due to spawn vs fork
        rescale_df = copy.copy(df[~df["_id"].isin(chosen)])
        rescale_df.loc[:, "T_trans"] = rescale_df["T"].map(weighting)
        rescale_df.loc[:, "T_pdf"] = rescale_df["T_trans"] / sum(rescale_df["T_trans"])

        # Generate CDF
        rescale_df.reset_index(inplace=True, drop=True)
        xk = rescale_df.index
        pk = rescale_df["T_pdf"]
        rv = rv_discrete(values=(xk, pk))
        ids = rescale_df["_id"]

        del rescale_df

        return rv, ids

    def _gen_df_from_tanimoto(self, mol_info, t_fp, min_T, n_cores):
        # Construct target df

        def chunks(lst, size):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), size):
                yield lst[i : i + size]

        then = time.time()
        print("-- Calculating Fingerprints and Tanimoto Values.")
        # target fingerprint dataframe
        t_df = pd.DataFrame(t_fp, columns=["fp"])

        # Calculate Tanimoto for each compound and drop T < min_T
        partial_T_calc = partial(_calc_max_T, t_df, min_T)

        df = pd.DataFrame()
        for mol_chunk in chunks(mol_info, 10000):

            # Construct targets to sample df
            temp_df = pd.DataFrame(mol_chunk, columns=["_id", "SMILES"])
            df = df.append(_parallelize_dataframe(temp_df, partial_T_calc, n_cores))

        # Reset index for CDF calculation
        df.reset_index(inplace=True, drop=True)
        print(f"-- Completed Tanimoto Calculation in {time.time() - then}")

        return df


def _parallelize_dataframe(df, func, n_cores=1):
    """Parallelize mapping a function to a dataframe.

    Applies a function to a dataframe in parallel by chunking it up over
    the specified number of cores.
    """
    # Require minimum number of compounds to parallelize
    if len(df) <= n_cores * 4:
        n_cores = 1

    if n_cores > 1:
        df_split = np.array_split(df, n_cores)
        pool = multiprocessing.Pool(n_cores)
        df = pd.concat(pool.map(func, df_split))
        pool.close()
        pool.join()
    else:
        df = func(df)
    return df


def _calc_max_T(t_df, min_T, df):
    """Calculate maximum tanimoto.

    Generate the tanimoto to use to generate the PMF to sample from.
    For each compound a list of tanimoito values are obtained by a generated
    compound to every target compound and the max is taken.
    """
    df["fp"] = df["SMILES"].map(get_fp)

    df["T"] = None
    fp = None
    for i in range(len(df)):
        fp = df["fp"].iloc[i]
        df["T"].iloc[i] = max(t_df["fp"].map(lambda x: FingerprintSimilarity(x, fp)))
    # Filter out low Tanimoto
    df = df[df["T"] > min_T]

    return df


# End Tanimoto Sampling Filter
###############################################################################

###############################################################################
# Metabolomics data filter


class MetabolomicsFilter(Filter):
    def __init__(
        self,
        filter_name,
        met_data_name,
        met_data_path,
        possible_adducts,
        mass_tolerance,
    ):
        """Load metabolomics data into a MetabolomicsDataset object."""

        self._filter_name = filter_name
        self.met_data_name = met_data_name

        if met_data_path:
            self.met_df = pd.read_csv(met_data_path).fillna("")
        else:
            self.met_df = None

        self.possible_adducts = possible_adducts
        self.mass_tolerance = mass_tolerance

        class Options:
            """Just need an object with an adducts property to pass to
            MetabolomicsDataset for initialization."""

            def __init__(self, adducts):
                self.adducts = adducts

        options = Options(possible_adducts)
        self.metabolomics_dataset = MetabolomicsDataset(self.met_data_name, options)

        # Load Metabolomics peaks
        for _, row in self.met_df.iterrows():

            smiles = row["Predicted Structure (smiles)"]
            if smiles:
                smiles = CanonSmiles(smiles)

                mol = MolFromSmiles(smiles)
                mol = utils.neutralise_charges(mol)
                inchi_key = MolToInchiKey(mol)
            else:
                mol = None
                inchi_key = None

            peak = Peak(
                name=row["Peak ID"],
                r_time=row["Retention Time"],
                mz=row["Aggregate M/Z"],
                charge=row["Polarity"].capitalize(),
                inchi_key=inchi_key,
            )

            if inchi_key:
                self.metabolomics_dataset.known_peaks.append(peak)
            else:
                self.metabolomics_dataset.unk_peaks.append(peak)

        # Calculate possible peak masses, they get saved to object
        self.metabolomics_dataset.enumerate_possible_masses(self.mass_tolerance)

    @property
    def filter_name(self):
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe, num_workers):
        """Choose compounds to expand based on whether they are found in a
        metabolomics dataset.

        Parameters
        ----------
        self : Pickaxe
            Instance of Pickaxe class.
        """
        # Get compounds eligible for expansion in the current generation
        compounds_to_check = []
        set_unreactive = True

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:

                cpd["Matched_Peak_IDs"] = []
                # Check for targets and only react if terminal
                if pickaxe.react_targets or not pickaxe.targets:
                    compounds_to_check.append(cpd)
                else:
                    for t_id in pickaxe.targets:
                        if "C" + t_id[1:] != cpd["_id"]:
                            compounds_to_check.append(cpd)
                            set_unreactive = False
                            break

                    if set_unreactive:
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                    else:
                        set_unreactive = True

        # Get compounds to keep
        cpd_info = [(cpd["_id"], cpd["SMILES"]) for cpd in compounds_to_check]

        possible_ranges = self.metabolomics_dataset.possible_ranges
        mass_matched_ids = self._filter_by_mass(pickaxe, cpd_info, possible_ranges)

        # Get compounds to remove
        ids = set(i[0] for i in cpd_info)
        cpds_remove_set = ids - mass_matched_ids

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set

    def _filter_by_mass(self, pickaxe, cpd_info, possible_ranges):
        """Check to see if compound masses each lie in any possible mass ranges."""
        matched_ids = set()
        for cpd in cpd_info:
            cpd_exact_mass = ExactMolWt(MolFromSmiles(cpd[1]))
            for possible_range in possible_ranges:
                if possible_range[0] < cpd_exact_mass < possible_range[1]:
                    c_id = cpd[0]
                    peak_id = possible_range[2]
                    matched_ids.add(c_id)
                    pickaxe.compounds[c_id]["Matched_Peak_IDs"].append(peak_id)
        return matched_ids

    def pre_print(self):
        print(f"Filtering compounds based on match with metabolomics data.")

    def post_print(self, pickaxe, n_total, n_filtered, time_sample):
        print(
            (
                f"{n_filtered} of {n_total} compounds selected after "
                f"Metabolomics filtering of generation {pickaxe.generation}"
                f"--took {round(time.time() - time_sample, 2)}s.\n"
            )
        )


# End metabolomics data filter
###############################################################################

###############################################################################
# Hard cutoff filters -- e.g. Tanimoto and MCS metric


class TanimotoFilter(Filter):
    """Calculate Tanimoto Coefficient for fingerprints between predicted and
    target compounds and filter out predicted compounds with low Tanimoto, set
    by specifying crit_tani."""

    def __init__(self, filter_name, crit_tani, increasing_tani):
        self._filter_name = filter_name
        self.crit_tani = crit_tani
        self.increasing_tani = increasing_tani

    @property
    def filter_name(self):
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe, num_workers=1):
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a tanimoto similarity score to a target
        compound greater than or equal to the crit_tani, for expansion.
        """

        if not pickaxe.target_fps:
            print("No targets to filter for. Can't expand.")
            return None

        # Set up variables required for filtering
        # Tanimoto Threshold
        if isinstance(self.crit_tani, list):
            crit_tani = self.crit_tani[pickaxe.generation]
        else:
            crit_tani = self.crit_tani

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = []

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:

                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    compounds_to_check.append(cpd)
                else:
                    # TODO this is not efficient
                    for t_id in pickaxe.targets:
                        if "C" + t_id[1:] != cpd["_id"]:
                            compounds_to_check.append(cpd)
                        else:
                            pickaxe.compounds[cpd["_id"]]["Expand"] = False

        # Run the filtering code to get a list of compounds to ignore
        print(
            f"Filtering Generation {pickaxe.generation} "
            f"with Tanimoto > {crit_tani}."
        )
        # Get input to filter code, c_id and smiles (to be
        # turned into fingerprint)
        cpd_info = [(cpd["_id"], cpd["SMILES"]) for cpd in compounds_to_check]
        this_gen_crit_tani = self.crit_tani[pickaxe.generation]
        cpd_filters = self._filter_by_tani_helper(
            cpd_info, pickaxe.target_fps, num_workers, this_gen_crit_tani
        )

        # Process filtering results
        cpds_remove_set = set()
        for c_id, current_tani in cpd_filters:
            # Check if tani is increasing
            if (
                self.increasing_tani
                and current_tani >= pickaxe.compounds[c_id]["last_tani"]
            ):
                pickaxe.compounds[c_id]["last_tani"] = current_tani
            if current_tani < this_gen_crit_tani:
                pickaxe.compounds[c_id]["Expand"] = False
                cpds_remove_set.add(c_id)

        return cpds_remove_set

    def _filter_by_tani_helper(
        self, compounds_info, target_fps, num_workers, this_crit_tani
    ):
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(0.1 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total * 100)} percent complete")

        # compound_info = [(smiles, id)]
        cpds_to_filter = list()
        compare_target_fps_partial = partial(
            self._compare_target_fps, target_fps, this_crit_tani
        )

        if num_workers > 1:
            # Set up parallel computing of compounds to
            chunk_size = max([round(len(compounds_info) / (num_workers * 4)), 1])
            pool = multiprocessing.Pool(num_workers)
            for i, res in enumerate(
                pool.imap_unordered(
                    compare_target_fps_partial, compounds_info, chunk_size
                )
            ):
                # If the result of comparison is false, compound is not expanded
                # Default value for a compound is True, so no need to
                # specify expansion
                if res:
                    cpds_to_filter.append(res)
                print_progress(i, len(compounds_info), "Tanimoto filter progress:")

        else:
            for i, cpd in enumerate(compounds_info):
                res = compare_target_fps_partial(cpd)
                if res:
                    cpds_to_filter.append(res)
                print_progress(i, len(compounds_info), "Tanimoto filter progress:")
        print("Tanimoto filter progress: 100 percent complete")

        return cpds_to_filter

    def _compare_target_fps(self, target_fps, this_crit_tani, compound_info):
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
            max_tani = 0
            for fp2 in target_fps:
                tani = AllChem.DataStructs.FingerprintSimilarity(fp1, fp2)
                if tani >= this_crit_tani:
                    return (compound_info[0], tani)
                elif tani >= max_tani:
                    max_tani = tani
            return (compound_info[0], max_tani)
            # TODO what except to use here?
        except:  # noqa
            return (compound_info[0], -1)

    def preprint(self, pickaxe):
        print(
            (
                "Filtering out compounds with maximum tanimoto match"
                f" < {self.crit_tani[pickaxe.generation]}"
            )
        )

    def post_print(self, pickaxe, n_total, n_filtered, time_sample):
        print(
            (
                f"{n_filtered} of {n_total} compounds selected after "
                f"Tanimoto filtering of generation {pickaxe.generation} "
                f"at cutoff of {self.crit_tani[pickaxe.generation]}. "
                f"--took {round(time.time() - time_sample, 2)}s.\n"
            )
        )


class MCSFilter(Filter):
    """Filter out compounds based on Maximum Common Substructure (MCS) with
    target compounds."""

    def __init__(self, filter_name, crit_mcs):
        self._filter_name = filter_name
        self.crit_mcs = crit_mcs

    @property
    def filter_name(self):
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe, num_workers=1):
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a tanimoto similarity score to a target
        compound greater than or equal to the crit_tani, for expansion.
        """

        if not pickaxe.target_fps:
            print("No targets to filter for. Can't expand.")
            return None

        # Set up variables required for filtering
        # Tanimoto Threshold
        if isinstance(self.crit_mcs, list):
            crit_mcs = self.crit_mcs[pickaxe.generation]
        else:
            crit_mcs = self.crit_mcs

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = []

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:

                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    compounds_to_check.append(cpd)
                else:
                    for t_id in pickaxe.targets:
                        if "C" + t_id[1:] != cpd["_id"]:
                            compounds_to_check.append(cpd)
                        else:
                            pickaxe.compounds[cpd["_id"]]["Expand"] = False

        # Run the filtering code to get a list of compounds to ignore
        print(f"Filtering Generation {pickaxe.generation}" " with MCS > {crit_mcs}.")
        # Get input to filter code, c_id and smiles
        # (to be turned into fingerprint)
        cpd_info = [(cpd["_id"], cpd["SMILES"]) for cpd in compounds_to_check]
        this_gen_crit_mcs = crit_mcs
        cpd_filters = self._filter_by_mcs_helper(
            cpd_info, pickaxe.target_smiles, num_workers, this_gen_crit_mcs
        )

        # Process filtering results
        keep_ids = [cpd[0] for cpd in cpd_filters]
        cpds_remove_set = set()
        for c_id, _ in cpd_info:
            if c_id not in keep_ids:
                pickaxe.compounds[c_id]["Expand"] = False
                cpds_remove_set.add(c_id)

        return cpds_remove_set

    def _filter_by_mcs_helper(
        self, compounds_info, target_smiles, num_workers, this_crit_mcs, retro=False
    ):
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(0.1 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total * 100)} percent complete")

        # compound_info = [(smiles, id)]
        cpds_to_filter = list()
        compare_target_mcs_partial = partial(
            self._compare_target_mcs, target_smiles, retro
        )

        if num_workers > 1:
            # Set up parallel computing of compounds to
            chunk_size = max([round(len(compounds_info) / (num_workers * 4)), 1])
            pool = multiprocessing.Pool(num_workers)
            for i, res in enumerate(
                pool.imap_unordered(
                    compare_target_mcs_partial,
                    compounds_info,
                    this_crit_mcs,
                    chunk_size,
                )
            ):

                if res:
                    cpds_to_filter.append(res)
                print_progress(
                    i,
                    len(compounds_info),
                    "Maximum Common Substructure filter progress:",
                )

        else:
            for i, cpd in enumerate(compounds_info):
                res = compare_target_mcs_partial(cpd, this_crit_mcs)
                if res:
                    cpds_to_filter.append(res)
                print_progress(
                    i,
                    len(compounds_info),
                    "Maximum Common Substructure filter progress:",
                )

        print("Maximum Common Substructure filter progress:" " 100 percent complete")
        return cpds_to_filter

    def _compare_target_mcs(self, target_smiles, retro, compound_info, this_crit_mcs):
        """Compare target MCS.

        Helper function to allow parallel computation of MCS filtering.
        Works with _filter_by_tani_helper

        Returns cpd_id if a the compound is similar enough to a target.

        """

        def get_mcs_overlap(mol, target_mol):
            mcs_out = mcs.FindMCS(
                [mol, target_mol], matchValences=False, ringMatchesRingOnly=False
            )

            if not mcs_out.canceled:
                ss_atoms = mcs_out.numAtoms
                ss_bonds = mcs_out.numBonds
                t_atoms = target_mol.GetNumAtoms()
                t_bonds = target_mol.GetNumBonds()

                mcs_overlap = (ss_atoms + ss_bonds) / (t_bonds + t_atoms)
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
                if mcs_overlap >= this_crit_mcs:
                    return (compound_info[0], mcs_overlap)
        # TODO what except to use here?
        except:  # noqa
            return (compound_info[0], -1)

    def preprint(self, pickaxe):
        print(
            (
                "Filtering out compounds with maximum MCS match"
                f" < {self.crit_mcs[pickaxe.generation]}"
            )
        )

    def post_print(self, pickaxe, n_total, n_filtered, time_sample):
        print(
            (
                f"{n_filtered} of {n_total} compounds selected after "
                f"MCS filtering of generation {pickaxe.generation} "
                f"at cutoff of {self.crit_mcs[pickaxe.generation]}. "
                f"--took {round(time.time() - time_sample, 2)}s.\n"
            )
        )


#
#
# End filters
###############################################################################
