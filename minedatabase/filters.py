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
from typing import Callable, List, Tuple

import numpy as np
import pandas as pd
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from mordred import Calculator, descriptors
from rdkit.Chem import AddHs, AllChem, CanonSmiles
from rdkit.Chem import rdFMCS as mcs
from rdkit.Chem.AllChem import RDKFingerprint
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.DataStructs import FingerprintSimilarity
from scipy.stats import rv_discrete

from minedatabase import utils

from minedatabase.metabolomics import MetabolomicsDataset, Peak
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import Chunks, get_fp


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")

###############################################################################
# ABC for all Filter Subclasses


class Filter(metaclass=abc.ABCMeta):
    """Abstract base class used to generate filters.

    The Filter class provides the framework for interaction with pickaxe expansions.
    Each filter subclass must inherit properties from the Filter class.
    All subclasses must implement properties and methods decorated with
    @abc.abstractmethod. Feel free to override other non-private methods as
    well, such as _pre_print() and _post_print().
    """

    @property
    @abc.abstractmethod
    def filter_name(self):
        """Obtain name of filter."""
        pass

    @abc.abstractmethod
    def _choose_cpds_to_filter(self, pickaxe, processes):
        """Return list of compounds to remove from pickaxe object."""
        pass

    def apply_filter(self, pickaxe: Pickaxe, processes: int = 1, print_on: bool = True):
        """Apply filter from Pickaxe object.

        Parameters
        ----------
        pickaxe : Pickaxe
            The Pickaxe object to filter.
        processes : int
            The number of processes to use, by default 1.
        print_on : bool
            Whether or not to print filtering results.
        """
        time_sample = time.time()

        if print_on:
            n_total = self._get_n(pickaxe, "total")
            self._pre_print_header(pickaxe)
            self._pre_print()

        compound_ids_to_check = self._choose_cpds_to_filter(pickaxe, processes)

        if compound_ids_to_check:
            self._apply_filter_results(pickaxe, compound_ids_to_check)

        if print_on:
            n_filtered = self._get_n(pickaxe, "filtered")
            self._post_print(pickaxe, n_total, n_filtered, time_sample)
            self._post_print_footer(pickaxe)

    def _pre_print_header(self, pickaxe):
        """Print header before filtering."""
        print("----------------------------------------")
        print(f"Filtering Generation {pickaxe.generation}\n")

    def _pre_print(self):
        """Print filter being applied."""
        print(f"Applying filter: {self.filter_name}")

    def _post_print(self, pickaxe, n_total, n_filtered, time_sample):
        """Print results of filtering."""
        print(
            f"{n_filtered} of {n_total} compounds remain after applying "
            f"filter: {self.filter_name}"
            f"--took {round(time.time() - time_sample, 2)}s.\n"
        )

    def _post_print_footer(self, pickaxe):
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
            3. Start with "C"
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
    """Filter that samples randomly from weighted tanimoto.

    TanimotoSamplingFilter takes a distribution of tanimoto similarity scores and uses
    inverse CDF sampling to select N compounds for further expansion. Each compound
    is assigned a tanimoto similarity score that corresponds to the maximum tanimoto
    score of the set of similarity scores obtained by comparing that compound to each
    target. These scores can also be weighted by a specified function to bias higher
    or lower tanimoto scores.

    Parameters
    ----------
    sample_size : int
        Number of compounds to sample.
    weight : Callable
        Function to weight the tanimoto similarity score with.

    Attributes
    ----------
    sample_size : int
        Number of compounds to sample.
    weight : Callable
        Function to weight the tanimoto similarity score with.
    """

    def __init__(self, sample_size: int, weight: Callable = None):
        self._filter_name = "Tanimoto Sampling Filter"
        self.sample_size = sample_size
        self.sample_weight = weight

    @property
    def filter_name(self):
        return self._filter_name

    def _pre_print(self):
        """Print before filtering."""
        print(
            (
                f"Sampling {self.sample_size} Compounds Based on a "
                f"Weighted Tanimoto Distribution"
            )
        )

    def _post_print(self, pickaxe, n_total, n_filtered, time_sample):
        """Print after filtering."""
        print(
            (
                f"{n_filtered} of {n_total} "
                "compounds selected after "
                f"Tanimoto Sampling of generation {pickaxe.generation}"
                f"--took {time.time() - time_sample}s.\n"
            )
        )

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int):
        """
        Samples N compounds to expand based on the weighted Tanimoto distribution.

        Parameters
        ----------
        pickaxe : Pickaxe
            Pickaxe object to filter
        processes : int
            Number of processes to use.

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

        sampled_ids = self._sample_by_tanimoto(
            cpd_info,
            pickaxe.target_fps,
            self.sample_size,
            min_T=0.15,
            weighting=self.sample_weight,
            max_iter=None,
            processes=processes,
        )

        # Get compounds to remove
        ids = set(i[0] for i in cpd_info)
        cpds_remove_set = ids - sampled_ids

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set

    def _sample_by_tanimoto(
        self,
        mol_info: List[Tuple[str, str]],
        t_fp: RDKFingerprint,
        n_cpds: int = None,
        min_T: float = 0.05,
        weighting: Callable = None,
        max_iter: int = None,
        processes: int = 1,
    ) -> List[str]:
        """Smple compounds by weighted tanimoto coefficient.

        Use inverse cumulative distrbution function (CDF) sampling to select
        compounds based on a weighted tanimoto coefficient distribution.

        Parameters
        ----------
        mol_info : List[Tuple[str, str]]
            A list consisting of (compound_id, SMILES).
        t_fp : List[RDKFingerprint]
            Target fingerprints to compare compounds to.
        n_cpds : int, optional
            Number of compounds to select for sampling, by default None.
        min_T : float, optional
            Minimum Tanimoto similarity to be considered for sampling, by default 0.05.
        weighting : Callable, optional
            Function that accepts a Tanimoto similarity score and returns
            a float, by default None.
        max_iter : int, optional
            The maximum number of iterations before regenerating the CDF for sampling
            , by default None.
        processes : int, optional
            Number of processes to use, by default 1.

        Returns
        -------
        List[str]
            The compound ids to expand.
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
        df = self._gen_df_from_tanimoto(
            mol_info, t_fp, min_T=min_T, processes=processes
        )
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

    def _gen_rv_from_df(
        self,
        df: pd.DataFrame,
        chosen: List = [],
        weighting: Callable = None,
    ) -> rv_discrete:
        """Generate a scipy.rv object to sample from the inverse CDF

        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing the data to sample
        chosen : List, optional
            Compound ids that have already been chosen, by default []
        weighting : Callable, optional
            Function to weight the Tanimoto distribution by, by default None

        Returns
        -------
        rv_discrete
            scipy.rv object to sample from.
        """
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

    def _gen_df_from_tanimoto(
        self,
        mol_info: List[Tuple[str, str]],
        t_fp: List[RDKFingerprint],
        min_T: float = 0.05,
        processes: int = 1,
    ):
        """Generate a dataframe from tanimoto

        Parameters
        ----------
        mol_info : List[Tuple[str, str]]
            A list consisting of (compound_id, SMILES)
        t_fp : List[RDKFingerprint]
            Target fingerprints to compare compounds to.
        min_T : float, optional
            Minimum Tanimoto similarity to be considered for sampling, by default 0.05
        processes : int, optional
            Number of processes to use, by default 1
        """

        then = time.time()
        print("-- Calculating Fingerprints and Tanimoto Values.")
        # target fingerprint dataframe
        t_df = pd.DataFrame(t_fp, columns=["fp"])

        # Calculate Tanimoto for each compound and drop T < min_T
        partial_T_calc = partial(_calc_max_T, t_df, min_T)

        df = pd.DataFrame()
        for mol_chunk in Chunks(mol_info, 10000):
            # Construct targets to sample df
            temp_df = pd.DataFrame(mol_chunk, columns=["_id", "SMILES"])
            df = df.append(_parallelize_dataframe(temp_df, partial_T_calc, processes))

        # Reset index for CDF calculation
        df.reset_index(inplace=True, drop=True)
        print(f"-- Completed Tanimoto Calculation in {time.time() - then}")

        return df


def _parallelize_dataframe(
    df: pd.DataFrame,
    func: Callable,
    processes: int = 1,
) -> pd.DataFrame:
    """Parallelize mapping a function to a dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to apply function to.
    func : Callable
        Function to map onto dataframe.
    processes : int
        Number of processes to use, by default 1.
    """
    # Require minimum number of compounds to parallelize
    if len(df) <= processes * 4:
        processes = 1

    if processes > 1:
        df_split = np.array_split(df, processes)
        pool = multiprocessing.Pool(processes)
        df = pd.concat(pool.map(func, df_split))
        pool.close()
        pool.join()
    else:
        df = func(df)
    return df


def _calc_max_T(t_df: pd.DataFrame, min_T: float, df: pd.DataFrame) -> pd.DataFrame:
    """Calculate maximum tanimoto.

    Generate the tanimoto to use to generate the PMF to sample from.
    For each compound a list of tanimoito values are obtained by a generated
    compound to every target compound and the max is taken.

    Parameters
    ----------
    t_df : pd.Dataframe
        Dataframe containing the target fingerprints.
    min_T : float
        The minimum tanimoto similarity score needed to consider a compound.
    df : pd.DataFrame
        Dataframe to calculate the max tanimoto for.
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
    """A metabolomics filter short description.

    A longer description.

    Parameters
    ----------
    please_jon : float
        Fill me in

    Attributes
    ----------
    I_believe_in_you : bool
        Do I believe in you?, by default True.
    """

    def __init__(
        self,
        filter_name,
        met_data_name,
        met_data_path,
        possible_adducts,
        mass_tolerance,
        rt_predictor=None,
        rt_threshold=None,
        rt_important_features=None,
    ):
        """Load metabolomics data into a MetabolomicsDataset object."""

        self._filter_name = filter_name
        self.met_data_name = met_data_name

        self.rt_predictor = rt_predictor
        self.rt_threshold = rt_threshold
        self.rt_important_features = rt_important_features

        if self.rt_predictor and self.rt_threshold:
            self.filter_by_rt = True
            self.fp_calculator = Calculator(descriptors, ignore_3D=False)
        else:
            self.filter_by_rt = False
            self.fp_calculator = None

        if met_data_path:
            self.met_df = pd.read_csv(met_data_path).fillna("")
        else:
            self.met_df = None

        self.possible_adducts = possible_adducts
        self.mass_tolerance = mass_tolerance

        options = MetabolomicsOptions(possible_adducts)
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

    def _choose_cpds_to_filter(self, pickaxe, processes):
        """Choose compounds to expand based on whether they are found in a
        metabolomics dataset.

        Parameters
        ----------
        self : Pickaxe
            Instance of Pickaxe class.
        """
        if pickaxe.generation == 0:
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

                cpd["Matched_Peak_IDs"] = []
                cpd["Matched_Adducts"] = []

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

        filter_by_mass_and_rt_partial = partial(
            self._filter_by_mass_and_rt, possible_ranges
        )

        mass_matched_ids = set()
        cpd_met_dict = {}

        if processes > 1:
            # Set up parallel computing
            chunk_size = max([round(len(cpd_info) / (processes * 4)), 1])
            pool = multiprocessing.Pool(processes)

            for res in pool.imap_unordered(
                filter_by_mass_and_rt_partial, cpd_info, chunk_size
            ):
                if res[0]:
                    this_cpd_id = res[0]
                    mass_matched_ids.add(this_cpd_id)
                    this_cpd_met_dict = res[1]
                    cpd_met_dict[this_cpd_id] = this_cpd_met_dict

        else:
            for cpd in cpd_info:
                res = filter_by_mass_and_rt_partial(cpd)
                if res[0]:
                    mass_matched_ids.add(res[0])
                    cpd_met_dict[res[0]] = res[1]

        for c_id in mass_matched_ids:
            pickaxe.compounds[c_id]["Matched_Peak_IDs"] += cpd_met_dict[c_id][
                "Matched_Peak_IDs"
            ]
            pickaxe.compounds[c_id]["Matched_Adducts"] += cpd_met_dict[c_id][
                "Matched_Adducts"
            ]
            pickaxe.compounds[c_id]["Predicted_RT"] = cpd_met_dict[c_id]["Predicted_RT"]

        # Get compounds to remove
        ids = set(i[0] for i in cpd_info)
        cpds_remove_set = ids - mass_matched_ids

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set

    def _filter_by_mass_and_rt(self, possible_ranges, cpd_info):
        """Check to see if compound masses  (and optionally, retention time)
        each lie in any possible mass ranges."""
        c_id_if_matched = None
        cpd_dict = {"Predicted_RT": None, "Matched_Peak_IDs": [], "Matched_Adducts": []}

        cpd_exact_mass = ExactMolWt(MolFromSmiles(cpd_info[1]))
        predicted_rt = None
        for possible_range in possible_ranges:
            if possible_range[0] < cpd_exact_mass < possible_range[1]:
                c_id = cpd_info[0]
                smiles = cpd_info[1]
                peak_id = possible_range[2]
                adduct = possible_range[3]

                if self.filter_by_rt:
                    if not predicted_rt:
                        predicted_rt = self._predict_rt(smiles)
                    if not predicted_rt:
                        # sometimes can't predict RT due to missing vals in fingerprint
                        continue

                    expt_rt = self.metabolomics_dataset.get_rt(peak_id)
                    if not expt_rt:
                        raise ValueError(f"No retention time found for peak, {peak_id}")

                    cpd_dict["Predicted_RT"] = predicted_rt
                    if abs(expt_rt - predicted_rt) > self.rt_threshold:
                        continue  # if outside threshold, don"t add to matched peaks

                c_id_if_matched = c_id
                cpd_dict["Matched_Peak_IDs"].append(peak_id)
                cpd_dict["Matched_Adducts"].append(adduct)

        return c_id_if_matched, cpd_dict

    def _predict_rt(self, smiles):
        """Predict Retention Time from SMILES string using provided predictor."""
        mol = MolFromSmiles(smiles)
        mol = AddHs(mol)

        fp = self.fp_calculator(mol)
        # Transform dict into array of values (fingerprint)
        if self.rt_important_features:
            fp = np.array(
                [fp[feature] for feature in self.rt_important_features]
            ).reshape(1, -1)

        def validate_np_val(val):
            """Make sure value is numeric, not NaN, and not infinity."""
            if isinstance(val, float) and not np.isnan(val) and not np.isinf(val):
                return True
            return False

        if all([validate_np_val(val) for val in fp[0]]):
            predicted_rt = self.rt_predictor.predict(fp)[0]
        else:
            return None

        return predicted_rt

    def _pre_print(self):
        print(f"Filtering compounds based on match with metabolomics data.")

    def _post_print(self, pickaxe, n_total, n_filtered, time_sample):
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
    """A filter that uses tanimoto similarity score to determine compounds to expand.

    TanimotoFilter applies a strict cutoff to to the tanimoto similarity score of
    compounds to determine which compounds to expand.

    Parameters
    ----------
    crit_tani : float
        The tanimoto similarity score threshold.
    increasing_tani : bool
        Whether or not to only keep compounds whos tanimoto score is higher than its
        parent.

    Attributes
    ----------
    crit_tani : float
        The tanimoto similarity score threshold.
    increasing_tani : bool
        Whether or not to only keep compounds whos tanimoto score is higher than its
        parent.
    """

    def __init__(self, crit_tani: float, increasing_tani: bool):
        self._filter_name = "Tanimoto Cutoff"
        self.crit_tani = crit_tani
        self.increasing_tani = increasing_tani

    @property
    def filter_name(self):
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int = 1):
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
        if type(self.crit_tani) in [list, tuple]:
            if len(self.crit_tani) - 1 < pickaxe.generation:
                crit_tani = self.crit_tani[-1]
            else:
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
        if type(self.crit_tani) in [list, tuple]:
            if len(self.crit_tani) - 1 < pickaxe.generation:
                this_gen_crit_tani = self.crit_tani[-1]
            else:
                this_gen_crit_tani = self.crit_tani[pickaxe.generation]
        else:
            this_gen_crit_tani = self.crit_tani
        cpd_filters = self._filter_by_tani_helper(
            cpd_info, pickaxe.target_fps, processes, this_gen_crit_tani
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
        self, compounds_info, target_fps, processes, this_crit_tani
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

        if processes > 1:
            # Set up parallel computing of compounds to
            chunk_size = max([round(len(compounds_info) / (processes * 4)), 1])
            pool = multiprocessing.Pool(processes)
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
        if type(self.crit_tani) in [list, tuple]:
            print_tani = self.crit_tani[pickaxe.generation]
        else:
            print_tani = self.crit_tani
        print(f"Filtering out compounds with maximum tanimoto match < {print_tani}")

    def _post_print(self, pickaxe, n_total, n_filtered, time_sample):
        if type(self.crit_tani) in [list, tuple]:
            if len(self.crit_tani) - 1 < pickaxe.generation:
                print_tani = self.crit_tani[-1]
            else:
                print_tani = self.crit_tani[pickaxe.generation]
        else:
            print_tani = self.crit_tani
        print(
            (
                f"{n_filtered} of {n_total} compounds selected after "
                f"Tanimoto filtering of generation {pickaxe.generation} "
                f"at cutoff of {print_tani}. "
                f"--took {round(time.time() - time_sample, 2)}s.\n"
            )
        )


class MCSFilter(Filter):
    """A filter that uses MCS score to determine compounds to expand.

    MCSFilter applies a strict cutoff to to the MCS score of
    compounds to determine which compounds to expand.

    Parameters
    ----------
    crit_mcs: float
        The maximum common substructure similarity score threshold.

    Attributes
    ----------
    crit_mcs : float
        The maximum common substructure similarity score threshold.
    """

    def __init__(self, crit_mcs: float):
        self._filter_name = "MCS Cutoff"
        self.crit_mcs = crit_mcs

    @property
    def filter_name(self):
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe, processes=1):
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a tanimoto similarity score to a target
        compound greater than or equal to the crit_tani, for expansion.
        """

        if not pickaxe.target_fps:
            print("No targets to filter for. Can't expand.")
            return None

        # Set up variables required for filtering
        # MCS Threshold
        if type(self.crit_mcs) in [list, tuple]:
            if len(self.crit_mcs) - 1 < pickaxe.generation:
                crit_mcs = self.crit_mcs[-1]
            else:
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
            cpd_info, pickaxe.target_smiles, processes, this_gen_crit_mcs
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
        self, compounds_info, target_smiles, processes, this_crit_mcs, retro=False
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

        if processes > 1:
            # Set up parallel computing of compounds to
            chunk_size = max([round(len(compounds_info) / (processes * 4)), 1])
            pool = multiprocessing.Pool(processes)
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
        if type(self.crit_mcs) in [list, tuple]:
            if len(self.crit_mcs) - 1 < pickaxe.generation:
                crit_mcs = self.crit_mcs[-1]
            else:
                crit_mcs = self.crit_mcs[pickaxe.generation]
        else:
            crit_mcs = self.crit_mcs
        print(f"Filtering out compounds with maximum MCS match < {crit_mcs}")

    def _post_print(self, pickaxe, n_total, n_filtered, time_sample):
        if type(self.crit_mcs) in [list, tuple]:
            if len(self.crit_mcs) - 1 < pickaxe.generation:
                crit_mcs = self.crit_mcs[-1]
            else:
                crit_mcs = self.crit_mcs[pickaxe.generation]
        else:
            crit_mcs = self.crit_mcs
        print(
            (
                f"{n_filtered} of {n_total} compounds selected after "
                f"MCS filtering of generation {pickaxe.generation} "
                f"at cutoff of {crit_mcs}. "
                f"--took {round(time.time() - time_sample, 2)}s.\n"
            )
        )


#
#
# End filters
###############################################################################
