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
from typing import Callable, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import periodictable
from rdkit import DataStructs
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from mordred import Calculator, descriptors
from rdkit.Chem import AddHs, AllChem, CanonSmiles, RDKFingerprint
from rdkit.Chem import rdFMCS as mcs
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.DataStructs import FingerprintSimilarity
from scipy.stats import rv_discrete
from sklearn.ensemble import RandomForestRegressor
from minedatabase import pickaxe

from minedatabase.metabolomics import MetabolomicsDataset, Peak
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import Chunks, get_fp, neutralise_charges


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
    def filter_name(self) -> str:
        """Obtain name of filter."""
        pass

    @abc.abstractmethod
    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
        """Return list of compounds to remove from pickaxe object.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        processes : int
            The number of processes to use, by default 1.
        """
        pass

    def apply_filter(
        self, pickaxe: Pickaxe, processes: int = 1, print_on: bool = True
    ) -> None:
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

    def _pre_print_header(self, pickaxe: Pickaxe) -> None:
        """Print header before filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        """
        print("----------------------------------------")
        print(f"Filtering Generation {pickaxe.generation}\n")

    def _pre_print(self) -> None:
        """Print filter being applied."""
        print(f"Applying filter: {self.filter_name}")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
        """Print results of filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
            Unused here, but may be useful in your implementation.
        n_total : int
            Total number of compounds.
        n_filtered : int
            Number of compounds remaining after filtering.
        times_sample : float
            Time in seconds from time.time().
        """
        print(
            f"{n_filtered} of {n_total} compounds remain after applying "
            f"filter: {self.filter_name}"
            f"--took {round(time.time() - time_sample, 2)}s.\n"
        )

    def _post_print_footer(self, pickaxe: Pickaxe) -> None:
        """Print end of filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        """
        print(f"Done filtering Generation {pickaxe.generation}")
        print("----------------------------------------\n")

    def _get_n(self, pickaxe: Pickaxe, n_type: str) -> int:
        """Get current number of compounds to be filtered.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        n_type : str
            Whether to return "total" number of "filtered" number of compounds.

        Returns
        -------
        n : int
            Either the total or filtered number of compounds.
        """
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

    def _apply_filter_results(
        self, pickaxe: Pickaxe, compound_ids_to_check: List[str]
    ) -> None:
        """Apply filter results to Pickaxe object.

        Remove compounds and reactions that can be removed.
        For a compound to be removed it must:
            1. Not be flagged for expansion
            2. Not have a coproduct in a reaction marked for expansion
            3. Start with "C"

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network,
            this method modifies the Pickaxe object's compound documents.
        compound_ids_to_check : List[str]
            List of compound IDs to try to remove, if possible.
        """

        def should_delete_reaction(rxn_id: str) -> bool:
            """Whether we should delete reaction with supplied ID.

            Parameters
            ----------
            rxn_id : str
                ID of reaction.

            Returns
            -------
            bool
                True if we should delete, False otherwise.
            """
            products = pickaxe.reactions[rxn_id]["Products"]
            for _, c_id in products:
                if c_id.startswith("C") and c_id not in cpds_to_remove:
                    return False
            # Every compound isn't in cpds_to_remove
            return True

        def get_compounds_to_check_from_ids(
            pickaxe: Pickaxe, cpd_ids_to_check: List[str]
        ) -> List[Dict]:
            """Get compound documents from their IDs

            Parameters
            ----------
            pickaxe : Pickaxe
                Instance of Pickaxe being used to expand and filter the network.
            cpd_ids_to_check : List[str]
                List of compound IDs to get compound documents for.

            Returns
            -------
            cpds_to_check : List[Dict]
                List of compound documents.
            """
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


class SimilaritySamplingFilter(Filter):
    """Filter that samples randomly from weighted Tanimoto.

    TanimotoSamplingFilter takes a distribution of similarity scores and uses
    inverse CDF sampling to select N compounds for further expansion. Each compound
    is assigned a similarity score that corresponds to the maximum Tanimoto
    score of the set of similarity scores obtained by comparing that compound to each
    target. These scores can also be weighted by a specified function to bias higher
    or lower similarity scores scores.

    Types of fingerprints:
        1) RDKit (default)
        2) Morgan

    Types of similarity score:
        1) Tanimoto (default)
        2) Dice

    Parameters
    ----------
    sample_size : int
        Number of compounds to sample.
    weight : Callable
        Function to weight the Tanimoto similarity score with.
    fingerprint_method : str
        Method by which to calculate fingerprints. Options are RDKitFingerprint and
        Morgan, by default RDKitFingerprint.
    fingerprint_args : dict
        A dictionary of keyword arguments for the fingerprint functiono, by default None.
    similarity_method : str
        Method by which to calculate the similarity score. Options are Tanimoto and Dice,
        by default Tanimoto.

    Attributes
    ----------
    sample_size : int
        Number of compounds to sample.
    sample_weight : Callable
        Function to weight the Tanimoto similarity score with.
    fingerprint_method : str
        Fingerprint calculation method.
    fingerprint_args : dict
        Arguments for fingerprint method.
    similairty_method : str
        Simlarity calulation method.
    target_fps : List
        List of fingerprints.
    """

    def __init__(
        self,
        sample_size: int,
        weight: Callable = None,
        fingerprint_method: str = None,
        fingerprint_args: str = None,
        similarity_method: str = None,
    ) -> None:
        self._filter_name = "Tanimoto Sampling Filter"
        self.sample_size = sample_size
        self.sample_weight = weight

        self.fingerprint_method = fingerprint_method or "RDKit"
        self.fingerprint_args = fingerprint_args or {}
        self.similarity_method = similarity_method or "Tanimoto"

        self.target_fps = []

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _set_target_fps(self, pickaxe: Pickaxe):
        self.target_fps = []
        for smiles in pickaxe.target_smiles:
            mol = MolFromSmiles(smiles)
            if self.fingerprint_method == "Morgan":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, **self.fingerprint_args)
            else:
                fp = RDKFingerprint(mol, **self.fingerprint_args)

            self.target_fps.append(fp)

    def _pre_print(self) -> None:
        """Print before filtering."""
        print(
            (
                f"Sampling {self.sample_size} Compounds Based on a "
                f"Weighted Tanimoto Distribution"
            )
        )

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
        """Print after filtering."""
        print(
            (
                f"{n_filtered} of {n_total} "
                "compounds selected after "
                f"Tanimoto Sampling of generation {pickaxe.generation}"
                f"--took {time.time() - time_sample}s.\n"
            )
        )

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
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

        self._set_target_fps(pickaxe)
        if not self.target_fps:
            print("No targets to filter for. Can't expand.")
            return None

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = []
        set_unreactive = False

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
                        if "C" + t_id[1:] == cpd["_id"]:
                            set_unreactive = True
                            break

                    if set_unreactive:
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                        set_unreactive = False
                    else:
                        compounds_to_check.append(cpd)

        # Get compounds to keep
        cpd_info = [(cpd["_id"], cpd["SMILES"]) for cpd in compounds_to_check]

        sampled_ids = self._sample_by_similarity(
            cpd_info,
            self.target_fps,
            self.sample_size,
            min_T=0.15,
            weighting=self.sample_weight,
            max_iter=None,
            processes=processes,
        )

        # TODO remove me
        print("num sampled = ", len(sampled_ids))
        # Get compounds to remove
        ids = set(i[0] for i in cpd_info)
        cpds_remove_set = ids - sampled_ids

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set

    def _sample_by_similarity(
        self,
        mol_info: List[Tuple[str, str]],
        t_fp: List[AllChem.RDKFingerprint],
        n_cpds: int = None,
        min_T: float = 0.05,
        weighting: Callable = None,
        max_iter: int = None,
        processes: int = 1,
    ) -> List[str]:
        """Sample compounds by weighted Tanimoto coefficient.

        Use inverse cumulative distrbution function (CDF) sampling to select
        compounds based on a weighted Tanimoto coefficient distribution.

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
        """Generate a scipy.rv object to sample from the inverse CDF.

        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing the data to sample.
        chosen : List, optional
            Compound ids that have already been chosen, by default [].
        weighting : Callable, optional
            Function to weight the Tanimoto distribution by, by default None.

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
        t_fp: List[AllChem.RDKFingerprint],
        min_T: float = 0.05,
        processes: int = 1,
    ) -> pd.DataFrame:
        """Generate a dataframe from Tanimoto

        Parameters
        ----------
        mol_info : List[Tuple[str, str]]
            A list consisting of (compound_id, SMILES).
        t_fp : List[RDKFingerprint]
            Target fingerprints to compare compounds to.
        min_T : float, optional
            Minimum Tanimoto similarity to be considered for sampling, by default 0.05.
        processes : int, optional
            Number of processes to use, by default 1.
        """

        then = time.time()
        print("-- Calculating Fingerprints and Tanimoto Values.")
        # target fingerprint dataframe
        t_df = pd.DataFrame(t_fp, columns=["fp"])

        # Calculate Tanimoto for each compound and drop T < min_T
        partial_T_calc = partial(
            _calc_max_T,
            t_df,
            min_T,
            self.fingerprint_method,
            self.fingerprint_args,
            self.similarity_method,
        )

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

    Returns
    -------
    df : pd.DataFrame
        New dataframe after having function applied to it in parallel.
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


def _calc_max_T(
    t_df: pd.DataFrame,
    min_T: float,
    fingerprint_method: str,
    fingerprint_args: str,
    similarity_method: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate maximum Tanimoto.

    Generate the Tanimoto to use to generate the PMF to sample from.
    For each compound a list of tanimoito values are obtained by a generated
    compound to every target compound and the max is taken.

    Parameters
    ----------
    t_df : pd.Dataframe
        Dataframe containing the target fingerprints.
    min_T : float
        The minimum Tanimoto similarity score needed to consider a compound.
    fingerprint_method: str
        Which fingerprint method to use, suppoorts RDKit and Morgan, by default RDKit.
    fingerprint_args: dict
        Keyword arguments to pass to fingerprint function, by default empty dict.
    similarity_method: str
        Which similarity method to use. Supports Tanimotoo and Dice.
    df : pd.DataFrame
        Dataframe to calculate the max Tanimoto for.

    Returns
    -------
    df : pd.Dataframe
        New dataframe with max Tanimoto values calculated.
    """

    def fingerprint(fingerprint_method, keyword_dict, smi):
        mol = AllChem.MolFromSmiles(smi)
        if fingerprint_method == "Morgan":
            return AllChem.GetMorganFingerprintAsBitVect(mol, **keyword_dict)
        else:
            return RDKFingerprint(mol, **keyword_dict)

    def similarity(similarity_method, fp1, fp2):
        if similarity_method == "dice":
            return FingerprintSimilarity(fp1, fp2, metric=DataStructs.DiceSimilarity)
        else:
            return FingerprintSimilarity(fp1, fp2)

    fingerprint_partial = partial(fingerprint, fingerprint_method, fingerprint_args)

    df["fp"] = df["SMILES"].map(fingerprint_partial)

    df["T"] = None
    fp = None
    for i in range(len(df)):
        fp = df["fp"].iloc[i]
        df["T"].iloc[i] = max(
            t_df["fp"].map(lambda x: similarity(similarity_method, x, fp))
        )
    # Filter out low Tanimoto
    df = df[df["T"] > min_T]

    return df


# End Tanimoto Sampling Filter
###############################################################################

###############################################################################
# Metabolomics data filter


class MetabolomicsFilter(Filter):
    """Filters out compounds that don't align with a metabolomics dataset.

    This filter compares the masses (and optionally, predicted retention times)
    of MINE compounds against peak masses (and retention times) in a
    metabolomics dataset. Tolerances for mass (in Da) and retention times
    (in units consistent with dataset) are specified by the user. If a
    compound's mass (and predicted retention time, if desired) does not match
    that for any peak in the dataset, it is filtered out.

    Parameters
    ----------
    filter_name : str
        Name of the filter, should be unique.
    met_data_name : str
        Name of the metabolomics dataset.
    met_data_path : str
        Path to metabolomics data CSV file with list of peak masses/RTs/etc.
    possible_adducts : List[str]
        List of possible adducts, see data/adducts for options.
    mass_tolerance : float
        Mass tolerance for peak matching in daltons.
    rt_predictor : sklearn.ensemble.RandomForestRegressor, optional
        Random forest regression model that takes a subset of a compound's 2D
        mordred fingerprint values (specified by rt_important_features) as
        input, defaults to None.
    rt_threshold : float, optional
        Retention time tolerance for peak matching in whatever units are used
        in the metabolomics dataset (e.g. seconds, minutes, etc.), defaults to
        None.
    rt_important_features : List[str], optional
        List of mordred descriptors to use as input into rt_predictor, make
        sure that the order is the same as how the model was trained, defaults
        to None.

    Attributes
    ----------
    filter_by_rt : Bool
        Whether the filter will filter by both mass and retention time (RT).
    fp_calculator : mordred.calculator.Calculator
        Calculator loaded with provided mordred descriptors.
    met_df : pd.DataFrame
        Dataframe containing metabolomics peak data.
    metabolomics_dataset : minedatabase.metabolomics.MetabolomicsDataset
        Instance of MetabolomicsDataset with loaded metabolomics data.
    """

    def __init__(
        self,
        filter_name: str,
        met_data_name: str,
        met_data_path: str,
        possible_adducts: List[str],
        mass_tolerance: float,
        rt_predictor: RandomForestRegressor = None,
        rt_threshold: float = None,
        rt_important_features: List[str] = None,
    ) -> None:
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

        self.metabolomics_dataset = MetabolomicsDataset(
            name=self.met_data_name,
            adducts=self.possible_adducts,
            tolerance=self.mass_tolerance,
        )
        self.metabolomics_dataset.known_peaks = []
        self.metabolomics_dataset.unknown_peaks = []

        # Load Metabolomics peaks
        for _, row in self.met_df.iterrows():

            smiles = row["Predicted Structure (smiles)"]
            if smiles:
                smiles = CanonSmiles(smiles)

                mol = MolFromSmiles(smiles)
                mol = neutralise_charges(mol)
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
                self.metabolomics_dataset.unknown_peaks.append(peak)

        # Calculate possible peak masses, they get saved to object
        self.metabolomics_dataset.enumerate_possible_masses(self.mass_tolerance)

    @property
    def filter_name(self) -> str:
        """Return filter name.

        Returns
        -------
        str
            Name of filter.
        """
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
        """Choose compounds to expand based on whether they are found in a
        metabolomics dataset.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe class.
        processes : int
            Number of processes (uses parallelization if > 1).

        Returns
        -------
        cpds_remove_set : Set[str]
            Set of IDs for compounds to try to remove from the expansion.
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

    def _filter_by_mass_and_rt(
        self,
        possible_ranges: List[Tuple[float, float, str, str]],
        cpd_info: List[Tuple[str]],
    ) -> Tuple[Optional[str], Dict]:
        """Check to see if compound masses  (and optionally, retention time)
        each lie in any possible mass ranges.

        Parameters
        ----------
        possible_ranges : List[Tuple[float, float, str, str]]
            Possible mass ranges based on peak masses and tolerance.
        cpd_info : List[Tuple[str]]
            Tuple of compound ID, SMILES, peak ID, and adduct name.

        Returns
        -------
        c_id_if_matched : str, optional
            Contains the compound ID if a hit is found, None by default.
        cpd_dict : Dict
            Contains predicted retention time, matched peak IDs (if any), and
            matched adduct names (if any).
        """
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

    def _predict_rt(self, smiles: str) -> Optional[float]:
        """Predict Retention Time from SMILES string using provided predictor.

        Parameters
        ----------
        smiles : str
            SMILES string of input compound.

        Returns
        -------
        predicted_rt : Optional[float]
            Predicted retention time, None if errors occur during prediction,
            for example if certain features of the input compound that are
            required for the prediction cannot be calculated.
        """
        mol = MolFromSmiles(smiles)
        mol = AddHs(mol)

        fp = self.fp_calculator(mol)
        # Transform dict into array of values (fingerprint)
        if self.rt_important_features:
            fp = np.array(
                [fp[feature] for feature in self.rt_important_features]
            ).reshape(1, -1)

        def validate_np_val(val: float) -> bool:
            """Make sure value is numeric, not NaN, and not infinity.

            Parameters
            ----------
            val : float
                Value to check.

            Returns
            -------
            bool
                True if input value is numeric, False otherwise.
            """
            if isinstance(val, float) and not np.isnan(val) and not np.isinf(val):
                return True
            return False

        if all([validate_np_val(val) for val in fp[0]]):
            predicted_rt = self.rt_predictor.predict(fp)[0]
        else:
            return None

        return predicted_rt

    def _pre_print(self) -> None:
        print(f"Filtering compounds based on match with metabolomics data.")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
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
# Cutoff filters -- e.g. Tanimoto, MCS metric, and molecular weight


class SimilarityFilter(Filter):
    """A filter that uses a similarity score to determine compounds to expand.

    SimilarityFilter applies a strict cutoff to to the similarity score of
    compounds to determine which compounds to expand.

    Types of fingerprints:
        1) RDKit (default)
        2) Morgan

    Types of similarity score:
        1) Tanimoto (default)
        2) Dice

    Parameters
    ----------
    crit_tani : float
        The Tanimoto similarity score threshold.
    increasing_tani : bool
        Whether or not to only keep compounds whos Tanimoto score is higher than its
        parent.
    fingerprint_method : str
        Function to calculate fingerprint with. Accepts a SMILES, by default uses
        RDKFingerprints
    fingerprint_args : str
        Arguments for the fingerprint_method, by default None.
    similarity_method : str
        Function to calculate similarity with. Accepts two fingerprints
        and returns a similarity score.

    Attributes
    ----------
    crit_tani : float
        The Tanimoto similarity score threshold.
    increasing_tani : bool
        Whether or not to only keep compounds whos Tanimoto score is higher than its
        parent.
    fingerprint_method : str
        Fingerprint calculation method.
    fingerprint_args : dict
        Arguments for fingerprint method.
    similairty_method : str
        Simlarity calulation method.
    target_fps : List
        List of fingerprints.
    """

    def __init__(
        self,
        crit_tani: float,
        increasing_tani: bool,
        fingerprint_method: str = "RDKit",
        fingerprint_args: dict = None,
        similarity_method: str = "Tanimoto",
    ) -> None:
        self._filter_name = "Tanimoto Cutoff"
        self.crit_tani = crit_tani
        self.increasing_tani = increasing_tani

        self.fingerprint_method = fingerprint_method
        self.fingerprint_args = fingerprint_args or dict()
        self.similarity_method = similarity_method

        self.target_fps = []

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _set_target_fps(self, pickaxe: Pickaxe):
        self.target_fps = []
        for smiles in pickaxe.target_smiles:
            mol = MolFromSmiles(smiles)
            if self.fingerprint_method == "Morgan":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, **self.fingerprint_args)
            else:
                fp = RDKFingerprint(mol, **self.fingerprint_args)

            self.target_fps.append(fp)

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a Tanimoto similarity score to a target
        compound greater than or equal to the crit_tani, for expansion.
        """
        self._set_target_fps(pickaxe)
        if not self.target_fps:
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
            cpd_info, self.target_fps, processes, this_gen_crit_tani
        )

        # Process filtering results
        cpds_remove_set = set()
        for c_id, current_tani in cpd_filters:
            # Check if tani is increasing
            if self.increasing_tani:
                if current_tani >= pickaxe.compounds[c_id]["last_tani"]:
                    pickaxe.compounds[c_id]["last_tani"] = current_tani
                else:
                    pickaxe.compounds[c_id]["Expand"] = False
                    cpds_remove_set.add(c_id)
                    continue

            if current_tani < this_gen_crit_tani:
                pickaxe.compounds[c_id]["Expand"] = False
                cpds_remove_set.add(c_id)

        return cpds_remove_set

    def _filter_by_tani_helper(
        self,
        compounds_info: List[Tuple[str, str]],
        target_fps: List[AllChem.RDKFingerprint],
        processes: int,
        this_crit_tani: float,
    ) -> List[Tuple[str, float]]:
        def print_progress(done: int, total: int, section: str) -> None:
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

    def _compare_target_fps(
        self,
        target_fps: List[AllChem.RDKFingerprint],
        this_crit_tani: float,
        compound_info: Tuple[str, str],
    ) -> Tuple[str, float]:
        # do finger print loop here
        """
        Helper function to allow parallel computation of Tanimoto filtering.
        Works with _filter_by_tani_helper.

        Returns cpd_id if a the compound is similar enough to a target.
        """
        # Generate the fingerprint of a compound and compare to the fingerprints
        # of the targets
        def fingerprint(fingerprint_method, keyword_dict, smi):
            mol = AllChem.MolFromSmiles(smi)
            if fingerprint_method == "Morgan":
                return AllChem.GetMorganFingerprintAsBitVect(mol, **keyword_dict)
            else:
                return RDKFingerprint(mol, **keyword_dict)

        def similarity(similarity_method, fp1, fp2):
            if similarity_method == "dice":
                return FingerprintSimilarity(
                    fp1, fp2, metric=DataStructs.DiceSimilarity
                )
            else:
                return FingerprintSimilarity(fp1, fp2)

        fingerprint_partial = partial(
            fingerprint, self.fingerprint_method, self.fingerprint_args
        )

        try:
            fp1 = fingerprint_partial(compound_info[1])
            max_tani = 0
            for fp2 in target_fps:
                tani = similarity(self.similarity_method, fp1, fp2)
                if tani >= this_crit_tani:
                    return (compound_info[0], tani)
                elif tani >= max_tani:
                    max_tani = tani
            return (compound_info[0], max_tani)
            # TODO what except to use here?
        except:  # noqa
            return (compound_info[0], -1)

    def preprint(self, pickaxe: Pickaxe) -> None:
        if type(self.crit_tani) in [list, tuple]:
            print_tani = self.crit_tani[pickaxe.generation]
        else:
            print_tani = self.crit_tani
        print(f"Filtering out compounds with maximum Tanimoto match < {print_tani}")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
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

    def __init__(self, crit_mcs: float) -> None:
        self._filter_name = "MCS Cutoff"
        self.crit_mcs = crit_mcs

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a Tanimoto similarity score to a target
        compound greater than or equal to the crit_tani, for expansion.
        """

        if not self.target_fps:
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
        self,
        compounds_info: List[Tuple[str, str]],
        target_smiles: List[str],
        processes: int,
        this_crit_mcs: float,
        retro: bool = False,
    ) -> List[Tuple[str, float]]:
        def print_progress(done: int, total: int, section: str) -> None:
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

    def _compare_target_mcs(
        self,
        target_smiles: List[str],
        retro: bool,
        compound_info: Tuple[str, str],
        this_crit_mcs: float,
    ) -> Tuple[str, float]:
        """Compare target MCS.

        Helper function to allow parallel computation of MCS filtering.
        Works with _filter_by_mcs_helper.

        Returns cpd_id if a the compound is similar enough to a target.

        """

        def get_mcs_overlap(mol, target_mol) -> float:
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

    def preprint(self, pickaxe: Pickaxe) -> None:
        if type(self.crit_mcs) in [list, tuple]:
            if len(self.crit_mcs) - 1 < pickaxe.generation:
                crit_mcs = self.crit_mcs[-1]
            else:
                crit_mcs = self.crit_mcs[pickaxe.generation]
        else:
            crit_mcs = self.crit_mcs
        print(f"Filtering out compounds with maximum MCS match < {crit_mcs}")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
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


class MWFilter(Filter):
    """A filter that removes compounds not within a MW range.

    This filter specifies a minimum and maximum molecular weight to create a range.
    Specifying None for either value will yield an unbounded MW range on that end.

    For example, specifying min_MW = None and max_MW = 1000 will give compounds less
    than or equal to 1000 g/mol.

    Parameters
    ----------
    min_MW : Union[float, None]
        Minimum MW in g/mol, by default None.
    max_MW : Union[float, None]
        Maximum MW in g/mol, by default None.

    Attributes
    ----------
    min_MW : Union[float, None]
        Minimum MW in g/mol.
    max_MW : Union[float, None]
        Maximum MW in g/mol.
    """

    def __init__(
        self,
        min_MW: Union[float, None] = None,
        max_MW: Union[float, None] = None,
    ) -> None:
        self._filter_name = "Molecular Weight"

        self.min_MW = min_MW or -1
        self.max_MW = max_MW

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Check the compounds against the MW constraints and return
        compounds to filter.
        """

        def MW_is_good(cpd):
            cpd_MW = periodictable.formula(cpd["Formula"]).mass
            return self.min_MW < cpd_MW and cpd_MW < self.max_MW

        def is_target(cpd, pickaxe):
            for t_id in pickaxe.targets:
                if "C" + t_id[1:] == cpd["_id"]:
                    return True
            return False

        cpds_remove_set = set()

        print(
            f"Filtering Generation {pickaxe.generation} "
            f"with {self.min_MW} < MW < {self.max_MW}."
        )

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:
                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    if not MW_is_good(cpd):
                        cpds_remove_set.add(cpd["_id"])
                else:
                    if is_target(cpd, pickaxe):
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                    else:
                        if not MW_is_good(cpd):
                            cpds_remove_set.add(cpd["_id"])

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set


class AtomicCompositionFilter(Filter):
    """Filter that removes compounds not within specific atomic composition ranges.

    This filter checks to see if the atomic composition of a compound is
    within a specified range. As an example, to only keep compounds with carbon between 4-7 and
    oxygen between 0-4 the following input would be used:

    atomic_composition = {"C": [4, 7], "O": [0, 4]}

    Parameters
    ----------
    atomic_composition : Dict
        A dictionary containing ranges for elemental composition. Of form
        {"Atom Symbol": [min, max]}, by default None.

    Attributes
    ----------
    atomic_composition : Dict
        A dictionary containing ranges for elemental composition.
    """

    def __init__(
        self,
        atomic_composition_constraints: Dict = None,
    ) -> None:
        self._filter_name = "Atomic Composition"

        self.atomic_composition_constraints = atomic_composition_constraints

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _choose_cpds_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Check the compounds against the atomic composition constraints and return
        compounds to filter.
        """

        def composition_is_good(cpd):
            atom_count = cpd["atom_count"]
            for atom in atom_count:
                atom_range = self.atomic_composition_constraints.get(atom)
                if atom_range:
                    atom_min = atom_range[0] or 0
                    atom_max = atom_range[1] or 10 ** 5

                    if not (
                        atom_min <= atom_count[atom] and atom_count[atom] <= atom_max
                    ):
                        return False

            return True

        def is_target(cpd, pickaxe):
            for t_id in pickaxe.targets:
                if "C" + t_id[1:] == cpd["_id"]:
                    return True
            return False

        cpds_remove_set = set()

        print(
            f"Filtering Generation {pickaxe.generation} "
            f"with atomic composition {self.atomic_composition_constraints}."
        )

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:
                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    if not composition_is_good(cpd):
                        cpds_remove_set.add(cpd["_id"])
                else:
                    if is_target(cpd, pickaxe):
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                    else:
                        if not composition_is_good(cpd):
                            cpds_remove_set.add(cpd["_id"])

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set


#
#
# End filters
###############################################################################
