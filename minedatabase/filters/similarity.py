from typing import Callable, List

from minedatabase.filters.base_filter import Filter
from minedatabase.pickaxe import Pickaxe


"""Definitions of filters for pickaxe.

Use this module to define your own filter classes. All filter classes must
subclass Filter. See Filter docstring for more information.

To use any filter, import it in pickaxe_run.py, initialize it, and append it
to the .filters property of your pickaxe object.
"""

import copy
import multiprocessing
import time
from functools import partial
from typing import Callable, List, Set, Tuple

import numpy as np
import pandas as pd
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from rdkit import DataStructs
from rdkit.Chem import AllChem, RDKFingerprint
from rdkit.Chem import rdFMCS as mcs
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.DataStructs import FingerprintSimilarity
from scipy.stats import rv_discrete

from minedatabase.utils import Chunks


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")

###############################################################################
# Tanimoto Sampling Filter


class SimilaritySamplingFilter(Filter):
    """Filter that samples randomly from weighted similarity scores.

    SimilaritySamplingFilter takes a distribution of similarity scores and uses
    inverse CDF sampling to select N compounds for further expansion. Each compound
    is assigned a similarity score that corresponds to the maximum similarity score
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
        Function to weight the similarity score with.
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
        Function to weight the similarity score with.
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
        self._filter_name = "Similarity Sampling Filter"
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
        for smiles in pickaxe.target_smiles:
            mol = MolFromSmiles(smiles)
            if self.fingerprint_method == "Morgan":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, **self.fingerprint_args)
            else:
                fp = RDKFingerprint(mol)

            self.target_fps.append(fp)

    def _pre_print(self) -> None:
        """Print before filtering."""
        print(
            (
                f"Sampling {self.sample_size} Compounds Based on a "
                f"Weighted Similarity Distribution"
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
                f"Similarity Sampling of generation {pickaxe.generation}"
                f"--took {time.time() - time_sample}s.\n"
            )
        )

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
        """
        Samples N compounds to expand based on the weighted Similarity distribution.

        Parameters
        ----------
        pickaxe : Pickaxe
            Pickaxe object to filter
        processes : int
            Number of processes to use.
        """

        print(f"Filtering Generation {pickaxe.generation}" " via Similarity Sampling.")

        if not self.target_fps:
            self._set_target_fps(pickaxe)
        if not self.target_fps:
            print("No targets to filter for. Can't expand.")
            return [], []

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
            min_S=0.15,
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

        return cpds_remove_set, []

    def _sample_by_similarity(
        self,
        mol_info: List[Tuple[str, str]],
        t_fp: List[AllChem.RDKFingerprint],
        n_cpds: int = None,
        min_S: float = 0.05,
        weighting: Callable = None,
        max_iter: int = None,
        processes: int = 1,
    ) -> List[str]:
        """Sample compounds by weighted similarity coefficient.

        Use inverse cumulative distrbution function (CDF) sampling to select
        compounds based on a weighted similarity coefficient distribution.

        Parameters
        ----------
        mol_info : List[Tuple[str, str]]
            A list consisting of (compound_id, SMILES).
        t_fp : List[RDKFingerprint]
            Target fingerprints to compare compounds to.
        n_cpds : int, optional
            Number of compounds to select for sampling, by default None.
        min_S : float, optional
            Minimum similarity similarity to be considered for sampling, by default 0.05.
        weighting : Callable, optional
            Function that accepts a similarity score and returns
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
        df = self.gen_df_from_similarity(
            mol_info, t_fp, min_S=min_S, processes=processes
        )
        if len(df) <= n_cpds:
            ids = set(df["_id"])
            print(
                f"-- After filtering by minimum similarity ({min_S}) "
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
        nCDF = 0

        collisions_before_recalc = n_cpds / 4
        collisions = 0

        while len(chosen_ids) != n_cpds:
            # Collisions merit reset
            if collisions == collisions_before_recalc:
                collisions = 0
                nCDF += 1
                rv, ids = self._gen_rv_from_df(df, chosen_ids, weighting=weighting)

            sampled_idx = ids.iloc[rv.rvs(size=1)[0]]
            if sampled_idx in chosen_ids:
                # collision
                collisions += 1
            else:
                chosen_ids.add(sampled_idx)
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
            Function to weight the Similarity distribution by, by default None.

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

    def gen_df_from_similarity(
        self,
        mol_info: List[Tuple[str, str]],
        t_fp: List[AllChem.RDKFingerprint],
        min_S: float = 0.05,
        processes: int = 1,
    ) -> pd.DataFrame:
        """Generate a dataframe from similarity

        Parameters
        ----------
        mol_info : List[Tuple[str, str]]
            A list consisting of (compound_id, SMILES).
        t_fp : List[RDKFingerprint]
            Target fingerprints to compare compounds to.
        min_S : float, optional
            Minimum similarity similarity to be considered for sampling, by default 0.05.
        processes : int, optional
            Number of processes to use, by default 1.
        """

        then = time.time()
        print("-- Calculating Fingerprints and Similarity Values.")
        # target fingerprint dataframe
        t_df = pd.DataFrame(t_fp, columns=["fp"])

        # Calculate similarity for each compound and drop T < min_S
        partial_T_calc = partial(
            _calc_max_S,
            t_df,
            min_S,
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
        print(f"-- Completed Similarity Calculation in {time.time() - then}")

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


def _calc_max_S(
    t_df: pd.DataFrame,
    min_S: float,
    fingerprint_method: str,
    fingerprint_args: str,
    similarity_method: str,
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Calculate maximum similarity.

    Generate the similarity to use to generate the PMF to sample from.
    For each compound a list of similarity values are obtained by a generated
    compound to every target compound and the max is taken.

    Parameters
    ----------
    t_df : pd.Dataframe
        Dataframe containing the target fingerprints.
    min_S : float
        The minimum similarity similarity score needed to consider a compound.
    fingerprint_method: str
        Which fingerprint method to use, suppoorts RDKit and Morgan, by default RDKit.
    fingerprint_args: dict
        Keyword arguments to pass to fingerprint function, by default empty dict.
    similarity_method: str
        Which similarity method to use. Supports Tanimoto and Dice.
    df : pd.DataFrame
        Dataframe to calculate the max similarity for.

    Returns
    -------
    df : pd.Dataframe
        New dataframe with max similarity values calculated.
    """

    def fingerprint(fingerprint_method, keyword_dict, smi):
        mol = AllChem.MolFromSmiles(smi)
        if fingerprint_method == "Morgan":
            return AllChem.GetMorganFingerprintAsBitVect(mol, **keyword_dict)
        else:
            return RDKFingerprint(mol)

    def calc_similarity(similarity_method, fp1, fp2):
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
            t_df["fp"].map(lambda x: calc_similarity(similarity_method, x, fp))
        )
    # Filter out low similarity
    df = df[df["T"] > min_S]

    return df


# End similarity Sampling Filter
###############################################################################

###############################################################################
# Cutoff filters -- e.g. similarity, MCS metric, and molecular weight


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
    crit_similarity : float
        The similarity similarity score threshold.
    increasing_similarity : bool
        Whether or not to only keep compounds whos similarity score is higher than its
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
    crit_similarity : float
        The similarity similarity score threshold.
    increasing_similarity : bool
        Whether or not to only keep compounds whos similarity score is higher than its
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
        crit_similarity: float,
        increasing_similarity: bool,
        fingerprint_method: str = "RDKit",
        fingerprint_args: dict = None,
        similarity_method: str = "Tanimoto",
    ) -> None:
        self._filter_name = "Similarity Cutoff"
        self.crit_similarity = crit_similarity
        self.increasing_similarity = increasing_similarity

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
                fingerprint_args = self.fingerprint_args["radius"]
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, fingerprint_args)
            else:
                fp = RDKFingerprint(mol)

            self.target_fps.append(fp)

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a similarity score to a target
        compound greater than or equal to the crit_similarity, for expansion.
        """
        if not self.target_fps:
            self._set_target_fps(pickaxe)
        if not self.target_fps:
            print("No targets to filter for. Can't expand.")
            return [], []

        # Set up variables required for filtering
        # similarity Threshold
        if type(self.crit_similarity) in [list, tuple]:
            if len(self.crit_similarity) - 1 < pickaxe.generation:
                crit_similarity = self.crit_similarity[-1]
            else:
                crit_similarity = self.crit_similarity[pickaxe.generation]
        else:
            crit_similarity = self.crit_similarity

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
            f"with similarity > {crit_similarity}."
        )
        # Get input to filter code, c_id and smiles (to be
        # turned into fingerprint)
        cpd_info = [(cpd["_id"], cpd["SMILES"]) for cpd in compounds_to_check]
        if type(self.crit_similarity) in [list, tuple]:
            if len(self.crit_similarity) - 1 < pickaxe.generation:
                this_gen_crit_similarity = self.crit_similarity[-1]
            else:
                this_gen_crit_similarity = self.crit_similarity[pickaxe.generation]
        else:
            this_gen_crit_similarity = self.crit_similarity
        cpd_filters = self._filter_by_similarity_helper(
            cpd_info, self.target_fps, processes, this_gen_crit_similarity
        )

        # Process filtering results
        cpds_remove_set = set()
        for c_id, current_similarity in cpd_filters:
            # Check if similarity is increasing
            if self.increasing_similarity:
                if current_similarity >= pickaxe.compounds[c_id]["last_similarity"]:
                    pickaxe.compounds[c_id]["last_similarity"] = current_similarity
                else:
                    pickaxe.compounds[c_id]["Expand"] = False
                    cpds_remove_set.add(c_id)
                    continue

            if current_similarity < this_gen_crit_similarity:
                pickaxe.compounds[c_id]["Expand"] = False
                cpds_remove_set.add(c_id)

        return cpds_remove_set, []

    def _filter_by_similarity_helper(
        self,
        compounds_info: List[Tuple[str, str]],
        target_fps: List[AllChem.RDKFingerprint],
        processes: int,
        this_crit_similarity: float,
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
            self._compare_target_fps, target_fps, this_crit_similarity
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
                print_progress(i, len(compounds_info), "Similarity filter progress:")

        else:
            for i, cpd in enumerate(compounds_info):
                res = compare_target_fps_partial(cpd)
                if res:
                    cpds_to_filter.append(res)
                print_progress(i, len(compounds_info), "Similarity filter progress:")
        print("Similarity filter progress: 100 percent complete")

        return cpds_to_filter

    def _compare_target_fps(
        self,
        target_fps: List[AllChem.RDKFingerprint],
        this_crit_similarity: float,
        compound_info: Tuple[str, str],
    ) -> Tuple[str, float]:
        # do finger print loop here
        """
        Helper function to allow parallel computation of Similarity filtering.
        Works with _filter_by_similarity_helper.

        Returns cpd_id if a the compound is similar enough to a target.
        """
        # Generate the fingerprint of a compound and compare to the fingerprints
        # of the targets
        def fingerprint(fingerprint_method, keyword_dict, smi):
            mol = AllChem.MolFromSmiles(smi)
            if fingerprint_method == "Morgan":
                return AllChem.GetMorganFingerprintAsBitVect(mol, **keyword_dict)
            else:
                return RDKFingerprint(mol)

        def calc_similarity(similarity_method, fp1, fp2):
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
            max_similarity = 0
            for fp2 in target_fps:
                similarity = calc_similarity(self.similarity_method, fp1, fp2)
                if similarity >= this_crit_similarity:
                    return (compound_info[0], similarity)
                elif similarity >= max_similarity:
                    max_similarity = similarity
            return (compound_info[0], max_similarity)
            # TODO what except to use here?
        except:  # noqa
            return (compound_info[0], -1)

    def preprint(self, pickaxe: Pickaxe) -> None:
        if type(self.crit_similarity) in [list, tuple]:
            print_similarity = self.crit_similarity[pickaxe.generation]
        else:
            print_similarity = self.crit_similarity
        print(
            f"Filtering out compounds with maximum Similarity match < {print_similarity}"
        )

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
        if type(self.crit_similarity) in [list, tuple]:
            if len(self.crit_similarity) - 1 < pickaxe.generation:
                print_similarity = self.crit_similarity[-1]
            else:
                print_similarity = self.crit_similarity[pickaxe.generation]
        else:
            print_similarity = self.crit_similarity
        print(
            (
                f"{n_filtered} of {n_total} compounds selected after "
                f"Similarity filtering of generation {pickaxe.generation} "
                f"at cutoff of {print_similarity}. "
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

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a similarity score to a target
        compound greater than or equal to the crit_similarity, for expansion.
        """

        if not pickaxe.target_smiles:
            print("No targets to filter for. Can't expand.")
            return [], []

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

        return cpds_remove_set, []

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


if __name__ == "__main__":
    pass
