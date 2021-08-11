import multiprocessing
import time
from functools import partial
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from mordred import Calculator, descriptors
from rdkit.Chem import AddHs, CanonSmiles, MolFromSmiles
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdmolfiles import MolFromSmiles
from sklearn.ensemble import RandomForestRegressor

from minedatabase.filters.base_filter import Filter
from minedatabase.metabolomics import MetabolomicsDataset, Peak
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import neutralise_charges


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")


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
                charge=row["Polarity"],
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

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
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
        # Function to determine if a reaction has a match in it
        def reaction_makes_match(pickaxe, rxn_id, matches):
            products = pickaxe.reactions[rxn_id]["Products"]
            for _, c_id in products:
                if c_id in matches:
                    return True
            return False

        if pickaxe.generation == 0:
            return None, None

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

        possible_ranges = (
            self.metabolomics_dataset.possible_ranges["+"]
            + self.metabolomics_dataset.possible_ranges["-"]
        )

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

        rxns_remove_set = set()
        # Remove any reaction that does not produce a t
        for cpd_id in cpds_remove_set:
            for rxn_id in pickaxe.compounds[cpd_id].get("Product_of", []):
                # Reaction doesn't make a met match
                if not reaction_makes_match(pickaxe, rxn_id, ids):
                    rxns_remove_set.update([rxn_id])

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set, rxns_remove_set

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


if __name__ == "__main__":
    pass
