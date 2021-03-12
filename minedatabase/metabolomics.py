"""Provides functionality to interact with metabolomics datasets and search
MINE databases for metabolomics hits."""
from __future__ import annotations

import math
import os
import re
import time
import xml.etree.ElementTree as ET
from ast import literal_eval
from typing import Callable, Generator, List, Optional, Set, Tuple

import numpy as np
import pymongo

import minedatabase
from minedatabase.databases import MINE
from minedatabase.utils import score_compounds


MINEDB_DIR = os.path.dirname(minedatabase.__file__)


class MetabolomicsDataset:
    """A class containing all the information for a metabolomics data set."""

    # pylint: disable=redefined-outer-name
    def __init__(
        self,
        name: str,
        adducts: List[str] = None,
        known_peaks: List[Peak] = None,
        unknown_peaks: List[Peak] = None,
        native_set: Set[str] = set(),
        ppm: bool = False,
        tolerance: float = 0.001,
        halogens: bool = False,
        verbose: bool = False
    ):
        """Metabolomics Dataset initialization.

        Parameters
        ----------
        name : str
            Name of metabolomics dataset
        adducts : [str], optional
            List of adduct names, e.g. ["[M-H]-", "[M+H]+"] (defaults to all)
            See minedatabase/data/adducts/{Negative/Positive} Adducts full.txt
            for a full list of possible adducts
        known_peaks : [Peak], optional
            List of Peak objects annotated with ID and associated data, by
            default None
        unknown_peaks : [Peak], optional
            List of Peak objects annotated with associated data, by default
            None
        native_set : set(str), optional
            Set of compound IDs native to the organism that generated the
            dataset (e.g. IDs from model), by default set()
        ppm : bool, optional
            If True, tolerance is set in parts per million, and if false
            (default), tolerance is set in Daltons
        tolerance : float, optional
            Mass tolerance for hits, by default 0.001
        halogens : bool, optional
            Filters out compounds containing halogens if True, by default False
        verbose : bool, optional
            Prints more info to stdout if True, by default False
        """
        # Load adducts
        dtype = "np.str, np.float, np.float"
        pos_fp = os.path.join(MINEDB_DIR, "data/adducts/Positive Adducts full.txt")
        neg_fp = os.path.join(MINEDB_DIR, "data/adducts/Negative Adducts full.txt")
        all_pos_adducts = self._read_adduct_file(pos_fp)
        all_neg_adducts = self._read_adduct_file(neg_fp)

        if adducts:
            self.pos_adducts = list(filter(lambda x: x[0] in adducts, all_pos_adducts))
            self.neg_adducts = list(filter(lambda x: x[0] in adducts, all_neg_adducts))
        else:
            self.pos_adducts = all_pos_adducts
            self.neg_adducts = all_neg_adducts

        # MongoDB projection for compound search
        self.hit_projection = {
            "Formula": 1,
            "MINE_id": 1,
            "NP_likeness": 1,
            "Names": 1,
            "SMILES": 1,
            "Inchikey": 1,
            "Generation": 1,
            "Pos_CFM_spectra": 1,
            "Neg_CFM_spectra": 1,
            "Sources": 1,
        }

        # Load peak data and initialize other attributes
        self.name = name
        self.known_peaks = known_peaks
        self.unknown_peaks = unknown_peaks
        self.native_set = native_set
        self.ppm = ppm
        self.tolerance = tolerance
        self.halogens = halogens
        self.verbose = verbose

        self.total_formulas = 0
        self.total_hits = 0
        self.matched_peaks = 0
        self.possible_masses = []
        self.possible_ranges = []

    def __str__(self) -> str:
        """Give string representation.

        Returns
        -------
        self.name : str
            Name of the metabolomics dataset
        """
        return self.name

    def _read_adduct_file(self, filepath) -> List[Tuple]:
        """Read specified adduct file.

        Parameters
        ----------
        filepath : str
            Path to adduct file

        Returns
        -------
        adducts : List[Tuple]
            A list of (str, float, float) tuples of form
            ('adduct name', m/z multiplier, adduct mass change)
        """
        adducts = []
        with open(filepath, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                adduct = line.strip().split('\t')
                adduct[0] = adduct[0].strip()
                adduct[1] = float(adduct[1])
                adduct[2] = float(adduct[2])
                adducts.append(tuple(adduct))
        return adducts

    def enumerate_possible_masses(self, tolerance: float):
        """Generate all possible masses from unknown peaks and list of
        adducts. Saves these mass ranges to self.possible_ranges.

        Parameters
        ----------
        tolerance : float
            Mass tolerance in Daltons
        """
        possible_masses = []
        possible_ranges = []
        for peak in self.unknown_peaks:
            for adduct in self.pos_adducts + self.neg_adducts:
                possible_mass = (peak.mz - adduct[2]) / adduct[1]
                possible_masses.append(possible_mass)
                possible_ranges.append(
                    (
                        possible_mass - tolerance,
                        possible_mass + tolerance,
                        peak.name,
                        adduct[1],
                    )
                )

        self.possible_masses = np.array(set(possible_masses))
        self.possible_ranges = possible_ranges

    def get_rt(self, peak_id: str) -> Optional[float]:
        """Return retention time for peak with given ID. If not found, returns
        None.

        Parameters
        ----------
        peak_id : str
            ID of peak as listed in dataset

        Returns
        -------
        rt : float, optional
            Retention time of peak with given ID, None if not found
        """
        rt = None
        for peak in self.unknown_peaks + self.known_peaks:
            if peak_id == peak.name:
                rt = peak.r_time
                break
        return rt

    def find_db_hits(self, peak: Peak, db: MINE, adducts: List[Tuple[str, float, float]]):
        """This function searches the database for matches of a peak given
        adducts and updates the peak object with that information.

        Parameters
        ----------
        peak : Peak
            Peak object to query against MINE compound database
        db : MINE
            MINE database to query
        adducts : List[Tuple[str, float, float]]
            List of adducts. Each adduct contains three values in a tuple:
            (adduct name, mass multiplier, ion mass)
        """
        # find nominal mass for a given m/z for each adduct and the max and
        # min values for db
        potential_masses = (peak.mz - adducts["f2"]) / adducts["f1"]
        if self.ppm:
            precision = (self.tolerance / 100000.0) * potential_masses
        else:
            precision = self.tolerance * 0.001  # convert to mDa
        upper_bounds = potential_masses + precision
        lower_bounds = potential_masses - precision

        # search database for hits in the each adducts mass range that have no
        # innate charge.
        for i, adduct in enumerate(adducts):
            # build the query by adding the optional terms
            query_terms = [
                {"Mass": {"$gte": float(lower_bounds[i])}},
                {"Mass": {"$lte": float(upper_bounds[i])}},
                {"Charge": 0},
            ]
            if adduct["f0"] == "[M]+":
                query_terms[2] = {"Charge": 1}
            for compound in db.compounds.find(
                {"$and": query_terms}, self.hit_projection
            ):
                # Filters out halogens if the flag is enabled by moving to the
                # next compound before the current compound is counted or
                # stored.
                if not self.halogens:
                    if re.search("F[^e]|Cl|Br", compound["Formula"]):
                        continue

                # update the total hits for the peak and make a note if the
                # compound is in the native_set
                peak.total_hits += 1

                if compound["_id"] in self.native_set:
                    peak.native_hit = True
                    compound["native_hit"] = True

                if compound["Generation"] < peak.min_steps:
                    peak.min_steps = compound["Generation"]

                peak.formulas.add(compound["Formula"])
                compound["adduct"] = adduct["f0"]
                compound["peak_name"] = peak.name
                peak.isomers.append(compound)

    def annotate_peaks(self, db: MINE):
        """This function iterates through the unknown peaks in the dataset and
        searches the database for compounds that match a peak m/z given the
        adducts permitted. Statistics on the annotated data set are printed.

        Parameters
        ----------
        db : MINE
            MINE database
        """
        for i, peak in enumerate(self.unknown_peaks):

            positive = (
                peak.charge == "+"
                or peak.charge == "Positive"
                or (peak.charge and isinstance(peak.charge, bool))
            )
            negative = (
                peak.charge == "-"
                or peak.charge == "Negative"
                or (not peak.charge and isinstance(peak.charge, bool))
            )

            if positive:
                self.find_db_hits(peak, db, self.pos_adducts)
            elif negative:
                self.find_db_hits(peak, db, self.neg_adducts)
            else:
                raise ValueError(
                    "Invalid compound charge specification. "
                    'Please use "+" or "Positive" for '
                    'positive ions and "-" or "Negative" for '
                    f"negative ions. (charge = {peak.charge})"
                )

            if peak.total_hits > 0:
                self.matched_peaks += 1
                self.total_hits += peak.total_hits
                self.total_formulas += len(peak.formulas)
            if self.verbose:
                pct_done = int(float(i) / float(len(self.unknown_peaks)) * 100)
                print(f"{pct_done} percent of peaks processed")


# Scoring functions appear before the Peak class because dot_product method is
# default object for Peak.score_isomers

def dot_product(x: List[tuple], y: List[tuple], epsilon: float = 0.01) -> float:
    """Calculate the dot product of two spectra, allowing for some variability
    in mass-to-charge ratios

    Parameters
    ----------
    x : List[tuple]
        First spectra m/z values
    y : List[tuple]
        Second spectra m/z values
    epsilon : float, optional
        Mass tolerance in Daltons, by default 0.01

    Returns
    -------
    dot_prod : float
        Dot product of x and y
    """
    z = 0
    n_v1 = 0
    n_v2 = 0

    for int1, int2 in _approximate_matches(x, y, epsilon):
        z += int1 * int2
        n_v1 += int1 * int1
        n_v2 += int2 * int2

    dot_prod = z / (math.sqrt(n_v1) * math.sqrt(n_v2))
    return dot_prod


def jaccard(x: List[tuple], y: List[tuple], epsilon: float = 0.01) -> float:
    """Calculate the Jaccard Index of two spectra, allowing for some
    variability in mass-to-charge ratios

    Parameters
    ----------
    x : List[tuple]
        First spectra m/z values
    y : List[tuple]
        Second spectra m/z values
    epsilon : float, optional
        Mass tolerance in Daltons, by default 0.01

    Returns
    -------
    jaccard_index : float
        Jaccard Index of x and y
    """
    intersect = 0

    for val1, val2 in _approximate_matches(x, y, epsilon):
        if val1 and val2:
            intersect += 1

    jaccard_index = intersect / float((len(x) + len(y) - intersect))
    return jaccard_index


def _approximate_matches(list1: List[tuple], list2: List[tuple],
                         epsilon: float = 0.01) -> Generator:
    """Takes two list of tuples and searches for matches of tuples first value
    within the supplied epsilon. Emits tuples with the tuples second values
    where found. if a value in one dist does not match the other list, it is
    emitted alone but with a 0 as the other value.

    Parameters
    ----------
    list1 : list
        First list of tuples
    list2 : list
        Second list of tuples
    epsilon : float, optional
        Maximum difference, by default 0.01

    Yields
    -------
    Generator
        Generator that yields found matches
    """
    list1.sort()
    list2.sort()
    list1_index = 0
    list2_index = 0

    while list1_index < len(list1) or list2_index < len(list2):
        if list1_index == len(list1):
            yield (0, list2[list2_index][1])
            list2_index += 1
            continue
        if list2_index == len(list2):
            yield (list1[list1_index][1], 0)
            list1_index += 1
            continue

        list1_element = list1[list1_index][0]
        list2_element = list2[list2_index][0]

        difference = abs(list1_element - list2_element)

        if difference < epsilon:
            yield (list1[list1_index][1], list2[list2_index][1])
            list1_index += 1
            list2_index += 1
        elif list1_element < list2_element:
            yield (list1[list1_index][1], 0)
            list1_index += 1
        elif list2_element < list1_element:
            yield (0, list2[list2_index][1])
            list2_index += 1


class Peak:
    """A class holding information about an metabolomics peak"""

    def __init__(
        self,
        name : str,
        r_time : float,
        mz : float,
        charge : str,
        inchi_key : str = None,
        ms2 : List[(float, float)] = None
    ):
        """Peak object which contains peak metadata as well as mass, retention
        time, spectra, and any MINE database hits.

        Parameters
        ----------
        name : str
            Name or ID of the peak
        r_time : float
            Retention time of the peak
        mz : float
            Mass-to-charge ratio (m/z) of the peak
        charge : str
            Charge of the peak, "+" or "-"
        inchi_key : str, optional
            InChI key of the peak, if already identified, by default None
        ms2 : List[float], optional
            MS2 spectra m/z values for this peak, by default None
        """
        self.name = name
        self.r_time = float(r_time)
        self.mz = float(mz)
        self.charge = charge
        self.inchi_key = inchi_key
        self.ms2peaks = ms2

        self.isomers = []
        self.formulas = set()
        self.total_hits = 0
        self.native_hit = False
        self.min_steps = 99

    def __str__(self) -> str:
        """String representation of the peak.

        Returns
        -------
        str
            Name of the peak
        """
        return self.name

    def score_isomers(
        self,
        metric: Callable[[list, list], float] = dot_product,
        energy_level: int = 20,
        tolerance: float = 0.005,
    ):
        """Scores and sorts isomers based on mass spectra data.

        Calculates the cosign similarity score between the provided ms2 peak
        list and pre-calculated CFM-spectra and sorts the isomer list
        according to this metric.

        Parameters
        ----------
        metric : function, optional
            The scoring metric to use for the spectra. Function must accept 2
            lists of (mz, intensity) tuples and return a score, by default dot_product
        energy_level : int, optional
            The Fragmentation energy level to use. May be 10,
            20 or 40., by default 20
        tolerance : float, optional
            The precision to use for matching m/z in mDa, by default 0.005

        Raises
        ------
        ValueError
            Empty ms2 peak
        """
        if not self.ms2peaks:
            raise ValueError("The ms2 peak list is empty")
        if self.charge:
            spec_key = "Pos_CFM_spectra"
        else:
            spec_key = "Neg_CFM_spectra"

        for i, hit in enumerate(self.isomers):
            if spec_key in hit:
                hit_spec = hit[spec_key][f"{energy_level} V"]
                score = metric(self.ms2peaks, hit_spec, epsilon=tolerance)
                rounded_score = round(score * 1000)
                self.isomers[i]["Spectral_score"] = rounded_score
                del hit[spec_key]
            else:
                self.isomers[i]["Spectral_score"] = None
        self.isomers.sort(key=lambda x: x["Spectral_score"], reverse=True)



def get_KEGG_comps(
    db: MINE,
    kegg_db: pymongo.database.Database,
    model_ids: List[str]) -> set:
    """Get KEGG IDs from KEGG MINE database for compounds in model(s).

    Parameters
    ----------
    db : MINE
        MINE Mongo database
    kegg_db : pymongo.database.Database
        Mongo database with annotated organism metabolomes from KEGG
    model_ids : List[str]
        List of organism identifiers from KEGG

    Returns
    -------
    set
        MINE IDs of compounds that are linked to a KEGG ID in at least one of
        the organisms in model_ids
    """
    kegg_ids, _ids = set(), set()
    for model_id in model_ids:
        comp_ids = kegg_db.models.find_one({"_id": model_id}, {"Compounds": 1})[
            "Compounds"
        ]
        for comp_id in comp_ids:
            kegg_ids.add(comp_id)
    for kegg_id in kegg_ids:
        for comp in db.compounds.find({"DB_links.KEGG": kegg_id}, {"_id": 1}):
            _ids.add(comp["_id"])
    return _ids


def read_adduct_names(filepath: str) -> List[str]:
    """Read adduct names from text file at specified path into a list.

    Parameters
    ----------
    filepath : str
        Path to adduct text file.

    Returns
    -------
    adducts : list
        Names of adducts in text file.

    Notes
    -----
    Not used in this codebase but used by MINE-Server to validate adduct input.
    """

    with open(filepath) as infile:
        adducts = [line.split(" \t")[0] for line in infile if not line[0] == "#"]

    return adducts


def read_mgf(input_string: str, charge: bool) -> List[Peak]:
    """Parse mgf metabolomics data file.

    Parameters
    ----------
    input_string : str
        Metabolomics input data file
    charge : bool
        True if positive, False if negative

    Returns
    -------
    peaks : List[Peak]
        A list of Peak objects
    """
    peaks = []
    ms2 = []
    for line in input_string.split("\n"):
        sl = line.strip(" \r\n").split("=")
        if sl[0] == "PEPMASS":
            mass = sl[1]
        elif sl[0] == "TITLE":
            name = sl[1]
        elif sl[0] == "RTINSECONDS":
            r_time = sl[1]
        elif sl[0] == "END IONS":
            peaks.append(Peak(name, r_time, mass, charge, "False", ms2=ms2))
            ms2 = []
        else:
            try:
                mz, i = sl[0].split("\t")
                ms2.append((float(mz), float(i)))
            except ValueError:
                continue
    return peaks


def read_msp(input_string: str, charge: bool) -> List[Peak]:
    """Parse msp metabolomics data file.

    Parameters
    ----------
    input_string : str
        Metabolomics input data file
    charge : bool
        True if positive, False if negative

    Returns
    -------
    peaks : List[Peak]
        A list of Peak objects
    """
    peaks = []
    for spec in input_string.strip().split("\n\n"):
        ms2 = []
        inchikey = "False"
        r_time = 0
        name = "N/A"
        for line in spec.split("\n"):
            sl = line.split(": ")
            sl[0] = sl[0].replace(" ", "").replace("/", "").upper()
            if sl[0] == "PRECURSORMZ":
                mass = sl[1]
            elif sl[0] == "NAME":
                name = sl[1]
            # elif sl[0] == "RETENTIONTIME":
            # r_time = sl[1]
            # elif sl[0] == "IONMODE":
            # charge = sl[1].capitalize()
            elif sl[0] == "INCHIKEY":
                inchikey = sl[1]
            elif line and line[0].isdigit():
                try:
                    row = re.split("[\t ]", line)
                    ms2.append((float(row[0]), float(row[1])))
                except ValueError:
                    continue
        peaks.append(Peak(name, r_time, mass, charge, inchikey, ms2=ms2))
    return peaks


def read_mzxml(input_string: str, charge: bool) -> List[Peak]:
    """Parse mzXML metabolomics data file.

    Parameters
    ----------
    input_string : str
        Metabolomics input data file
    charge : bool
        True if positive, False if negative

    Returns
    -------
    List[Peak]
        A list of Peak objects
    """
    peaks = []
    root = ET.fromstring(input_string)
    prefix = root.tag.strip("mzXML")

    for scan in root.findall(f".//{prefix}scan"):
        # somewhat counter intuitively we will get the peak info from the
        # second fragments precursor info.
        if scan.attrib["msLevel"] == "2":
            precursor = scan.find(f"./{prefix}precursorMz")
            mz = precursor.text
            r_time = scan.attrib["retentionTime"][2:-1]
            name = f"{mz} @ {r_time}"
            charge = scan.attrib["polarity"]
            peaks.append(Peak(name, r_time, mz, charge, "False"))

    return peaks


class Struct:
    """convert key-value pairs into object-attribute pairs."""

    def __init__(self, **entries):
        self.__dict__.update(entries)


def ms_adduct_search(
    db: MINE,
    keggdb: pymongo.database.Database,
    text: str,
    text_type: str,
    ms_params) -> List:
    """Search for compound-adducts matching precursor mass.

    Parameters
    ----------
    db : MINE
        Contains compound documents to search.
    keggdb : pymongo.database.Database
        Contains models with associated compound documents.
    text : str
        Text as in metabolomics datafile for specific peak.
    text_type : str
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        text, assumes m/z values are separated by newlines (and set text_type
        to "form").
    ms_params : dict
        Specifies search settings, using the following key-value pairs:
        ------------------------
        Required Key-Value Pairs
        ------------------------
        "tolerance": float specifying tolerance for m/z, in mDa by default.
            Can specify in ppm if "ppm" key's value is set to True.
        "charge": bool (1 for positive, 0 for negative).
        ------------------------
        Optional Key-Value Pairs
        ------------------------
        "adducts": list of adducts to use. If not specified, uses all adducts.
        "models": List of model _ids. If supplied, score compounds higher if
            present in model. ["eco"] by default (E. coli).
        "ppm": bool specifying whether "tolerance" is in mDa or ppm. Default
            value for ppm is False (so tolerance is in mDa by default).
        "kovats": length 2 tuple specifying min and max kovats retention index
            to filter compounds (e.g. (500, 1000)).
        "logp": length 2 tuple specifying min and max logp to filter compounds
            (e.g. (-1, 2)).
        "halogens": bool specifying whether to filter out compounds containing
            F, Cl, or Br. Filtered out if set to True. False by default.

    Returns
    -------
    ms_adduct_output : list
        Compound JSON documents matching ms adduct query.
    """
    print(
        f"<MS Adduct Search: TextType={text_type}, Text={text}, Parameters={ms_params}>"
    )
    name = text_type + time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())

    if isinstance(ms_params, dict):
        ms_params = Struct(**ms_params)

    dataset = MetabolomicsDataset(name, ms_params)
    ms_adduct_output = []

    if text_type == "form":
        for mz in text.split("\n"):
            dataset.unk_peaks.append(Peak(mz, 0, float(mz), ms_params.charge, "False"))
    elif text_type == "mgf":
        dataset.unk_peaks = read_mgf(text, ms_params.charge)
    elif text_type == "mzXML" or text_type == "mzxml":
        dataset.unk_peaks = read_mzxml(text, ms_params.charge)
    elif text_type == "msp":
        dataset.unk_peaks = read_msp(text, ms_params.charge)
    else:
        raise IOError(f"{text_type} files not supported")

    if not ms_params.models:
        ms_params.models = ["eco"]

    dataset.native_set = get_KEGG_comps(db, keggdb, ms_params.models)
    dataset.annotate_peaks(db)

    for peak in dataset.unk_peaks:
        for hit in peak.isomers:
            if "CFM_spectra" in hit:
                del hit["CFM_spectra"]
            ms_adduct_output.append(hit)

    ms_adduct_output = score_compounds(
        db, ms_adduct_output, ms_params.models[0], parent_frac=0.75, reaction_frac=0.25
    )

    return ms_adduct_output


def ms2_search(
    db: MINE,
    keggdb: pymongo.database.Database,
    text: str,
    text_type: str,
    ms_params):
    """Search for compounds matching MS2 spectra.

    Parameters
    ----------
    db : MINE
        Contains compound documents to search.
    keggdb : pymongo.database.Database
        Contains models with associated compound documents.
    text : str
        Text as in metabolomics datafile for specific peak.
    text_type : str
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        text, assumes m/z values are separated by newlines (and set text_type
        to "form").
    ms_params : dict
        Specifies search settings, using the following key-value pairs:
        ------------------------
        Required Key-Value Pairs
        ------------------------
        "tolerance": float specifying tolerance for m/z, in mDa by default.
            Can specify in ppm if "ppm" key's value is set to True.
        "charge": bool (1 for positive, 0 for negative).
        "energy_level": int specifying fragmentation energy level to use. May
            be 10, 20, or 40.
        "scoring_function": str describing which scoring function to use. Can
            be either "jaccard" or "dot product".
        ------------------------
        Optional Key-Value Pairs
        ------------------------
        "adducts": list of adducts to use. If not specified, uses all adducts.
        "models": List of model _ids. If supplied, score compounds higher if
            present in model.
        "ppm": bool specifying whether "tolerance" is in mDa or ppm. Default
            value for ppm is False (so tolerance is in mDa by default).
        "kovats": length 2 tuple specifying min and max kovats retention index
            to filter compounds (e.g. (500, 1000)).
        "logp": length 2 tuple specifying min and max logp to filter compounds
            (e.g. (-1, 2)).
        "halogens": bool specifying whether to filter out compounds containing
            F, Cl, or Br. Filtered out if set to True. False by default.

    Returns
    -------
    ms_adduct_output : list
        Compound JSON documents matching ms2 search query.
    """
    print(f"<MS Adduct Sea" "rch: TextType={text_type}, Parameters={ms_params}>")
    name = text_type + time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())

    if isinstance(ms_params, dict):
        ms_params = Struct(**ms_params)

    dataset = MetabolomicsDataset(name, ms_params)
    ms_adduct_output = []

    if text_type == "form":
        split_form = [x.split() for x in text.strip().split("\n")]
        ms2_data = [(float(mz), float(i)) for mz, i in split_form[1:]]
        peak = Peak(
            split_form[0][0],
            0,
            float(split_form[0][0]),
            ms_params.charge,
            "False",
            ms2=ms2_data,
        )
        dataset.unk_peaks.append(peak)
    elif text_type == "mgf":
        dataset.unk_peaks = read_mgf(text, ms_params.charge)
    elif text_type == "mzXML":
        dataset.unk_peaks = read_mzxml(text, ms_params.charge)
    elif text_type == "msp":
        dataset.unk_peaks = read_msp(text, ms_params.charge)
    else:
        raise IOError(f"{text_type} files not supported")

    if not ms_params.models:
        ms_params.models = ["eco"]

    dataset.native_set = get_KEGG_comps(db, keggdb, ms_params.models)
    dataset.annotate_peaks(db)

    for peak in dataset.unk_peaks:

        if ms_params.scoring_function == "jaccard":
            if not ms_params.ppm:
                peak.score_isomers(
                    metric=jaccard,
                    energy_level=ms_params.energy_level,
                    tolerance=float(ms_params.tolerance) / 1000,
                )
            else:
                peak.score_isomers(metric=jaccard, energy_level=ms_params.energy_level)
        elif ms_params.scoring_function == "dot product":
            if not ms_params.ppm:
                peak.score_isomers(
                    metric=dot_product,
                    energy_level=ms_params.energy_level,
                    tolerance=float(ms_params.tolerance) / 1000,
                )
            else:
                peak.score_isomers(
                    metric=dot_product, energy_level=ms_params.energy_level
                )
        else:
            raise ValueError(
                'ms_params["scoring_function"] must be either '
                '"jaccard" or "dot product".'
            )

        for hit in peak.isomers:
            ms_adduct_output.append(hit)

        if ms_params.models:
            ms_adduct_output = score_compounds(
                db,
                ms_adduct_output,
                ms_params.models[0],
                parent_frac=0.75,
                reaction_frac=0.25,
            )

    return ms_adduct_output


def spectra_download(
    db: MINE,
    mongo_query: str = None,
    parent_filter: str = None,
    putative: bool = True) -> str:
    """Download one or more spectra for compounds matching a given query.

    Parameters
    ----------
    db : MINE
        Contains compound documents to search.
    mongo_query : str, optional (default: None)
        A valid Mongo query as a literal string. If None, all compound spectra
        are returned.
    parent_filter : str, optional (default: None)
        If set to a metabolic model's Mongo _id, only get spectra for compounds
        in or derived from that metabolic model.
    putative : bool, optional (default: True)
        If False, only find known compounds (i.e. in Generation 0). Otherwise,
        finds both known and predicted compounds.

    Returns
    -------
    spectral_library : str
        Text of all matching spectra, including headers and peak lists.
    """

    def print_peaklist(peaklist):
        text = [f"Num Peaks: {len(peaklist)}"]
        for x in peaklist:
            text.append(f"{x[0]} {x[1]}")
        text.append("")
        return text

    spectral_library = []
    msp_projection = {
        "MINE_id": 1,
        "Names": 1,
        "Mass": 1,
        "Generation": 1,
        "Inchikey": 1,
        "Formula": 1,
        "SMILES": 1,
        "Sources": 1,
        "Pos_CFM_spectra": 1,
        "Neg_CFM_spectra": 1,
    }

    if mongo_query:
        query_dict = literal_eval(mongo_query)
    else:
        query_dict = {}

    if not putative:
        query_dict["Generation"] = 0

    if parent_filter:
        model = db.models.find_one({"_id": parent_filter})

        if not model:
            raise ValueError("Invalid Model specified")

        parents = model["Compound_ids"]
        query_dict["$or"] = [
            {"_id": {"$in": parents}},
            {"Sources.Compound": {"$in": parents}},
        ]

    results = db.compounds.find(query_dict, msp_projection)

    for compound in results:

        # add header
        header = []
        if "Names" in compound and len(compound["Names"]) > 0:
            header.append(f"Name: {compound['Names'][0]}")
            for alt in compound["Names"][1:]:
                header.append(f"Synonym: {alt}")
        else:
            header.append("Name: MINE Compound {compound{'MINE_id']")

        for k, v in compound.items():
            if k not in {"Names", "Pos_CFM_spectra", "Neg_CFM_spectra"}:
                header.append(f"{k}: {v}")

        header.append("Instrument: CFM-ID")

        # add peak lists
        if "Pos_CFM_spectra" in compound:
            for energy, spec in compound["Pos_CFM_spectra"].items():
                # ??? not sure what James meant with this conditional:
                # he had spec_type as a argument to this function
                # if not spec_type or [True, int(energy[:2])] in spec_type:
                spectral_library += header
                spectral_library += [f"Ionization: Positive", "Energy: {energy}"]
                spectral_library += print_peaklist(spec)

        if "Neg_CFM_spectra" in compound:
            for energy, spec in compound["Neg_CFM_spectra"].items():
                # ??? not sure what James meant with this conditional:
                # if not spec_type or [False, int(energy[:2])] in spec_type:
                spectral_library += header
                spectral_library += [f"Ionization Mode: Negative", "Energy: {energy}"]
                spectral_library += print_peaklist(spec)

    spectral_library = "\n".join(spectral_library)

    return spectral_library
