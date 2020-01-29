"""This script reads a set of unknowns in a CSV, MGF, mzXML, or excel file and
queries a mongo database for potential matches. The matches are written back
to a CSV, HTML doc, or printed to STDOUT.

USAGE - metabolomics.py [options] run_name file_of_Unknowns
"""

import math
import pickle
import re
import sys
import time
import xml.etree.ElementTree as ET
from optparse import OptionParser
from ast import literal_eval

import pandas as pd
import numpy as np

from app import app
from minedatabase.databases import establish_db_client
from minedatabase.utils import score_compounds


class MetabolomicsDataset:
    """A class containing all the information for a metabolomics data set"""
    # pylint: disable=redefined-outer-name
    def __init__(self, name, options):
        self.name = name
        self.options = options  # Object containing all options from optparser

        dtype = np.dtype("S20, f8, f8")
        all_pos_adducts = np.loadtxt(app.config['POS_ADDUCT_PATH'],
                                     dtype=dtype)
        all_neg_adducts = np.loadtxt(app.config['NEG_ADDUCT_PATH'],
                                     dtype=dtype)

        if hasattr(options, 'adducts'):
            pos_adducts = filter(lambda x: x[0] in options.adducts,
                                 all_pos_adducts)
            self.pos_adducts = np.array(pos_adducts, dtype=dtype)
            neg_adducts = filter(lambda x: x[0] in options.adducts,
                                 all_neg_adducts)
            self.neg_adducts = np.array(neg_adducts, dtype=dtype)
        else:
            if hasattr(options, 'positive_adduct_file'):
                self.pos_adducts = np.loadtxt(options.positive_adduct_file,
                                              dtype=dtype)
            else:
                self.pos_adducts = all_pos_adducts

            if hasattr(options, 'negative_adduct_file'):
                self.neg_adducts = np.loadtxt(options.negative_adduct_file,
                                              dtype=dtype)
            else:
                self.neg_adducts = all_neg_adducts

        if hasattr(options, 'kovats'):
            self.min_kovats = options.kovats[0]
            self.max_kovats = options.kovats[1]

        if hasattr(options, 'logP'):
            self.min_logp = options.logP[0]
            self.max_logp = options.logP[1]

        self.hit_projection = {'Formula': 1, 'MINE_id': 1, 'logP': 1,
                               'minKovatsRI': 1, 'maxKovatsRI': 1,
                               'NP_likeness': 1, 'Names': 1, 'SMILES': 1,
                               'Inchikey': 1, 'Generation': 1,
                               'Pos_CFM_spectra': 1, 'Neg_CFM_spectra': 1,
                               "Sources": 1}
        self.known_peaks = []  # contains Peak objects for knowns
        self.unk_peaks = []  # contains Peak objects for unknowns
        self.clusters = []  # tuples of formula and list of matching peaks
        self.known_set = set()  # contains InChI key of all known compounds
        self.native_set = set()  # contains _ids of compounds in model set
        self.total_formulas = 0
        self.total_hits = 0
        self.matched_peaks = 0

    def __str__(self):
        return self.name

    def find_db_hits(self, peak, db, adducts):
        """This function searches the database for matches of a peak given
        adducts and updates a peak and dataset file with that information."""
        # find nominal mass for a given m/z for each adducts and the max and
        # min values for db
        potential_masses = (peak.mz - adducts['f2']) / adducts['f1']
        if self.options.ppm:
            precision = (self.options.tolerance / 100000.) * potential_masses
        else:
            precision = self.options.tolerance * 0.001
        upper_bounds = potential_masses + precision
        lower_bounds = potential_masses - precision

        # search database for hits in the each adducts mass range that have no
        # innate charge.
        for i, adduct in enumerate(adducts):
            # build the query by adding the optional terms
            query_terms = [{"Mass": {"$gte": float(lower_bounds[i])}},
                           {"Mass": {"$lte": float(upper_bounds[i])}},
                           {'Charge': 0}]
            if hasattr(self, 'min_logp'):
                query_terms += [{"logP": {"$gte": self.min_logp}},
                                {"logP": {"$lte": self.max_logp}}]
            if hasattr(self, 'min_kovats'):
                query_terms += [{"maxKovatsRI": {"$gte": self.min_kovats}},
                                {"minKovatsRI": {"$lte": self.max_kovats}}]
            if adduct['f0'] == '[M]+':
                query_terms[2] = {'Charge': 1}
            for compound in db.compounds.find({"$and": query_terms},
                                              self.hit_projection):
                # Filters out halogens if the flag is enabled by moving to the
                # next compound before the current compound is counted or
                # stored.
                if not self.options.halogens:
                    if re.search('F[^e]|Cl|Br', compound['Formula']):
                        continue

                # update the total hits for the peak and make a note if the
                # compound is in the native_set
                peak.total_hits += 1
                if compound['_id'] in self.native_set:
                    peak.native_hit = True
                    compound['native_hit'] = True
                if compound['Generation'] < peak.min_steps:
                    peak.min_steps = compound['Generation']
                peak.formulas.add(compound['Formula'])
                compound['adduct'] = adduct['f0']
                compound['peak_name'] = peak.name
                peak.isomers.append(compound)

    def annotate_peaks(self, db):
        """This function iterates the through the unknown peaks in the dataset
        searches the database for compounds that match a peak m/z given the
        adducts permitted. Statistics on the annotated data set are printed"""
        for i, peak in enumerate(self.unk_peaks):

            positive = peak.charge == '+' or peak.charge == 'Positive' \
                or (peak.charge and isinstance(peak.charge, bool))
            negative = peak.charge == '-' or peak.charge == 'Negative' \
                or (not peak.charge and isinstance(peak.charge, bool))

            if positive:
                self.find_db_hits(peak, db, self.pos_adducts)
            elif negative:
                self.find_db_hits(peak, db, self.neg_adducts)
            else:
                raise ValueError("Invalid compound charge specification. "
                                 "Please use \"+\" or \"Positive\" for "
                                 "positive ions and \"-\" or \"Negative\"for "
                                 "negative ions.")

            if peak.total_hits > 0:
                self.matched_peaks += 1
                self.total_hits += peak.total_hits
                self.total_formulas += len(peak.formulas)
            if self.options.verbose:
                pct_done = int(float(i) / float(len(self.unk_peaks)) * 100)
                print("%s percent of peaks processed" % pct_done)

        if self.options.verbose:
            print("Proposed matches for %s of %s peaks"
                  % (self.matched_peaks, len(self.unk_peaks)))

            try:
                print("Average hits per peak: %s"
                      % (float(self.total_hits) / float(self.matched_peaks)))
                print("Average formulas per peak: %s"
                      % (float(self.total_formulas)
                         / float(self.matched_peaks)))
            except ZeroDivisionError:
                pass


def dot_product(x, y, epsilon=0.01):
    """Calculate the dot_product of two spectra"""
    z = 0
    n_v1 = 0
    n_v2 = 0

    for int1, int2 in approximate_matches(x, y, epsilon):
        z += int1 * int2
        n_v1 += int1 * int1
        n_v2 += int2 * int2

    return z / (math.sqrt(n_v1) * math.sqrt(n_v2))


def jaccard(x, y, epsilon=0.01):
    """Calculate the Jaccard Index of two spectra"""
    intersect = 0

    for val1, val2 in approximate_matches(x, y, epsilon):
        if val1 and val2:
            intersect += 1

    return intersect / float((len(x) + len(y) - intersect))


def approximate_matches(list1, list2, epsilon=0.01):
    """
    Takes two list of tuples and searches for matches of tuples first value
    within the supplied epsilon. Emits tuples with the tuples second values
    where found. if a value in one dist does not match the other list, it is
    emitted alone.
    :param list1: first list of tuples
    :type list1: list
    :param list2: second list of tuples
    :type list2: list
    :param epsilon: maximum difference in
    :type epsilon: float
    :return: second values of tuples
    :rtype: generator
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
        else:
            raise AssertionError('Unexpected else taken')


class Peak:
    """A class holding information about an unknown peak"""
    def __init__(self, name, r_time, mz, charge, inchi_key, ms2=[]):
        self.name = name
        self.r_time = float(r_time)  # retention time
        self.mz = float(mz)  # mass to charge ratio
        self.charge = charge  # polarity of charge
        self.inchi_key = inchi_key  # Inchikey, if known
        self.ms2peaks = ms2
        self.isomers = []
        self.formulas = set()
        self.total_hits = 0
        self.native_hit = False
        self.min_steps = 99

    def __str__(self):
        return self.name

    def score_isomers(self, metric=dot_product, energy_level=20,
                      tolerance=0.005):
        """
        Calculates the cosign similarity score between the provided ms2 peak
        list and pre-calculated CFM-spectra and sorts the isomer list
        according to this metric.
        :param metric: The scoring metric to use for the spectra. Function
         must accept 2 lists of (mz,intensity) tuples and return a score.
         Defaults to dot_product.
        :type energy_level: function
        :param energy_level: The Fragmentation energy level to use. May be 10,
         20 or 40. Defaults to 20
        :type energy_level: int
        :param tolerance: The precision to use for matching m/z in mDa
        :type energy_level: float
        :return:
        :rtype:
        """
        if not self.ms2peaks:
            raise ValueError('The ms2 peak list is empty')
        if self.charge:
            spec_key = "Pos_CFM_spectra"
        else:
            spec_key = "Neg_CFM_spectra"

        for i, hit in enumerate(self.isomers):
            if spec_key in hit:
                hit_spec = hit[spec_key]['%s V' % energy_level]
                score = metric(self.ms2peaks, hit_spec, epsilon=tolerance)
                rounded_score = round(score * 1000)
                self.isomers[i]['Spectral_score'] = rounded_score
                del hit[spec_key]
            else:
                self.isomers[i]['Spectral_score'] = None
        self.isomers.sort(key=lambda x: x['Spectral_score'], reverse=True)


def get_modelseed_comps(kb_db, model_ids):
    """Get MongoIDs from KBase database for compounds in model(s)."""
    comp_ids, _ids = set(), set()
    for model_id in model_ids:
        comp_ids = kb_db.models.find_one({'_id': model_id},
                                         {'Compounds': 1})['Compounds']
        for comp_id in comp_ids:
            comp_ids.add(comp_id)
    for comp_id in comp_ids:
        comp_id = int(comp_id.split('.')[1])
        for comp in kb_db.compounds.find({'Model_SEED': comp_id}, {'_id': 1}):
            _ids.add(comp['_id'])
    return _ids


def get_KEGG_comps(db, kegg_db, model_ids):
    """Get KEGG IDs from KEGG MINE database for compounds in model(s)."""
    kegg_ids, _ids = set(), set()
    for model_id in model_ids:
        comp_ids = kegg_db.models.find_one({'_id': model_id},
                                           {'Compounds': 1})['Compounds']
        for comp_id in comp_ids:
            kegg_ids.add(comp_id)
    for kegg_id in kegg_ids:
        for comp in db.compounds.find({'DB_links.KEGG': kegg_id}, {'_id': 1}):
            _ids.add(comp['_id'])
    return _ids


def sort_nplike(dic_list):
    """Sort in descending order of natural product likeness."""
    return sorted(dic_list, key=lambda x: float(x['NP_likeness']),
                  reverse=True)


def read_adduct_names(filepath):
    """Read adduct names from text file at specified path into a list.

    Parameters
    ----------
    filepath : str
        Path to adduct text file.

    Returns
    -------
    adducts : list
        Names of adducts in text file.
    """

    with open(filepath) as infile:
        adducts = [line.split(' \t')[0] for line in infile
                   if not line[0] == '#']

    return adducts


def read_mgf(input_string, charge):
    """Parse mgf metabolomics data file."""
    peaks = []
    ms2 = []
    for line in input_string.split('\n'):
        sl = line.strip(' \r\n').split('=')
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
                mz, i = sl[0].split('\t')
                ms2.append((float(mz), float(i)))
            except ValueError:
                continue
    return peaks


def read_msp(input_string, charge):
    """Parse msp metabolomics data file."""
    peaks = []
    for spec in input_string.strip().split('\n\n'):
        ms2 = []
        inchikey = "False"
        r_time = 0
        for line in spec.split('\n'):
            sl = line.split(': ')
            sl[0] = sl[0].replace(' ', '').replace('/', '').upper()
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
                    row = re.split('[\t ]', line)
                    ms2.append((float(row[0]), float(row[1])))
                except ValueError:
                    continue
        peaks.append(Peak(name, r_time, mass, charge, inchikey, ms2=ms2))
    return peaks


def read_mzxml(input_string, charge):
    """Parse mzXML metabolomics data file."""
    peaks = []
    tree = ET.fromstring(input_string)
    root = tree.getroot()
    prefix = root.tag.strip('mzXML')

    for scan in root.findall('.//%sscan' % prefix):
        # somewhat counter intuitively we will get the peak info from the
        # second fragments precursor info.
        if scan.attrib['msLevel'] == '2':
            precursor = scan.find('./%sprecursorMz' % prefix)
            mz = precursor.text
            r_time = scan.attrib['retentionTime'][2:-1]
            name = "%s @ %s" % (mz, r_time)
            charge = scan.attrib['polarity']
            peaks.append(Peak(name, r_time, mz, charge, "False"))

    return peaks


class Struct:
    """convert key-value pairs into object-attribute pairs."""
    def __init__(self, **entries):
        self.__dict__.update(entries)


def ms_adduct_search(db, keggdb, text, text_type, ms_params):
    """Search for compound-adducts matching precursor mass.

    Parameters
    ----------
    db : Mongo DB
        Contains compound documents to search.
    keggdb : Mongo DB
        Contains models with associated compound documents.
    text : str
        Text as in metabolomics datafile for specific peak.
    text_type : str, optional (default: None)
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        None, assumes m/z values are separated by newlines.
    ms_params : dict
        Specifies search settings, using the following key-value pairs:
        ------------------------
        Required Key-Value Pairs
        ------------------------
        'tolerance': float specifying tolerance for m/z, in mDa by default.
            Can specify in ppm if 'ppm' key's value is set to True.
        'charge': bool (1 for positive, 0 for negative).
        ------------------------
        Optional Key-Value Pairs
        ------------------------
        'adducts': list of adducts to use. If not specified, uses all adducts.
        'models': List of model _ids. If supplied, score compounds higher if
            present in model.
        'ppm': bool specifying whether 'tolerance' is in mDa or ppm. Default
            value for ppm is False (so tolerance is in mDa by default).
        'kovats': length 2 tuple specifying min and max kovats retention index
            to filter compounds (e.g. (500, 1000)).
        'logp': length 2 tuple specifying min and max logp to filter compounds
            (e.g. (-1, 2)).
        'halogens': bool specifying whether to filter out compounds containing
            F, Cl, or Br. Filtered out if set to True. False by default.

    Returns
    -------
    ms_adduct_output : list
        Compound JSON documents matching ms adduct query.
    """
    print("<MS Adduct Search: TextType=%s, Text=%s, Parameters=%s>"
          % (text_type, text, ms_params))
    name = text_type + time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())

    if isinstance(ms_params, dict):
        ms_params = Struct(**ms_params)

    dataset = MetabolomicsDataset(name, ms_params)
    ms_adduct_output = []

    if text_type == 'form':
        for mz in text.split('\n'):
            dataset.unk_peaks.append(Peak(mz, 0, float(mz), ms_params.charge,
                                          "False"))
    elif text_type == 'mgf':
        dataset.unk_peaks = read_mgf(text, ms_params.charge)
    elif text_type == 'mzXML':
        dataset.unk_peaks = read_mzxml(text, ms_params.charge)
    elif text_type == 'msp':
        dataset.unk_peaks = read_msp(text, ms_params.charge)
    else:
        raise IOError('%s files not supported' % text_type)

    dataset.native_set = get_KEGG_comps(db, keggdb, ms_params.models)
    dataset.annotate_peaks(db)

    for peak in dataset.unk_peaks:
        for hit in peak.isomers:
            if 'CFM_spectra' in hit:
                del hit['CFM_spectra']
            ms_adduct_output.append(hit)

    if ms_params.models:
        ms_adduct_output = score_compounds(db, ms_adduct_output,
                                           ms_params.models[0],
                                           parent_frac=.75, reaction_frac=.25)

    return ms_adduct_output


def ms2_search(db, keggdb, text, text_type, ms_params):
    """Search for compounds matching MS2 spectra.

    Parameters
    ----------
    db : Mongo DB
        Contains compound documents to search.
    keggdb : Mongo DB
        Contains models with associated compound documents.
    text : str
        Text as in metabolomics datafile for specific peak.
    text_type : str, optional (default: None)
        Type of metabolomics datafile (mgf, mzXML, and msp are supported). If
        None, assumes m/z values are separated by newlines.
    ms_params : dict
        Specifies search settings, using the following key-value pairs:
        ------------------------
        Required Key-Value Pairs
        ------------------------
        'tolerance': float specifying tolerance for m/z, in mDa by default.
            Can specify in ppm if 'ppm' key's value is set to True.
        'charge': bool (1 for positive, 0 for negative).
        'energy_level': int specifying fragmentation energy level to use. May
            be 10, 20, or 40.
        'scoring_function': str describing which scoring function to use. Can
            be either 'jaccard' or 'dot product'.
        ------------------------
        Optional Key-Value Pairs
        ------------------------
        'adducts': list of adducts to use. If not specified, uses all adducts.
        'models': List of model _ids. If supplied, score compounds higher if
            present in model.
        'ppm': bool specifying whether 'tolerance' is in mDa or ppm. Default
            value for ppm is False (so tolerance is in mDa by default).
        'kovats': length 2 tuple specifying min and max kovats retention index
            to filter compounds (e.g. (500, 1000)).
        'logp': length 2 tuple specifying min and max logp to filter compounds
            (e.g. (-1, 2)).
        'halogens': bool specifying whether to filter out compounds containing
            F, Cl, or Br. Filtered out if set to True. False by default.

    Returns
    -------
    ms_adduct_output : list
        Compound JSON documents matching ms2 search query.
    """
    print("<MS Adduct Sea""rch: TextType=%s, Parameters=%s>"
          % (text_type, ms_params))
    name = text_type + time.strftime("_%d-%m-%Y_%H:%M:%S", time.localtime())

    if isinstance(ms_params, dict):
        ms_params = Struct(**ms_params)

    dataset = MetabolomicsDataset(name, ms_params)
    ms_adduct_output = []

    if text_type == 'form':
        split_form = [x.split() for x in text.strip().split('\n')]
        ms2_data = [(float(mz), float(i)) for mz, i in split_form[1:]]
        peak = Peak(split_form[0][0], 0, float(split_form[0][0]),
                    ms_params.charge, "False", ms2=ms2_data)
        dataset.unk_peaks.append(peak)
    elif text_type == 'mgf':
        dataset.unk_peaks = read_mgf(text, ms_params.charge)
    elif text_type == 'mzXML':
        dataset.unk_peaks = read_mzxml(text, ms_params.charge)
    elif text_type == 'msp':
        dataset.unk_peaks = read_msp(text, ms_params.charge)
    else:
        raise IOError('%s files not supported' % text_type)

    dataset.native_set = get_KEGG_comps(db, keggdb, ms_params.models)
    dataset.annotate_peaks(db)

    for peak in dataset.unk_peaks:

        if ms_params.scoring_function == 'jaccard':
            if not ms_params.ppm:
                peak.score_isomers(metric=jaccard,
                                   energy_level=ms_params.energy_level,
                                   tolerance=float(ms_params.tolerance) / 1000)
            else:
                peak.score_isomers(metric=jaccard,
                                   energy_level=ms_params.energy_level)
        elif ms_params.scoring_function == 'dot product':
            if not ms_params.ppm:
                peak.score_isomers(metric=dot_product,
                                   energy_level=ms_params.energy_level,
                                   tolerance=float(ms_params.tolerance) / 1000)
            else:
                peak.score_isomers(metric=dot_product,
                                   energy_level=ms_params.energy_level)
        else:
            raise ValueError("ms_params['scoring_function'] must be either "
                             "'jaccard' or 'dot product'.")

        for hit in peak.isomers:
            ms_adduct_output.append(hit)

        if ms_params.models:
            ms_adduct_output = score_compounds(db, ms_adduct_output,
                                               ms_params.models[0],
                                               parent_frac=.75,
                                               reaction_frac=.25)

    return ms_adduct_output


def spectra_download(db, mongo_query=None, parent_filter=False, putative=True):
    """Download one or more spectra for compounds matching a given query.

    Parameters
    ----------
    db : Mongo DB
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
        text = ["Num Peaks: %s" % len(peaklist)]
        for x in peaklist:
            text.append("%s %s" % (x[0], x[1]))
        text.append("")
        return text

    spectral_library = []
    msp_projection = {'MINE_id': 1, 'Names': 1, 'Mass': 1, 'Generation': 1,
                      'Inchikey': 1, 'Formula': 1, 'SMILES': 1, 'Sources': 1,
                      'Pos_CFM_spectra': 1, 'Neg_CFM_spectra': 1}

    if mongo_query:
        query_dict = literal_eval(mongo_query)
    else:
        query_dict = {}

    if not putative:
        query_dict['Generation'] = 0

    if parent_filter:
        model = db.models.find_one({"_id": parent_filter})

        if not model:
            raise ValueError('Invalid Model specified')

        parents = model["Compound_ids"]
        query_dict['$or'] = [{'_id': {'$in': parents}},
                             {'Sources.Compound': {'$in': parents}}]

    results = db.compounds.find(query_dict, msp_projection)

    for compound in results:

        # add header
        header = []
        if "Names" in compound and len(compound['Names']) > 0:
            header.append("Name: %s" % compound['Names'][0])
            for alt in compound['Names'][1:]:
                header.append("Synonym: %s" % alt)
        else:
            header.append("Name: MINE Compound %s" % compound['MINE_id'])

        for k, v in compound.items():
            if k not in {"Names", "Pos_CFM_spectra", "Neg_CFM_spectra"}:
                header.append("%s: %s" % (k, v))

        header.append("Instrument: CFM-ID")

        # add peak lists
        if 'Pos_CFM_spectra' in compound:
            for energy, spec in compound['Pos_CFM_spectra'].items():
                # ??? not sure what James meant with this conditional:
                # he had spec_type as a argument to this function
                # if not spec_type or [True, int(energy[:2])] in spec_type:
                spectral_library += header
                spectral_library += ["Ionization: Positive",
                                     "Energy: %s" % energy]
                spectral_library += print_peaklist(spec)

        if 'Neg_CFM_spectra' in compound:
            for energy, spec in compound['Neg_CFM_spectra'].items():
                # ??? not sure what James meant with this conditional:
                # if not spec_type or [False, int(energy[:2])] in spec_type:
                spectral_library += header
                spectral_library += ["Ionization Mode: Negative",
                                     "Energy: %s" % energy]
                spectral_library += print_peaklist(spec)

    spectral_library = "\n".join(spectral_library)

    return spectral_library


if __name__ == '__main__':
    # pylint: disable=invalid-name

    sys.exit()
    t_start = time.time()
    # This block handles user flags and arguments. For more information see the
    # optparse API documentation

    usage = "usage: %prog [options] run_name"
    parser = OptionParser(usage)
    parser.add_option("-c", "--charge", dest="unknowns_charge",
                      default='Positive', help='Charge for unknowns. Options '
                      'are "Positive" or "Negative" (case-sensitive).')
    parser.add_option("-p", "--positive_adducts", dest="positive_adduct_file",
                      default="Batch Adduct Query/Positive Adducts.txt",
                      help="The path to the desired positive adducts file")
    parser.add_option("-n", "--negative_adducts", dest="negative_adduct_file",
                      default="Batch Adduct Query/Negative Adducts.txt",
                      help="The path to the desired negative adducts file")
    parser.add_option("-d", "--database", dest="database", default="1GenKEGG",
                      help="The name of the database to search")
    parser.add_option("-m", "--modelSEED", dest="modelSEED", default="null",
                      help="The model SEED id of the organism for the data "
                      "set")
    parser.add_option("-k", "--known", dest="known_file", default="null",
                      help="The path to a file containing known or targeted "
                      "peaks")
    parser.add_option("-t", "--tolerance", dest="tolerance", type="float",
                      default=2, help="The m/z tolerance(precision) for peak "
                      "searching")
    parser.add_option("--ppm", dest="ppm", action="store_true", default=False,
                      help="Tolerance is in Parts Per Million")
    parser.add_option("-x", "--halogens", dest="halogens", action="store_true",
                      default=False, help="include compounds containing F, "
                      "Cl, and Br")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False, help="include compounds containing F, "
                      "Cl, and Br")

    (options, args) = parser.parse_args()
    # Ensure that the right number of arguments are passed to the script.
    if len(args) < 2:
        sys.exit("Missing required arguments Usage: [options] run_name, "
                 "path_to_unk_peaks")
    if len(args) > 2:
        sys.exit("Too many arguments submitted. Usage: [options] run_name, "
                 "path_to_unk_peaks")

    data = MetabolomicsDataset(args[0], options)
    unknowns_file = args[1]

    # if the user specifies known peaks import them for use in clustering and
    # ordering of isomers
    if options.known_file != 'null':
        print("Loading known peaks")
        data.known_peaks = list(pd.read_csv(options.known_file))
        for peak in data.known_peaks:
            data.known_set.add(peak.known)

    # detect unknown compounds file type and parse accordingly
    if ('.txt' or '.csv') in unknowns_file:
        data.unk_peaks = list(pd.read_csv(open(unknowns_file).read()))
    elif '.mgf' in unknowns_file:
        data.unk_peaks = read_mgf(open(unknowns_file).read(), options.charge)
    elif '.mzXML' in unknowns_file:
        data.unk_peaks = read_mzxml(open(unknowns_file).read(), options.charge)
    else:
        sys.exit("Unknown file type not recognised. Please use .mgf, .xlsx, "
                 ".txt or .csv file")

    client = establish_db_client()
    db = client[options.database]
    if not db.compounds.count():
        sys.exit('No compounds in supplied database')

    if options.modelSEED != 'null':
        print("Loading native compounds")
        kbase_db = client['KBase']
        data.native_set = get_modelseed_comps(kbase_db, [options.modelSEED])

    data.annotate_peaks(db)
    with open(data.name + '.pkl') as outfile:
        pickle.dump(data, outfile)

    t_end = time.time()
    print("BatchAdductQuery.py completed in %s seconds" % (t_end - t_start))
