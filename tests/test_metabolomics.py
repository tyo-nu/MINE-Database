"""Tests for metabolomics.py using pytest."""

import pytest
from pymongo.errors import ServerSelectionTimeoutError

from minedatabase.databases import MINE
from minedatabase.metabolomics import MetabolomicsDataset, Peak, dot_product, jaccard


# -------------------------------- Fixtures --------------------------------- #
@pytest.fixture()
def test_db():
    """Create a test MINE database. Created and torn down before and after each
    test it is used in."""
    try:
        testdb = MINE("mongotest")
    except ServerSelectionTimeoutError:
        print("No Mongo DB server detected")
    yield testdb


@pytest.fixture
def known_peaks():
    """Return a list of Peak instances with data (identified)."""
    peak1 = Peak(name='Test1Known', r_time=5.00, mz=867.1391, charge=True,
                 inchi_key='IRPOHFRNKHKIQA-UHFFFAOYSA-N')
    peak2 = Peak(name='Test2Known', r_time=8.00, mz=260.0297, charge=False,
                 inchi_key='HXXFSFRBOHSIMQ-FPRJBGLDSA-N')
    return [peak1, peak2]


@pytest.fixture
def unknown_peaks():
    """Return a list of Peak instances with data (unidentified)."""
    peak1 = Peak(name='Test1Unknown', r_time=3.00, mz=427.0294, charge=True,
                 ms2=[(10, 110), (20, 320), (25, 5)])
    peak2 = Peak(name='Test2Unknown', r_time=1.50, mz=180.0000, charge=False,
                 ms2=[(10, 105), (20, 50), (25, 90)])
    return [peak1, peak2]


@pytest.fixture
def native_cpd_set():
    """Return a set of native compound IDs."""
    return set(['C03e0b10e6490ce79a7b88cb0c4e17c2bf6204352',
                'C189c6f78587772ccbd0c2b9e118cf88ff35c316a',
                'C90ba1f8d4d84f6305539f6d05a74497d4d5dfe06'])


@pytest.fixture
def metabolomics_dataset(known_peaks, unknown_peaks, native_cpd_set):
    """Default instance of MetabolomicsDataset."""
    met_dataset = MetabolomicsDataset(name='Test Metabolomics Dataset',
                                      adducts=['[M-H]-', '[M+H]+'],
                                      known_peaks=known_peaks,
                                      unknown_peaks=unknown_peaks,
                                      native_set=native_cpd_set,
                                      ppm=False,
                                      tolerance=0.002,
                                      halogens=True,
                                      verbose=True)
    return met_dataset


# ------------------------ MetabolomicsDataset tests ------------------------ #
def test_metabolomics_dataset_adducts(metabolomics_dataset):
    """Test adducts of a metabolomics dataset.
    GIVEN an instance MetabolomicsDataset
    WHEN an it is initialized
    THEN make sure it has the correct adducts
    """
    assert metabolomics_dataset.pos_adducts == [('[M+H]+', 1.0, 1.007276)]
    assert metabolomics_dataset.neg_adducts == [('[M-H]-', 1.0, -1.007276)]


def test_metabolomics_dataset_str_representation(metabolomics_dataset):
    """Test MetabolomicsDataset object string representation.
    GIVEN an instance of MetabolomicsDataset
    WHEN the string representation is accessed (e.g. using str() or print())
    THEN make sure the metabolomics dataset name is returned
    """
    assert str(metabolomics_dataset) == metabolomics_dataset.name


def test_metabolomics_dataset_enumerate_possible_masses(metabolomics_dataset):
    """Make sure possible masses and ranges are calculated properly.
    GIVEN an instance of MetabolomicsDataset
    WHEN enumerate_possible_masses is used to get mass ranges
    THEN make sure those mass ranges are correctly calculated
    """
    metabolomics_dataset.enumerate_possible_masses(tolerance=0.001)
    # possible masses are each peak mass +/- each adduct
    possible_masses = set([178.992724, 181.007276, 426.022124, 428.036676])
    assert metabolomics_dataset.possible_masses == possible_masses
    # possible ranges are possible masses +/- tolerance
    possible_ranges = set([(426.02112400000004, 426.023124, 'Test1Unknown', 1.0),
                           (428.035676, 428.037676, 'Test1Unknown', 1.0),
                           (178.991724, 178.99372400000001, 'Test2Unknown', 1.0),
                           (181.00627599999999, 181.008276, 'Test2Unknown', 1.0)])
    assert set(metabolomics_dataset.possible_ranges) == possible_ranges


def test_metabolomics_dataset_get_rt(metabolomics_dataset):
    """Make sure we can get a Peak's retention time from its ID.
    GIVEN an instance of MetabolomicsDataset
    WHEN the retention time for a given peak in that dataset is requested
    THEN make sure the correct retention time is returned.
    """
    peak1_actual_rt = 5.00
    peak2_actual_rt = 8.00

    assert metabolomics_dataset.get_rt('Test1Known') == peak1_actual_rt
    assert metabolomics_dataset.get_rt('Test2Known') == peak2_actual_rt
    assert metabolomics_dataset.get_rt('InvalidID') is None



# -------------------------- Scoring function tests ------------------------- #
def test_dot_product():
    """Test dot product calculation.
    GIVEN two spectra
    WHEN their dot product is calculated
    THEN make sure it is the correct value
    """
    x = [(10, 100), (20, 300), (25, 10)]
    y = [(10.011, 10), (20.009, 300), (25, 100)]
    epsilon = 0.01

    calc_dot_product = dot_product(x, y, epsilon=epsilon)
    assert round(calc_dot_product, 5) == 0.90909


def test_jaccard():
    """Test jaccard index calculation.
    GIVEN two spectra
    WHEN their jaccard index is calculated
    THEN make sure it is the correct value
    """
    x = [(10, 100), (20, 300), (25, 10)]
    y = [(10.011, 10), (20.009, 300), (25, 100)]
    epsilon = 0.01

    calc_jaccard = jaccard(x, y, epsilon=epsilon)
    assert calc_jaccard == 0.5


# -------------------------------- Peak tests ------------------------------- #
def test_peak_str_representation(known_peaks):
    """Test Peak object string representation.
    GIVEN an instance of Peak
    WHEN the string representation is accessed (e.g. using str() or print())
    THEN make sure the peak name is returned
    """
    for peak in known_peaks:
        assert str(peak) == peak.name


def test_peak_score_isomers(unknown_peaks):
    """Test score_isomers method of Peak object.
    GIVEN an instance of Peak
    WHEN isomers are scored against mass spectra data
    THEN make sure they are in the correct order
    """
    # In reality pos and neg spectra would be different, but this is just for testing
    isomer1 = {'Pos_CFM_spectra': {'20 V': [(10, 100), (20, 300), (25, 10)]},
               'Neg_CFM_spectra': {'20 V': [(10, 100), (20, 300), (25, 10)]}}
    isomer2 = {'Pos_CFM_spectra': {'20 V': [(10, 10), (20, 55), (25, 100)]},
               'Neg_CFM_spectra': {'20 V': [(10, 10), (20, 55), (25, 100)]}}

    unknown_peak1, unknown_peak2 = unknown_peaks
    unknown_peak1.isomers = [isomer1, isomer2]
    unknown_peak2.isomers = [isomer1, isomer2]

    # Isomer order should get switched for peak2 but not for peak1
    unknown_peak1.score_isomers()
    unknown_peak2.score_isomers()

    assert unknown_peak1.isomers == [isomer1, isomer2]
    assert unknown_peak2.isomers == [isomer2, isomer1]
