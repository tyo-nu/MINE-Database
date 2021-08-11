"""Tests for metabolomics.py using pytest."""

from pathlib import Path

import pymongo
import pytest
from pymongo.errors import ServerSelectionTimeoutError

from minedatabase.metabolomics import (
    MetabolomicsDataset,
    Peak,
    dot_product,
    jaccard,
    read_mgf,
    read_msp,
    read_mzxml,
    spectra_download,
)


try:
    client = pymongo.MongoClient(ServerSelectionTimeoutMS=20)
    client.server_info()
    del client
    is_mongo = True
except ServerSelectionTimeoutError as err:
    is_mongo = False
valid_db = pytest.mark.skipif(not is_mongo, reason="No MongoDB Connection")

file_path = Path(__file__)
file_dir = file_path.parent

DATA_DIR = (file_dir / "../data/").resolve()

# -------------------------------- Fixtures --------------------------------- #
@pytest.fixture
def known_peaks():
    """Return a list of Peak instances with data (identified)."""
    peak1 = Peak(
        name="Test1Known",
        r_time=5.00,
        mz=867.1391,
        charge="+",
        inchi_key="IRPOHFRNKHKIQA-UHFFFAOYSA-N",
    )
    peak2 = Peak(
        name="Test2Known",
        r_time=8.00,
        mz=260.0297,
        charge="-",
        inchi_key="HXXFSFRBOHSIMQ-FPRJBGLDSA-N",
    )
    return [peak1, peak2]


@pytest.fixture
def unknown_peaks():
    """Return a list of Peak instances with data (unidentified)."""
    peak1 = Peak(
        name="Test1Unknown",
        r_time=3.00,
        mz=427.0294,
        charge="+",
        ms2=[(10, 110), (20, 320), (25, 5)],
    )
    peak2 = Peak(
        name="Test2Unknown",
        r_time=1.50,
        mz=180.0000,
        charge="-",
        ms2=[(10, 105), (20, 50), (25, 90)],
    )
    return [peak1, peak2]


@pytest.fixture
def native_cpd_set():
    """Return a set of native compound IDs."""
    return set(
        [
            "C03e0b10e6490ce79a7b88cb0c4e17c2bf6204352",
            "C189c6f78587772ccbd0c2b9e118cf88ff35c316a",
            "C90ba1f8d4d84f6305539f6d05a74497d4d5dfe06",
        ]
    )


@pytest.fixture
def metabolomics_dataset(known_peaks, unknown_peaks, native_cpd_set):
    """Default instance of MetabolomicsDataset."""
    met_dataset = MetabolomicsDataset(
        name="Test Metabolomics Dataset",
        adducts=["[M-H]-", "[M+H]+"],
        known_peaks=known_peaks,
        unknown_peaks=unknown_peaks,
        native_set=native_cpd_set,
        ppm=False,
        tolerance=0.002,
        halogens=True,
        verbose=True,
    )
    return met_dataset


# ------------------------ MetabolomicsDataset tests ------------------------ #
def test_metabolomics_dataset_adducts(metabolomics_dataset):
    """Test adducts of a metabolomics dataset.
    GIVEN an instance MetabolomicsDataset
    WHEN an it is initialized
    THEN make sure it has the correct adducts
    """
    assert metabolomics_dataset.pos_adducts == [("[M+H]+", 1.0, 1.007276)]
    assert metabolomics_dataset.neg_adducts == [("[M-H]-", 1.0, -1.007276)]


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
    possible_masses = {"+": set([426.022124]), "-": set([181.007276])}
    assert metabolomics_dataset.possible_masses == possible_masses
    # possible ranges are possible masses +/- tolerance
    possible_ranges = {
        "+": [(426.02112400000004, 426.023124, "Test1Unknown", "[M+H]+")],
        "-": [(181.00627599999999, 181.008276, "Test2Unknown", "[M-H]-")],
    }
    assert metabolomics_dataset.possible_ranges == possible_ranges


def test_metabolomics_dataset_get_rt(metabolomics_dataset):
    """Make sure we can get a Peak's retention time from its ID.
    GIVEN an instance of MetabolomicsDataset
    WHEN the retention time for a given peak in that dataset is requested
    THEN make sure the correct retention time is returned
    """
    peak1_actual_rt = 5.00
    peak2_actual_rt = 8.00

    assert metabolomics_dataset.get_rt("Test1Known") == peak1_actual_rt
    assert metabolomics_dataset.get_rt("Test2Known") == peak2_actual_rt
    assert metabolomics_dataset.get_rt("InvalidID") is None


# TODO Jon fix me!
@valid_db
def test_metabolomics_dataset_find_db_hits(test_db, metabolomics_dataset):
    """Search for expected metaoblomics hits in test database
    GIVEN a MINE database and metabolomics dataset
    WHEN that database is searched against a specific peak in that dataset
    THEN make sure the correct hits are returned
    """
    peak = Peak(name="test", r_time=1, mz=261.0369946, charge="+")
    adducts = [("[M+H]+", 1, 1.007276), ("[M-H]-", 1, -1.007276)]
    metabolomics_dataset.find_db_hits(peak, test_db, adducts)
    assert len(peak.isomers) == 1


@valid_db
def test_metabolomics_dataset_annotate_peaks(test_db, metabolomics_dataset):
    """Uses find_db_hits to try to annotate all unknown peaks in dataset
    GIVEN a metabolomics dataset and MINE db
    WHEN trying to find hits for every unknown peak in that dataset
    THEN make sure all peaks get the correct # of hits
    """
    if not test_db.compounds.find_one({"_id": "mass_test"}):
        test_db.compounds.insert_one(
            {
                "_id": "mass_test",
                "Mass": 181.007276,
                "Charge": 0,
                "Formula": "C6H12O2N",
                "Generation": 1,
            }
        )
    metabolomics_dataset.annotate_peaks(test_db)
    for peak in metabolomics_dataset.unknown_peaks:
        if peak.name == "Test1Unknown":
            assert len(peak.isomers) == 0
        elif peak.name == "Test2Unknown":
            assert len(peak.isomers) == 1


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
    isomer1 = {
        "Pos_CFM_spectra": {"20 V": [(10, 100), (20, 300), (25, 10)]},
        "Neg_CFM_spectra": {"20 V": [(10, 100), (20, 300), (25, 10)]},
    }
    isomer2 = {
        "Pos_CFM_spectra": {"20 V": [(10, 10), (20, 55), (25, 100)]},
        "Neg_CFM_spectra": {"20 V": [(10, 10), (20, 55), (25, 100)]},
    }

    unknown_peak1, unknown_peak2 = unknown_peaks
    unknown_peak1.isomers = [isomer1, isomer2]
    unknown_peak2.isomers = [isomer1, isomer2]

    # Isomer order should get switched for peak2 but not for peak1
    unknown_peak1.score_isomers()
    unknown_peak2.score_isomers()

    assert unknown_peak1.isomers == [isomer1, isomer2]
    assert unknown_peak2.isomers == [isomer2, isomer1]


# --------------------- Metabolomics File Parsing Tests --------------------- #


def test_read_mgf():
    """Test reading peak data from an MGF metabolomics file.
    GIVEN an MGF file
    WHEN that MGF file is parsed into a list of Peak objects
    THEN make sure those Peak objects are correct
    """
    mgf_path = DATA_DIR / "test_metabolomics/test.mgf"
    with open(mgf_path, "r") as infile:
        mgf_data = infile.read()

    peaks = read_mgf(mgf_data, charge="+")
    test_peak1 = peaks[0]
    test_peak2 = peaks[1]

    assert test_peak1.charge == "+"
    assert test_peak1.mz == 100.5
    assert test_peak1.ms2peaks == [(40.4, 2.3), (50.3, 12.3), (102.1, 10.0)]
    assert test_peak1.r_time == 420

    assert test_peak2.charge == "+"
    assert test_peak2.mz == 131.0
    assert test_peak2.ms2peaks == [(93.8, 12.4), (250.1, 1000.2)]
    assert test_peak2.r_time == 305


def test_read_msp():
    """Test reading peak data from an MSP metabolomics file.
    GIVEN an MSP file
    WHEN that MSP file is parsed into a list of Peak objects
    THEN make sure those Peak objects are correct
    """
    msp_path = DATA_DIR / "test_metabolomics/test.msp"
    with open(msp_path, "r") as infile:
        msp_data = infile.read()

    peaks = read_msp(msp_data, charge="+")
    test_peak1 = peaks[0]
    test_peak2 = peaks[1]

    assert test_peak1.charge == "+"
    assert test_peak1.mz == 100.5
    assert test_peak1.ms2peaks == [(40.4, 2.3), (50.3, 12.3), (102.1, 10.0)]
    assert test_peak1.r_time == 420

    assert test_peak2.charge == "+"
    assert test_peak2.mz == 131.0
    assert test_peak2.ms2peaks == [(93.8, 12.4), (250.1, 1000.2)]
    assert test_peak2.r_time == 305


def test_read_mzxml():
    """Test reading peak data from an mzXML metabolomics file.
    GIVEN an mzXML file
    WHEN that mzXML file is parsed into a list of Peak objects
    THEN make sure those Peak objects are correct
    """
    mzxml_path = DATA_DIR / "test_metabolomics/test.mzXML"
    with open(mzxml_path, "r") as infile:
        mzxml_data = infile.read()

    peaks = read_mzxml(mzxml_data, charge="+")
    test_peak1 = peaks[0]
    test_peak2 = peaks[1]

    assert test_peak1.charge == "+"
    assert test_peak1.mz == 100.5
    assert test_peak1.r_time == 420

    assert test_peak2.charge == "+"
    assert test_peak2.mz == 131.0
    assert test_peak2.r_time == 305


# ----------------------------- Spectra Download ---------------------------- #


@valid_db
def test_spectra_download(test_db):
    """Test download of spectra from MINE database.
    GIVEN a MINE database
    WHEN spectra are downloaded in string format
    THEN make sure the downloaded spectra is correct
    """
    if not test_db.compounds.find_one({"_id": "ms2_test"}):
        test_spectra = {
            "10V": [(1, 2), (3, 4), (5, 6)],
            "20V": [(6, 7)],
            "40V": [(100, 1000), (250, 200)],
        }
        test_db.compounds.insert_one(
            {
                "_id": "ms2_test",
                "Mass": 260,
                "Charge": 0,
                "Formula": "C12H12O2N",
                "Generation": 1,
                "Pos_CFM_spectra": test_spectra,
            }
        )
    cpd_query = '{"_id": "ms2_test"}'
    spectra = spectra_download(test_db, cpd_query)

    spectra_path = DATA_DIR / "test_metabolomics/test_spectra.txt"

    with open(spectra_path, "r") as infile:
        spectra_data = infile.read()

    assert spectra == spectra_data
