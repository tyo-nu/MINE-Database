"""Tests for pickaxe.py using pytest."""
from pathlib import Path

import pytest

from minedatabase.filters import (
    AtomicCompositionFilter,
    MCSFilter,
    MetabolomicsFilter,
    MWFilter,
    SimilarityFilter,
    SimilaritySamplingFilter,
)

file_path = Path(__file__)
file_dir = file_path.parent
DATA_DIR = (file_dir / "../data/").resolve()

# check for eQ
try:
    from equilibrator_api import Q_

    from minedatabase.filters.thermodynamics import ThermoFilter
    from minedatabase.thermodynamics import Thermodynamics

    thermo = Thermodynamics()
    try:
        thermo.load_thermo_from_postgres()
        loaded_db = True
    except:
        try:
            thermo.load_thermo_from_sqlite()
            loaded_db = True
        except:
            pass
except:
    pass


def test_similarity_cutoff_single(pk_target):
    """Test similarity cutoff filter"""
    tani_threshold = 0.5
    _filter = SimilarityFilter(
        crit_similarity=tani_threshold, increasing_similarity=False
    )
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 355
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


@pytest.mark.skip("Heisenbug Test")
def test_filter_after(pk_target):
    """Test similarity cutoff filter"""
    tani_threshold = 0.5
    _filter = SimilarityFilter(
        crit_similarity=tani_threshold, increasing_similarity=False
    )
    pk_target.filter_after_final_gen = True
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 257
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == False
    )


def test_similarity_cutoff_multi(pk_target):
    """Test similarity cutoff filter"""
    tani_threshold = [0, 0.3, 0.5]
    _filter = SimilarityFilter(
        crit_similarity=tani_threshold, increasing_similarity=False
    )
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 1094
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


def test_similarity_cutoff_multi_short_list(pk_target):
    """Test similarity filter when the tani_threshold is shorter than generations."""
    tani_threshold = [0.5]
    _filter = SimilarityFilter(
        crit_similarity=tani_threshold, increasing_similarity=False
    )
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 355
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


def test_similarity_no_targets(pk_target):
    pk_target.target_smiles = []
    tani_threshold = 0.5
    _filter = SimilarityFilter(
        crit_similarity=tani_threshold, increasing_similarity=False
    )

    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 1348
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


def test_similarity_sample_default_weight(pk_target):
    """Test similarity cutoff filter"""
    _filter = SimilaritySamplingFilter(sample_size=10, weight=None)
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452


def test_similarity_sample_user_weight(pk_target):
    """Test similarity cutoff filter"""

    def weight(T):
        return T ** 4

    _filter = SimilaritySamplingFilter(sample_size=10, weight=weight)
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452


def test_similarity_sample_morgan(pk_target):
    """Test overwriting defaults"""

    fingerprint_method = "Morgan"
    fingerprint_args = {"radius": 2}

    _filter = SimilaritySamplingFilter(
        sample_size=10,
        fingerprint_method=fingerprint_method,
        fingerprint_args=fingerprint_args,
    )
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452 and len(pk_target.compounds) > 100


def test_similarity_sample_dice(pk_target):
    """Test overwriting defaults"""

    fingerprint_method = "Morgan"
    fingerprint_args = {"radius": 2}
    similarity_method = "Dice"

    _filter = SimilaritySamplingFilter(
        sample_size=10,
        fingerprint_method=fingerprint_method,
        fingerprint_args=fingerprint_args,
        similarity_method=similarity_method,
    )

    pk_target.filters.append(_filter)
    pk_target.transform_all(processes=2, generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452 and len(pk_target.compounds) > 100


def test_similarity_sample_multiprocess(pk_target):
    """Test similarity cutoff filter"""

    def weight(T):
        return T ** 4

    _filter = SimilaritySamplingFilter(sample_size=10, weight=weight)
    pk_target.react_targets = True
    pk_target.filters.append(_filter)
    pk_target.transform_all(processes=2, generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452


def test_MCS_list(pk_target):
    """Test similarity cutoff filter"""
    MCS_threshold = [0.1, 0.5]
    _filter = MCSFilter(crit_mcs=MCS_threshold)
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 340


@pytest.mark.skipif(not loaded_db, reason="No eQuilibrator DB found.")
def test_thermo_phys(pk_target):
    """Test thermo cutoff for physiological"""
    _filter = ThermoFilter(physiological=True)
    pk_target.filters.append(_filter)
    pk_target.transform_all(generations=1)

    assert True


def test_met_filter_mass(pk_target):
    """Test MetabolomicsFilter output without RT predictor."""
    metabolomics_data_path = DATA_DIR / "test_metabolomics/test_metabolomics_data.csv"
    met_filter = MetabolomicsFilter(
        filter_name="test_metabolomics_filter",
        met_data_name="test_metabolomics_data",
        met_data_path=metabolomics_data_path,
        possible_adducts=["[M+H]+", "[M-H]-"],
        mass_tolerance=0.001,
    )
    pk_target.filters.append(met_filter)
    pk_target.transform_all(generations=2)

    gen1_cpds = [
        pk_target.compounds[cpd]
        for cpd in pk_target.compounds
        if pk_target.compounds[cpd]["Generation"] == 1
    ]

    assert len(gen1_cpds) == 1
    assert gen1_cpds[0]["Matched_Peak_IDs"] == ["Test3"]
