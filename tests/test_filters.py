"""Tests for pickaxe.py using pytest."""
import pytest
from rdkit.RDLogger import logger

from minedatabase.filters import (
    MCSFilter,
    TanimotoFilter,
    TanimotoSamplingFilter,
)
from minedatabase.pickaxe import Pickaxe


# Default to no errors
lg = logger()
lg.setLevel(4)


@pytest.fixture(scope="function")
def pk():
    """Generate pickaxe to expand and filter."""
    new_pk = Pickaxe(
        rule_list="./tests/data/test_filters/test_filter_rules.tsv",
        coreactant_list="./tests/data/test_filters/metacyc_coreactants.tsv",
        filter_after_final_gen=False,
    )
    new_pk.load_compound_set("./tests/data/test_filters/test_filter_compounds.csv")
    new_pk.load_targets("./tests/data/test_filters/test_filter_targets.csv")

    return new_pk


def test_tani_cutoff_single(pk):
    """Test tanimoto cutoff filter"""
    tani_threshold = 0.5
    filter = TanimotoFilter(
        filter_name="Tani", crit_tani=tani_threshold, increasing_tani=False
    )
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    assert len(pk.compounds) == 355
    assert pk.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"] == True


@pytest.mark.skip("Schrodingers Test")
def test_filter_after(pk):
    """Test tanimoto cutoff filter"""
    tani_threshold = 0.5
    filter = TanimotoFilter(
        filter_name="Tani", crit_tani=tani_threshold, increasing_tani=False
    )
    pk.filter_after_final_gen = True
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    assert len(pk.compounds) == 257
    assert pk.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"] == False


def test_tani_cutoff_multi(pk):
    """Test tanimoto cutoff filter"""
    tani_threshold = [0, 0.3, 0.5]
    filter = TanimotoFilter(
        filter_name="Tani", crit_tani=tani_threshold, increasing_tani=False
    )
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    assert len(pk.compounds) == 1094
    assert pk.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"] == True


def test_tani_cutoff_multi_short_list(pk):
    """Test tanimoto filter when the tani_threshold is shorter than generations."""
    tani_threshold = [0.5]
    filter = TanimotoFilter(
        filter_name="Tani", crit_tani=tani_threshold, increasing_tani=False
    )
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    assert len(pk.compounds) == 355
    assert pk.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"] == True


def test_tani_sample_default_weight(pk):
    """Test tanimoto cutoff filter"""
    filter = TanimotoSamplingFilter(filter_name="Tani", sample_size=10, weight=None)
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk.compounds) < 1452


def test_tani_sample_user_weight(pk):
    """Test tanimoto cutoff filter"""

    def weight(T):
        return T ** 4

    filter = TanimotoSamplingFilter(filter_name="Tani", sample_size=10, weight=weight)
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk.compounds) < 1452


def test_tani_sample_multiprocess(pk):
    """Test tanimoto cutoff filter"""

    def weight(T):
        return T ** 4

    filter = TanimotoSamplingFilter(filter_name="Tani", sample_size=10, weight=weight)
    pk.processes = 2
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk.compounds) < 1452


def test_MCS_list(pk):
    """Test tanimoto cutoff filter"""
    MCS_threshold = [0.1, 0.5]
    filter = MCSFilter(
        filter_name="MCS",
        crit_mcs=MCS_threshold,
    )
    pk.filters.append(filter)
    pk.transform_all(generations=2)

    assert len(pk.compounds) == 340
