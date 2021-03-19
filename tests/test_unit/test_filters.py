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


def test_tani_cutoff_single(pk_target):
    """Test tanimoto cutoff filter"""
    tani_threshold = 0.5
    filter = TanimotoFilter(crit_tani=tani_threshold, increasing_tani=False)
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 355
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


@pytest.mark.skip("Schrodingers Test")
def test_filter_after(pk_target):
    """Test tanimoto cutoff filter"""
    tani_threshold = 0.5
    filter = TanimotoFilter(crit_tani=tani_threshold, increasing_tani=False)
    pk_target.filter_after_final_gen = True
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 257
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == False
    )


def test_tani_cutoff_multi(pk_target):
    """Test tanimoto cutoff filter"""
    tani_threshold = [0, 0.3, 0.5]
    filter = TanimotoFilter(crit_tani=tani_threshold, increasing_tani=False)
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 1094
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


def test_tani_cutoff_multi_short_list(pk_target):
    """Test tanimoto filter when the tani_threshold is shorter than generations."""
    tani_threshold = [0.5]
    filter = TanimotoFilter(crit_tani=tani_threshold, increasing_tani=False)
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 355
    assert (
        pk_target.compounds["C779bfa0d747509f0499664b390657a336edec104"]["Expand"]
        == True
    )


def test_tani_sample_default_weight(pk_target):
    """Test tanimoto cutoff filter"""
    filter = TanimotoSamplingFilter(sample_size=10, weight=None)
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452


def test_tani_sample_user_weight(pk_target):
    """Test tanimoto cutoff filter"""

    def weight(T):
        return T ** 4

    filter = TanimotoSamplingFilter(sample_size=10, weight=weight)
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452


def test_tani_sample_multiprocess(pk_target):
    """Test tanimoto cutoff filter"""

    def weight(T):
        return T ** 4

    filter = TanimotoSamplingFilter(sample_size=10, weight=weight)
    pk_target.processes = 2
    pk_target.react_targets = True
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    # Filter must return less compounds than non-filter
    # Non-deterministic results, so no exact value can be used
    assert len(pk_target.compounds) < 1452


def test_MCS_list(pk_target):
    """Test tanimoto cutoff filter"""
    MCS_threshold = [0.1, 0.5]
    filter = MCSFilter(crit_mcs=MCS_threshold)
    pk_target.filters.append(filter)
    pk_target.transform_all(generations=2)

    assert len(pk_target.compounds) == 340
