"""Tests for rules.py using pytest."""


import json
from pathlib import Path

from minedatabase import pickaxe
from minedatabase.rules import BNICE, metacyc_generalized, metacyc_intermediate


file_path = Path(__file__)
file_dir = file_path.parent

with open((file_dir / "../data/test_rules/rules_to_assert.json"), "r") as f:
    rule_assert_dict = json.load(f)


def test_metacyc_generalized_full():
    rule_list, correactant_list, rule_name = metacyc_generalized()
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_generalized"
    assert len(pk.operators) == 1221
    assert len(pk.coreactants) == 45

    assert (
        pk.operators["rule0001"][1]["SMARTS"]
        == rule_assert_dict["Metacyc_generalized_rule0001"]
    )


def test_metacyc_generalized_specify_number():
    rule_list, correactant_list, rule_name = metacyc_generalized(n_rules=10)
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_generalized_10_rules"
    assert len(pk.operators) == 10
    assert len(pk.coreactants) == 45

    assert (
        pk.operators["rule0001"][1]["SMARTS"]
        == rule_assert_dict["Metacyc_generalized_rule0001"]
    )


def test_metacyc_generalized_specify_fraction():
    rule_list, correactant_list, rule_name = metacyc_generalized(fraction_coverage=0.9)
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_generalized_0,9_fraction_coverage"
    assert len(pk.operators) == 353
    assert len(pk.coreactants) == 45

    assert (
        pk.operators["rule0001"][1]["SMARTS"]
        == rule_assert_dict["Metacyc_generalized_rule0001"]
    )


def test_metacyc_exclude():
    rule_list, correactant_list, rule_name = metacyc_generalized(
        fraction_coverage=0.9, exclude_containing=["aromatic", "halogen"]
    )
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_generalized_0,9_fraction_coverage_with_exclusion"
    assert len(pk.operators) == 843
    assert len(pk.coreactants) == 45


def test_metacyc_include():
    rule_list, correactant_list, rule_name = metacyc_generalized(
        fraction_coverage=0.9, include_containing=["aromatic", "halogen"]
    )
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_generalized_0,9_fraction_coverage_with_inclusion"
    assert len(pk.operators) == 378
    assert len(pk.coreactants) == 45


def test_metacyc_intermediate():
    rule_list, correactant_list, rule_name = metacyc_intermediate()
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_intermediate"
    assert len(pk.operators) == 7358
    assert len(pk.coreactants) == 45

    assert (
        pk.operators["rule0001_0167"][1]["SMARTS"]
        == rule_assert_dict["Metacyc_intermediate_rule0001_0167"]
    )


def test_metacyc_intermediate_specify_number():
    rule_list, correactant_list, rule_name = metacyc_intermediate(n_rules=20)
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_intermediate_20_rules"
    assert len(pk.operators) == 20
    assert len(pk.coreactants) == 45


def test_metacyc_intermediate_specify_fraction():
    rule_list, correactant_list, rule_name = metacyc_intermediate(fraction_coverage=0.2)
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_intermediate_0,2_fraction_coverage"
    assert len(pk.operators) == 75
    assert len(pk.coreactants) == 45


def test_metacyc_intermediate_exclude():
    rule_list, correactant_list, rule_name = metacyc_intermediate(
        fraction_coverage=0.9, exclude_containing=["aromatic", "halogen"]
    )
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_intermediate_0,9_fraction_coverage_with_exclusion"
    assert len(pk.operators) == 5887
    assert len(pk.coreactants) == 45


def test_metacyc_intermediate_include():
    rule_list, correactant_list, rule_name = metacyc_intermediate(
        fraction_coverage=0.9, include_containing=["halogen"]
    )
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "Metacyc_intermediate_0,9_fraction_coverage_with_inclusion"
    assert len(pk.operators) == 84
    assert len(pk.coreactants) == 45


def test_BNICE():
    rule_list, correactant_list, rule_name = BNICE()
    pk = pickaxe.Pickaxe(rule_list=rule_list, coreactant_list=correactant_list)

    assert rule_name == "BNICE"
    assert len(pk.operators) == 250
    assert len(pk.coreactants) == 33

    assert pk.operators["1.1.1_01"][1]["SMARTS"] == rule_assert_dict["BNICE_1.1.1_01"]
