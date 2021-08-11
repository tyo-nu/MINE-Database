# Running feasibility test is an issue and crashes, do it here instead
from minedatabase.pickaxe import Pickaxe
from minedatabase.filters.feasibility import ReactionFeasibilityFilter
from pathlib import Path

cwd = Path(__file__).parent

data_path = cwd / "../../tests/data/test_filters"
pk = Pickaxe(
    rule_list=data_path / "test_filter_rules.tsv",
    coreactant_list=data_path / "metacyc_coreactants.tsv",
    filter_after_final_gen=True,
    quiet=True,
)
pk.load_compound_set(data_path / "test_filter_compounds.csv")
pk.load_targets(data_path / "test_filter_targets.csv")

_filter = ReactionFeasibilityFilter(use_unpredicted=False)
pk.filters.append(_filter)
pk.transform_all(generations=1)

len(pk.compounds) == 58
