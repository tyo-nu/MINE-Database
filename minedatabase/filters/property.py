from typing import Set, Union

import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from rdkit.Chem.AllChem import MolFromSmiles
from rdkit.Chem.Descriptors import ExactMolWt

from minedatabase.filters.base_filter import Filter
from minedatabase.pickaxe import Pickaxe


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")


class MWFilter(Filter):
    """A filter that removes compounds not within a MW range.

    This filter specifies a minimum and maximum molecular weight to create a range.
    Specifying None for either value will yield an unbounded MW range on that end.

    For example, specifying min_MW = None and max_MW = 1000 will give compounds less
    than or equal to 1000 g/mol.

    Parameters
    ----------
    min_MW : Union[float, None]
        Minimum MW in g/mol, by default None.
    max_MW : Union[float, None]
        Maximum MW in g/mol, by default None.

    Attributes
    ----------
    min_MW : Union[float, None]
        Minimum MW in g/mol.
    max_MW : Union[float, None]
        Maximum MW in g/mol.
    """

    def __init__(
        self,
        min_MW: Union[float, None] = None,
        max_MW: Union[float, None] = None,
    ) -> None:
        self._filter_name = "Molecular Weight"

        self.min_MW = min_MW or 0
        self.max_MW = max_MW or 100000

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Check the compounds against the MW constraints and return
        compounds to filter.
        """

        def MW_is_good(cpd):
            cpd_MW = ExactMolWt(MolFromSmiles(cpd["SMILES"]))
            return self.min_MW < cpd_MW and cpd_MW < self.max_MW

        def is_target(cpd, pickaxe):
            for t_id in pickaxe.targets:
                if "C" + t_id[1:] == cpd["_id"]:
                    return True
            return False

        cpds_remove_set = set()
        rxn_remove_set = set()

        print(
            f"Filtering Generation {pickaxe.generation} "
            f"with {self.min_MW} < MW < {self.max_MW}."
        )

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:
                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    if not MW_is_good(cpd):
                        cpds_remove_set.add(cpd["_id"])
                else:
                    if is_target(cpd, pickaxe):
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                    else:
                        if not MW_is_good(cpd):
                            cpds_remove_set.add(cpd["_id"])

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set, rxn_remove_set


class AtomicCompositionFilter(Filter):
    """Filter that removes compounds not within specific atomic composition ranges.

    This filter checks to see if the atomic composition of a compound is
    within a specified range. As an example, to only keep compounds with carbon between 4-7 and
    oxygen between 0-4 the following input would be used:

    atomic_composition = {"C": [4, 7], "O": [0, 4]}

    Parameters
    ----------
    atomic_composition : Dict
        A dictionary containing ranges for elemental composition. Of form
        {"Atom Symbol": [min, max]}, by default None.

    Attributes
    ----------
    atomic_composition : Dict
        A dictionary containing ranges for elemental composition.
    """

    def __init__(
        self,
        atomic_composition_constraints: dict = None,
    ) -> None:
        self._filter_name = "Atomic Composition"

        self.atomic_composition_constraints = atomic_composition_constraints

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Check the compounds against the atomic composition constraints and return
        compounds to filter.
        """

        def composition_is_good(cpd):
            atom_count = cpd["atom_count"]
            for atom in atom_count:
                atom_range = self.atomic_composition_constraints.get(atom)
                if atom_range:
                    atom_min = atom_range[0] or 0
                    atom_max = atom_range[1] or 10 ** 5

                    if not (
                        atom_min <= atom_count[atom] and atom_count[atom] <= atom_max
                    ):
                        return False

            return True

        def is_target(cpd, pickaxe):
            for t_id in pickaxe.targets:
                if "C" + t_id[1:] == cpd["_id"]:
                    return True
            return False

        cpds_remove_set = set()

        print(
            f"Filtering Generation {pickaxe.generation} "
            f"with atomic composition {self.atomic_composition_constraints}."
        )

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:
                # Check for targets and only react if terminal
                if pickaxe.react_targets:
                    if not composition_is_good(cpd):
                        cpds_remove_set.add(cpd["_id"])
                else:
                    if is_target(cpd, pickaxe):
                        pickaxe.compounds[cpd["_id"]]["Expand"] = False
                    else:
                        if not composition_is_good(cpd):
                            cpds_remove_set.add(cpd["_id"])

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set, []


if __name__ == "__main__":
    pass
