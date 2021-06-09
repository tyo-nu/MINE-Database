from typing import Set

import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from equilibrator_api import Q_

from minedatabase.filters.base_filter import Filter
from minedatabase.pickaxe import Pickaxe
from minedatabase.thermodynamics import Thermodynamics


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")


class ThermoFilter(Filter):
    """A filter that removes reactions and compounds with bad ∆Gr

    This filter allows for the specification of a pH, Ionic strength, pMg and using
    these to calculate ∆Gr. Reeactions are then filtered out based on ∆Gr.

    Parameters
    ----------
    eq_uri : str
        URI of ccache. Can be an sqlite or postgres URI. If no uri is given
        the default is the system default eQuilibrator uri.
    dg_max : float
        Maximum ∆Gr in kJ/mol, by default 0.
    pH : float
        pH of the expansion, by default 7
    ionic_strength : float
        ionic strength of the expansion, by default 0
    pMg : float
        pMg of the expansion, by default 3
    generation_list : list
        Generations to apply filter -- empty list filters all, by default empty list
    last_generation_only : bool
        Whether or not to only filter last generation, by default False

    Attributes
    ----------
    dg_max : Q_
    pH : Q_
        pH of the expansion
    ionic_strength : Q_
        ionic strength of the expansion
    pMg : Q_
        pMg of the expansion
    generation_list : list
        Generations to apply filter -- empty list filters all, by default empty list
    last_generation_only : bool
        Whether or not to only filter last generation, by default False
    """

    def __init__(
        self,
        eq_uri=None,
        dg_max=0,
        p_h=7,
        ionic_strength=0,
        p_mg=3,
        physiological=False,
        generation_list=[],
        last_generation_only=False,
    ) -> None:
        self._filter_name = "Thermodynamic Filter"
        self.dg_max = Q_(f"{dg_max}kJ/mol")
        self.p_h = Q_(f"{p_h}")
        self.ionic_strength = Q_(f"{ionic_strength}M")
        self.p_mg = Q_(f"{p_mg}")
        self.physiological = physiological
        self.generation_list = generation_list
        self.last_generation_only = last_generation_only

        self.thermo = Thermodynamics()
        if not eq_uri:
            eq_uri = ""
        if "post" in eq_uri:
            self.thermo.load_thermo_from_postgres(eq_uri)
        elif "sql" in eq_uri:
            self.thermo.load_thermo_from_sqlite(eq_uri)
        else:
            self.thermo.load_thermo_from_sqlite()

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _pre_print(self) -> None:
        """Print before filtering."""
        print(f"Filter out reactions with ∆Gr < {self.dg_max}")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
        """Print after filtering."""
        print(
            (
                f"{n_filtered} of {n_total} "
                "compounds selected after thermodynamic filtering in time_sample."
            )
        )

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Check the compounds against the MW constraints and return
        compounds to filter.
        """
        cpds_remove_set = set()
        rxns_remove_set = set()

        # TODO put these statements together
        # No reactions to filter for
        if len(pickaxe.reactions) == 0:
            print("No reactions to calculate ∆Gr for.")
            return cpds_remove_set, rxns_remove_set

        if self.last_generation_only and pickaxe.generation != self.generation:
            print("Not filtering for this generation using thermodynamics.")
            return cpds_remove_set, rxns_remove_set

        if self.generation_list and (self.generation - 1) not in self.generation_list:
            print("Not filtering for this generation using thermodynamics.")
            return cpds_remove_set, rxns_remove_set

        print(
            f"Filtering Generation {pickaxe.generation} "
            f"with ∆G <= {self.dg_max} at pH={self.p_h}, "
            f"I={self.ionic_strength}, pMg={self.p_mg}"
        )

        reactions_to_check = []
        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:
                reactions_to_check.extend(cpd["Product_of"])

        reactions_to_check = set(reactions_to_check)

        for rxn_id in reactions_to_check:
            if self.physiological:
                rxn_dg = self.thermo.physiological_dg_prime_from_rid(
                    r_id=rxn_id, pickaxe=pickaxe
                )
            else:
                rxn_dg = self.thermo.dg_prime_from_rid(
                    r_id=rxn_id,
                    pickaxe=pickaxe,
                    p_h=Q_(f"{self.p_h}"),
                    ionic_strength=Q_(f"{self.ionic_strength}"),
                    p_mg=Q_(f"{self.p_mg}"),
                )
            if rxn_dg >= self.dg_max:
                rxns_remove_set.add(rxn_id)

        return cpds_remove_set, rxns_remove_set
