import abc
import time
from copy import copy
from typing import List, Set

import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl

from minedatabase.pickaxe import Pickaxe


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")


class Filter(metaclass=abc.ABCMeta):
    """Abstract base class used to generate filters.

    The Filter class provides the framework for interaction with pickaxe expansions.
    Each filter subclass must inherit properties from the Filter class.
    All subclasses must implement properties and methods decorated with
    @abc.abstractmethod. Feel free to override other non-private methods as
    well, such as _pre_print() and _post_print().
    """

    @property
    @abc.abstractmethod
    def filter_name(self) -> str:
        """Obtain name of filter."""
        pass

    @abc.abstractmethod
    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
        """Return list of compounds to remove from pickaxe object.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        processes : int
            The number of processes to use, by default 1.
        generation : int
            Which generation the expansion is in.
        """
        pass

    def apply_filter(
        self,
        pickaxe: Pickaxe,
        processes: int = 1,
        generation: int = 0,
        print_on: bool = True,
    ) -> None:
        """Apply filter from Pickaxe object.

        Parameters
        ----------
        pickaxe : Pickaxe
            The Pickaxe object to filter.
        processes : int
            The number of processes to use, by default 1.
        print_on : bool
            Whether or not to print filtering results.
        """
        time_sample = time.time()

        self.generation = generation

        if print_on:
            n_total = self._get_n(pickaxe, "total")
            self._pre_print_header(pickaxe)
            self._pre_print()

        compound_ids_to_check, reaction_ids_to_check = self._choose_items_to_filter(
            pickaxe, processes
        )

        self._apply_filter_results(
            pickaxe, compound_ids_to_check, reaction_ids_to_check
        )

        if print_on:
            n_filtered = self._get_n(pickaxe, "filtered")
            self._post_print(pickaxe, n_total, n_filtered, time_sample)
            self._post_print_footer(pickaxe)

    def _pre_print_header(self, pickaxe: Pickaxe) -> None:
        """Print header before filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        """
        print("----------------------------------------")
        print(f"Filtering Generation {pickaxe.generation}\n")

    def _pre_print(self) -> None:
        """Print filter being applied."""
        print(f"Applying filter: {self.filter_name}")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
        """Print results of filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
            Unused here, but may be useful in your implementation.
        n_total : int
            Total number of compounds.
        n_filtered : int
            Number of compounds remaining after filtering.
        times_sample : float
            Time in seconds from time.time().
        """
        print(
            f"{n_filtered} of {n_total} compounds remain after applying "
            f"filter: {self.filter_name}"
            f"--took {round(time.time() - time_sample, 2)}s.\n"
        )

    def _post_print_footer(self, pickaxe: Pickaxe) -> None:
        """Print end of filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        """
        print(f"Done filtering Generation {pickaxe.generation}")
        print("----------------------------------------\n")

    def _get_n(self, pickaxe: Pickaxe, n_type: str) -> int:
        """Get current number of compounds to be filtered.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        n_type : str
            Whether to return "total" number of "filtered" number of compounds.

        Returns
        -------
        n : int
            Either the total or filtered number of compounds.
        """
        n = 0
        for cpd_dict in pickaxe.compounds.values():
            is_in_current_gen = cpd_dict["Generation"] == pickaxe.generation
            is_predicted_compound = cpd_dict["_id"].startswith("C")
            if is_in_current_gen and is_predicted_compound:
                if n_type == "total":
                    n += 1
                elif n_type == "filtered" and cpd_dict["Expand"]:
                    n += 1
        return n

    def _apply_filter_results(
        self,
        pickaxe: Pickaxe,
        compound_ids_to_check: List[str] = [],
        reaction_ids_to_delete: List[str] = [],
    ) -> None:
        """Apply filter results to Pickaxe object.

        Remove compounds and reactions that can be removed.
        For a compound to be removed it must:
            1. Not be flagged for expansion
            2. Not have a coproduct in a reaction marked for expansion
            3. Start with "C"

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network,
            this method modifies the Pickaxe object's compound documents.
        compound_ids_to_check : List[str]
            List of compound IDs to try to remove, if possible.
        """

        def should_delete_reaction(rxn_id: str) -> bool:
            """Returns whether or not a reaction can safely be deleted."""
            products = pickaxe.reactions[rxn_id]["Products"]
            for _, c_id in products:
                if c_id.startswith("C") and c_id not in cpds_to_remove:
                    return False
            # Every compound isn't in cpds_to_remove
            return True

        def remove_reaction(rxn_id):
            """Removes reaction and any resulting orphan compounds"""
            cpds_to_return = set()
            # Remove affiliations of reaction and check for orphans
            product_ids = [cpd[1] for cpd in pickaxe.reactions[rxn_id]["Products"]]
            for prod_id in product_ids:
                if prod_id.startswith("C"):
                    pickaxe.compounds[prod_id]["Product_of"].remove(rxn_id)
                    cpds_to_return.add(prod_id)
            compound_ids = [cpd[1] for cpd in pickaxe.reactions[rxn_id]["Reactants"]]
            for cpd_id in compound_ids:
                if cpd_id.startswith("C"):
                    pickaxe.compounds[cpd_id]["Reactant_in"].remove(rxn_id)
                    cpds_to_return.add(cpd_id)
            # Delete reaction itself
            del pickaxe.reactions[rxn_id]

            return cpds_to_return

        # Process reactions to delete
        # Loop through reactions to add compounds to check and to delete reactions
        if reaction_ids_to_delete:
            cpd_check_from_rxn = set()
            for rxn_id in reaction_ids_to_delete:
                cpd_check_from_rxn = cpd_check_from_rxn.union(remove_reaction(rxn_id))

            # Check for orphaned compounds due to reaction deletion
            while len(cpd_check_from_rxn) != 0:
                cpd_id = cpd_check_from_rxn.pop()
                # Orphan compound is one that has no reaction connecting it
                if cpd_id in pickaxe.compounds:
                    product_of = copy(pickaxe.compounds[cpd_id].get("Product_of", []))
                    # Delete if no reactions
                    if not product_of:
                        # Delete out reactions
                        reactant_in = copy(
                            pickaxe.compounds[cpd_id].get("Reactant_in", [])
                        )
                        for rxn_id in reactant_in:
                            cpd_check_from_rxn = cpd_check_from_rxn.union(
                                remove_reaction(rxn_id)
                            )
                        # Now delete compound
                        del pickaxe.compounds[cpd_id]

        # Go through compounds_ids_to_check and delete cpds/rxns as needed
        if compound_ids_to_check:
            cpds_to_remove = set()
            rxns_to_check = []

            compound_ids_to_check = set(compound_ids_to_check)
            for cpd_id in compound_ids_to_check:
                cpd_dict = pickaxe.compounds.get(cpd_id)
                if not cpd_dict:
                    continue

                if not cpd_dict["Expand"] and cpd_id.startswith("C"):
                    cpds_to_remove.add(cpd_id)

                    rxns_to_check.extend(pickaxe.compounds[cpd_id]["Product_of"])
                    rxns_to_check.extend(pickaxe.compounds[cpd_id]["Reactant_in"])

            rxns_to_check = set(rxns_to_check)
            # Function to check to see if should delete reaction
            # If reaction has compound that won't be deleted keep it
            # Check reactions for deletion
            for rxn_id in rxns_to_check:
                if should_delete_reaction(rxn_id):
                    for _, c_id in pickaxe.reactions[rxn_id]["Products"]:
                        if c_id.startswith("C"):
                            if rxn_id in pickaxe.compounds[c_id]["Product_of"]:
                                pickaxe.compounds[c_id]["Product_of"].remove(rxn_id)

                    for _, c_id in pickaxe.reactions[rxn_id]["Reactants"]:
                        if c_id.startswith("C"):
                            if rxn_id in pickaxe.compounds[c_id]["Reactant_in"]:
                                pickaxe.compounds[c_id]["Reactant_in"].remove(rxn_id)

                    del pickaxe.reactions[rxn_id]
                else:
                    # Reaction is dependent on compound that is flagged to be
                    # removed. Don't remove compound
                    products = pickaxe.reactions[rxn_id]["Products"]
                    cpds_to_remove -= set(i[1] for i in products)

                    # for _, c_id in products:
                    #     if c_id in cpds_to_remove:
                    #         cpds_to_remove -= {c_id}

            # Remove compounds and reactions if any found
            for cpd_id in cpds_to_remove:
                del pickaxe.compounds[cpd_id]


if __name__ == "__main__":
    pass
