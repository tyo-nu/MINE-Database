"""Basic thermodynamic calculations for pickaxe."""

from typing import Union

import pint
from equilibrator_api import (
    Q_,
    ComponentContribution,
    Reaction,
    default_physiological_ionic_strength,
    default_physiological_p_h,
    default_physiological_p_mg,
    default_physiological_temperature,
)
from equilibrator_api.phased_reaction import PhasedReaction
from equilibrator_assets.compounds import Compound
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache.compound_cache import CompoundCache
from pymongo import MongoClient
from sqlalchemy import create_engine

from minedatabase.pickaxe import Pickaxe


class Thermodynamics:
    """Class to calculate thermodynamics of Pickaxe runs.

    Thermodynamics allows for the calculation of:
        1) Standard ∆G' of formation
        2) Standard ∆G'o of reaction
        3) Physiological ∆G'm of reaction
        4) Adjusted ∆G' of reaction

    eQuilibrator objects can also be obtained from r_ids and c_ids.

    Parameters
    ----------
    mongo_uri: str
        URI of the mongo database.
    client: MongoClient
        Connection to Mongo.
    CC: ComponentContribution
        eQuilibrator Component Contribution object to calculate ∆G with.
    lc: LocalCompoundCache
        The local compound cache to generate eQuilibrator compounds from.
    """

    def __init__(
        self,
    ):
        # Mongo params
        self.mongo_uri = None
        self.client = None
        self._core = None

        # eQ params
        self.CC = ComponentContribution()
        self.lc = None
        self._water = None

    def load_mongo(self, mongo_uri: Union[str, None] = None):
        if mongo_uri:
            self.mongo_uri = mongo_uri
            self.client = MongoClient(mongo_uri)
        else:
            self.mongo_uri = "localhost:27017"
            self.client = MongoClient()

        self._core = self.client["core"]

    def _all_dbs_loaded(self):
        if self.client and self._core and self.lc:
            return True
        else:
            print("Load connection to Mongo and eQuilibrator local cache.")
            return False

    def _eq_loaded(self):
        if self.lc:
            return True
        else:
            print("Load eQulibrator local cache.")
            return False

    def _reset_CC(self):
        """reset CC back to defaults"""
        self.CC.p_h = default_physiological_p_h
        self.CC.p_mg = default_physiological_p_mg
        self.CC.temperature = default_physiological_temperature
        self.CC.ionic_strength = default_physiological_ionic_strength

    def load_thermo_from_postgres(
        self, postgres_uri: str = "postgresql:///eq_compounds"
    ) -> None:
        """Load a LocalCompoundCache from a postgres uri for equilibrator.

        Parameters
        ----------
        postgres_uri : str, optional
            uri of the postgres DB to use, by default "postgresql:///eq_compounds"
        """
        self.lc = LocalCompoundCache()
        self.lc.ccache = CompoundCache(create_engine(postgres_uri))

        self._water = self.lc.get_compounds("O")

    def load_thermo_from_sqlite(
        self, sqlite_filename: str = "compounds.sqlite"
    ) -> None:
        """Load a LocalCompoundCache from a sqlite file for equilibrator.

        compounds.sqlite can be generated through LocalCompoundCache's method
        generate_local_cache_from_default_zenodo

        Parameters
        ----------
        sqlite_filename: str
            filename of the sqlite file to load.
        """
        self.lc = LocalCompoundCache()
        self.lc.load_cache(sqlite_filename)

        self._water = self.lc.get_compounds("O")

    def get_eQ_compound_from_cid(
        self, c_id: str, pickaxe: Pickaxe = None, db_name: str = None
    ) -> Union[Compound, None]:
        """Get an equilibrator compound for a given c_id from the core.

        Attempts to retrieve a compound from the core or a specified db_name.

        Parameters
        ----------
        c_id : str
            compound ID for MongoDB lookup of a compound.
        pickaxe : Pickaxe
            pickaxe object to look for the compound in, by default None.
        db_name : str
            Database to look for compound in before core database, by default None.

        Returns
        -------
        equilibrator_assets.compounds.Compound
            eQuilibrator Compound
        """
        # Find locally in pickaxe
        compound_smiles = None
        if pickaxe:
            if c_id in pickaxe.compounds:
                compound_smiles = pickaxe.compounds[c_id]["SMILES"]
            else:
                return None

        # Find in mongo db
        elif self._all_dbs_loaded():
            if db_name:
                compound = self.client[db_name].compounds.find_one(
                    {"_id": c_id}, {"SMILES": 1}
                )
                if compound:
                    compound_smiles = compound["SMILES"]

            # No cpd smiles from database name
            if not compound_smiles:
                compound = self._core.compounds.find_one({"_id": c_id}, {"SMILES": 1})
                if compound:
                    compound_smiles = compound["SMILES"]

        # No compound_smiles at all
        if not compound_smiles or "*" in compound_smiles:
            return None
        else:
            eQ_compound = self.lc.get_compounds(
                compound_smiles, bypass_chemaxon=True, save_empty_compounds=True
            )
            return eQ_compound

    def standard_dg_formation_from_cid(
        self, c_id: str, pickaxe: Pickaxe = None, db_name: str = None
    ) -> Union[float, None]:
        """Get standard ∆Gfo for a compound.

        Parameters
        ----------
        c_id : str
            Compound ID to get the ∆Gf for.
        pickaxe : Pickaxe
            pickaxe object to look for the compound in, by default None.
        db_name : str
            Database to look for compound in before core database, by default None.

        Returns
        -------
        Union[float, None]
            ∆Gf'o for a compound, or None if unavailable.
        """
        eQ_cpd = self.get_eQ_compound_from_cid(c_id, pickaxe, db_name)
        if not eQ_cpd:
            return None
        dgf = self.CC.standard_dg_formation(eQ_cpd)
        dgf = dgf[0]

        return dgf

    def get_eQ_reaction_from_rid(
        self, r_id: str, pickaxe: Pickaxe = None, db_name: str = None
    ) -> Union[PhasedReaction, None]:
        """Get an eQuilibrator reaction object from an r_id.

        Parameters
        ----------
        r_id : str
            Reaction id to get object for.
        pickaxe : Pickaxe
            pickaxe object to look for the compound in, by default None.
        db_name : str
            Database to look for reaction in.

        Returns
        -------
        PhasedReaction
            eQuilibrator reactiono to calculate ∆Gr with.
        """
        if pickaxe:
            if r_id in pickaxe.reactions:
                reaction_info = pickaxe.reactions[r_id]
            else:
                return None
        elif db_name:
            mine = self.client[db_name]
            reaction_info = mine.reactions.find_one({"_id": r_id})
            if not reaction_info:
                return None
        else:
            return None

        reactants = reaction_info["Reactants"]
        products = reaction_info["Products"]

        lhs = " + ".join(f"{r[0]} {r[1]}" for r in reactants)
        rhs = " + ".join(f"{p[0]} {p[1]}" for p in products)
        reaction_string = " => ".join([lhs, rhs])

        compounds = set([r[1] for r in reactants])
        compounds.update(tuple(p[1] for p in products))

        eQ_compound_dict = {
            c_id: self.get_eQ_compound_from_cid(c_id, pickaxe, db_name)
            for c_id in compounds
        }

        if not all(eQ_compound_dict.values()):
            return None

        if "X73bc8ef21db580aefe4dbc0af17d4013961d9d17" not in compounds:
            eQ_compound_dict["water"] = self._water

        eq_reaction = Reaction.parse_formula(eQ_compound_dict.get, reaction_string)

        return eq_reaction

    def physiological_dg_prime_from_rid(
        self, r_id: str, pickaxe: Pickaxe = None, db_name: str = None
    ) -> Union[pint.Measurement, None]:
        """Calculate the ∆Gm' of a reaction.

        Parameters
        ----------
        r_id : str
            ID of the reaction to calculate.
        pickaxe : Pickaxe
            pickaxe object to look for the compound in, by default None.
        db_name : str
            MINE the reaction is found in.

        Returns
        -------
        pint.Measurement
            The calculated ∆G'm.
        """
        eQ_reaction = self.get_eQ_reaction_from_rid(r_id, pickaxe, db_name)
        if not eQ_reaction:
            return None
        dGm_prime = self.CC.physiological_dg_prime(eQ_reaction)

        return dGm_prime

    def standard_dg_prime_from_rid(
        self, r_id: str, pickaxe: Pickaxe = None, db_name: str = None
    ) -> Union[pint.Measurement, None]:
        """Calculate the ∆G'o of a reaction.

        Parameters
        ----------
        r_id : str
            ID of the reaction to calculate.
        pickaxe : Pickaxe
            pickaxe object to look for the compound in, by default None.
        db_name : str
            MINE the reaction is found in.

        Returns
        -------
        pint.Measurement
            The calculated ∆G'o.
        """
        eQ_reaction = self.get_eQ_reaction_from_rid(r_id, pickaxe, db_name)
        if not eQ_reaction:
            return None
        dG0_prime = self.CC.standard_dg_prime(eQ_reaction)

        return dG0_prime

    def dg_prime_from_rid(
        self,
        r_id: str,
        pickaxe: Pickaxe = None,
        db_name: str = None,
        p_h: Q_ = default_physiological_p_h,
        p_mg: Q_ = default_physiological_p_mg,
        ionic_strength: Q_ = default_physiological_ionic_strength,
    ) -> Union[pint.Measurement, None]:
        """Calculate the ∆G' of a reaction.

        Parameters
        ----------
        r_id : str
            ID of the reaction to calculate.
        pickaxe : Pickaxe
            pickaxe object to look for the compound in, by default None.
        db_name : str
            MINE the reaction is found in.
        p_h : Q_
            pH of system.
        p_mg: Q_
            pMg of the system.
        ionic_strength: Q_
            ionic strength of the system.

        Returns
        -------
        pint.Measurement
            The calculated ∆G'.
        """
        eQ_reaction = self.get_eQ_reaction_from_rid(r_id, pickaxe, db_name)
        if not eQ_reaction:
            return None

        self.CC.p_h = p_h
        self.CC.p_mg = p_mg
        self.CC.ionic_strength = ionic_strength

        dG_prime = self.CC.dg_prime(eQ_reaction)

        self._reset_CC()

        return dG_prime
