"""Databases.py: This file contains MINE database classes including database
loading and writing functions."""
import datetime
import multiprocessing
import os
import platform
import sys
from copy import deepcopy
from math import ceil
from shutil import move
from subprocess import call
from typing import List, Union

import pymongo
from pymongo.errors import ServerSelectionTimeoutError
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger

from minedatabase import utils


# from minedatabase.NP_Score import npscorer as nps

# nps_model = nps.readNPModel()

lg = logger()
lg.setLevel(4)


def establish_db_client(uri: str = None) -> pymongo.MongoClient:
    """Establish a connection to a mongo database given a URI.

    Uses the provided URI to connect to a mongoDB. If none is given
    the default URI is used when using pymongo.

    Parameters
    ----------
    uri : str, optional
        URI to connect to mongo DB, by default None.

    Returns
    -------
    pymongo.MongoClient
        Connection to the specified mongo instance.

    Raises
    ------
    IOError
        Attempt to connect to database timed out.
    """
    # If a connection string is given use that, otherwise go to local
    try:
        if uri:
            client = pymongo.MongoClient(uri, ServerSelectionTimeoutMS=5)
        else:
            client = pymongo.MongoClient(ServerSelectionTimeoutMS=5)
    except ServerSelectionTimeoutError:
        raise IOError(
            "Failed to load database client. Please verify that " "mongod is running"
        )
    return client


class MINE:
    """
    This class provides an interface to the MongoDB and some useful functions.

    Parameters
    ----------
    name : str
        Name of the database to work with.
    uri : str, optional
        uri of the mongo server, by default "mongodb://localhost:27017/".

    Attributes
    ----------
    client : pymongo.MongoClient
        client connection to the MongoDB.
    compounds : Collection
        Compounds collection.
    core_compounds : Collection
        Core compounds collection.
    meta_data : Collection
        Metadata collection.
    models : Collection
        Models collection.
    name : str
        Name of the database
    operators : Collection
        Operators collection.
    reactions : Collection
        Reactions collection.
    target_compounds : Collection
        Target compounds collection.
    uri : str
        MongoDB connection string.
    """

    def __init__(self, name: str, uri: str = "mongodb://localhost:27017/") -> None:
        self.client = establish_db_client(uri)
        self.uri = uri
        self._db = self.client[name]
        self._core_db = self.client.core
        self.core_compounds = self._core_db.compounds
        self.name = name
        self.meta_data = self._db.meta_data
        self.compounds = self._db.compounds
        self.target_compounds = self._db.target_compounds
        self.reactions = self._db.reactions
        self.operators = self._db.operators
        self.models = self._db.models
        self.reactant_in = self._db.reactant_in
        self.product_of = self._db.product_of
        # self.nps_model = nps.readNPModel()
        self._mass_cache = {}  # for rapid calculation of reaction mass change

    def add_reaction_mass_change(self, reaction: str = None) -> Union[float, None]:
        """Calculate the change in mass between reactant and product compounds.

        This is useful for discovering compounds in molecular networking. If no reaction
        is specified then mass change of each reaction in the database will be
        calculated.

        Parameters
        ----------
        reaction : str, optional
            Reaction ID to calculate the mass change for, by default None.

        Returns
        -------
        float, optional
            Mass change of specified reaction. None if masses not all found.
        """

        def _get_mass(_id):
            if _id not in self._mass_cache:
                try:
                    self._mass_cache["_id"] = self.compounds.find_one(
                        {"_id": _id}, {"Mass": 1}
                    )["Mass"]
                except (TypeError, KeyError):
                    raise ValueError(
                        f"A mass value for {_id} was not found in "
                        "compounds collection"
                    )
            return self._mass_cache["_id"]

        def _get_mass_change(rxn):
            if isinstance(rxn["Reactants"][0], dict):
                key = "c_id"
            else:
                key = 1
            start_mass = _get_mass(rxn["Reactants"][0][key])
            return [
                _get_mass(product[key]) - start_mass
                for product in rxn["Products"]
                if product[key][0] == "C"
            ]

        if reaction:
            return _get_mass_change(reaction)

        bulk = self.reactions.initialize_unordered_bulk_op()
        reactions = self.reactions.find({"Mass_Change": {"$exists": 0}})
        for rxn in reactions:
            try:
                bulk.find({"_id": rxn["_id"]}).update(
                    {"$set": {"Mass_Change": _get_mass_change(rxn)}}
                )
            except ValueError as e:
                print(e.args[0])
                print(f"Failed to parse rxn {rxn['_id']}")
        bulk.execute()

    def generate_image_files(
        self,
        path: str,
        query: dict = None,
        dir_depth: int = 0,
        img_type: str = "svg:-a,nosource,w500,h500",
        convert_r: bool = False,
    ) -> None:
        """Generates image files for compounds in database using ChemAxon's MolConvert.

        Parameters
        ----------
        path : str
            Target directory for image file.
        query : dict, optional
            Query to limit number of files generated, by default None.
        dir_depth : int, optional
            The number of directory levels to split the compounds
            into for files system efficiency. Ranges from 0 (all in top
            level directory) to the length of the file name (40 for MINE hashes),
            by default 0.
        img_type : str, optional
            Type of image file to be generated. See molconvert
            documentation for valid options, by default 'svg:-a,nosource,w500,h500'.
        convert_r : bool, optional
            Convert R in the smiles to *, by default False.
        """
        ids = []
        extension = img_type.split(":")[0]
        structure_file = os.path.join(path, "tmp.smiles")
        if not query:
            query = {}

        if not os.path.exists(path):
            if sys.platform == "linux":
                os.system(f"sudo mkdir {path} -m 777")
            else:
                os.mkdir(path)
        with open(structure_file, "w") as outfile:
            for comp in self.compounds.find(query, {"SMILES": 1}):
                if convert_r:
                    outfile.write(f"{comp['SMILES'].replace('R', '*')}\n")
                else:
                    outfile.write(f"{comp['SMILES']}\n")
                ids.append(comp["_id"])
        if platform.system() == "Windows":
            rc = call(
                [f"molconvert -mo '{path}/.{extension}' {img_type} '{structure_file}'"],
                shell=True,
            )
        else:
            rc = call(
                [f"molconvert -mo '{path}/.{extension}' {img_type} '{structure_file}'"],
                shell=True,
            )
        if rc:
            raise RuntimeError(f"molconvert returned {rc}")
        os.remove(structure_file)

        for i, _id in enumerate(ids):
            old = os.path.join(path, f"{i+1}.{extension}")
            new = path
            for j in range(0, dir_depth):
                new = os.path.join(new, _id[j])
            if not os.path.exists(new):
                os.makedirs(new)
            new = os.path.join(new, _id + "." + extension)
            if os.path.isfile(old):
                move(old, new)

    def build_indexes(self) -> None:
        """Build indexes for efficient querying of the database."""
        self.core_compounds.drop_indexes()
        self.core_compounds.create_index([("Mass", pymongo.ASCENDING)])
        self.core_compounds.create_index("Inchikey")
        self.core_compounds.create_index("MINES")
        self.compounds.create_index("Inchikey")
        self.compounds.create_index("Inchi")
        self.compounds.create_index("SMILES")
        self.reactions.create_index("Reactants.c_id")
        self.reactions.create_index("Products.c_id")
        self.meta_data.insert_one(
            {"Timestamp": datetime.datetime.now(), "Action": "Database indexes built"}
        )


# Functions to write data to MINE
# Reactions
def write_reactions_to_mine(
    reactions: List[dict], db: MINE, chunk_size: int = 10000
) -> None:
    """Write reactions to reaction collection of MINE.

    Parameters
    ----------
    reactions : List[dict]
        Dictionary of reactions to write.
    db : MINE
        MINE object to write reactions with.
    chunk_size : int, optional
        Size of chunks to break reactions into when writing, by default 10000.
    """
    n_rxns = len(reactions)
    for i, rxn_chunk in enumerate(utils.Chunks(reactions, chunk_size)):
        if i % 20 == 0:
            print(f"Writing Reactions: Chunk {i} of {int(n_rxns/chunk_size) + 1}")
        rxn_requests = [
            pymongo.InsertOne(utils.convert_sets_to_lists(rxn_dict))
            for rxn_dict in rxn_chunk
        ]

        db.reactions.bulk_write(rxn_requests, ordered=False)


# Compounds
def write_compounds_to_mine(
    compounds: List[dict], db: MINE, chunk_size: int = 10000, processes: int = 1
) -> None:
    """Write compounds to reaction collection of MINE.

    Parameters
    ----------
    compounds : List[dict]
        Dictionary of compounds to write.
    db : MINE
        MINE object to write compounds with.
    chunk_size : int, optional
        Size of chunks to break compounds into when writing, by default 10000.
    processes : int, optional
        Number of processors to use, by default 1.
    """
    n_cpds = len(compounds)
    if processes == 1:
        pool = None
    else:
        pool = multiprocessing.Pool(processes)

    for i, cpd_chunk in enumerate(utils.Chunks(compounds, chunk_size)):
        if i % 20 == 0:
            print(f"Writing Compounds: Chunk {i} of {int(n_cpds/chunk_size) + 1}")

        cpd_requests = []
        reactant_in_requests = []
        product_of_requests = []

        if pool:
            for res in pool.imap_unordered(_get_cpd_insert, cpd_chunk):
                cpd_request, reactant_in_request, product_of_request = res
                cpd_requests.append(cpd_request)
                reactant_in_requests.extend(reactant_in_request)
                product_of_requests.extend(product_of_request)
        else:
            for res in map(_get_cpd_insert, cpd_chunk):
                cpd_request, reactant_in_request, product_of_request = res
                cpd_requests.append(cpd_request)
                reactant_in_requests.extend(reactant_in_request)
                product_of_requests.extend(product_of_request)

        db.compounds.bulk_write(cpd_requests, ordered=False)
        if reactant_in_requests:
            db.reactant_in.bulk_write(reactant_in_requests, ordered=False)
        if product_of_requests:
            db.product_of.bulk_write(product_of_requests, ordered=False)

    if pool:
        pool.close()
        pool.join()


def _get_cpd_insert(cpd_dict: dict):
    output_keys = [
        "_id",
        "ID",
        "SMILES",
        "InChI_key",
        "Type",
        "Generation",
        "Expand",
        "Reactant_in",
        "Product_of",
        "Matched_Peak_IDs",
        "Matched_Adducts",
        "Predicted_RT",
    ]

    # create Reactant_in
    reactant_in_requests = []
    product_of_requests = []
    insert_dict = {
        key: cpd_dict.get(key) for key in output_keys if cpd_dict.get(key) != None
    }
    if "Reactant_in" in insert_dict:
        chunked_reactant_in = _get_reactant_in_insert(cpd_dict)
        insert_dict["Reactant_in"] = []
        for r_in_dict in chunked_reactant_in:
            reactant_in_requests.append(pymongo.InsertOne(r_in_dict))
            insert_dict["Reactant_in"].append(r_in_dict["_id"])

    # create Product_of
    if "Product_of" in insert_dict:
        chunked_product_of = _get_product_of_insert(cpd_dict)
        insert_dict["Product_of"] = []
        for p_of_dict in chunked_product_of:
            product_of_requests.append(pymongo.InsertOne(p_of_dict))
            insert_dict["Product_of"].append(p_of_dict["_id"])

    cpd_request = pymongo.InsertOne(insert_dict)
    return cpd_request, reactant_in_requests, product_of_requests


def _get_reactant_in_insert(compound: dict) -> List[dict]:
    """Write reactants_in, ensuring memory size isn't too big.

    MongoDB only allows < 16 MB entries. This function breaks large reactants_in
    up to ensure this doesn't happen.

    Parameters
    ----------
    compounds : List[dict]
        Dictionary of compounds to write.
    db : MINE
        MINE object to write compounds with.

    Returns
    -------
    List[dict]
        dicts of reactant_in to insert
    """

    # Get number of chunks reactant_in must be broken up into
    # 16 MB is the max for BSON, cut to 14 MB max just to be safe
    # Also some weirdness is that max_size must be 1 order of magnitude lower
    max_size = 1.4 * 10 ** 6
    r_in_size = sys.getsizeof(compound["Reactant_in"])
    chunks, rem = divmod(r_in_size, max_size)

    if rem:
        chunks += 1
    chunk_size = ceil(len(compound["Reactant_in"]) / chunks)

    # Generate the InsertOne requests
    r_in_chunks = utils.Chunks(compound["Reactant_in"], chunk_size, return_list=True)

    requests = []
    for i, r_in_chunk in enumerate(r_in_chunks):
        requests.append(
            {
                "_id": f"{compound['_id']}_{i}",
                "c_id": compound["_id"],
                "Reactant_in": r_in_chunk,
            }
        )

    return requests


def _get_product_of_insert(compound: dict) -> List[dict]:
    """Write reactants_in, ensuring memory size isn't too big.

    MongoDB only allows < 16 MB entries. This function breaks large reactants_in
    up to ensure this doesn't happen.

    Parameters
    ----------
    compounds : List[dict]
        Dictionary of compounds to write.
    db : MINE
        MINE object to write compounds with.

    Returns
    -------
    List[dict]
        dicts of product_of to insert
    """

    # Get number of chunks product_of must be broken up into
    # 16 MB is the max for BSON, cut to 14 MB max just to be safe
    # Also some weirdness is that max_size must be 1 order of magnitude lower
    max_size = 1.4 * 10 ** 6
    p_of_size = sys.getsizeof(compound["Product_of"])
    chunks, rem = divmod(p_of_size, max_size)
    if rem:
        chunks += 1
    chunk_size = ceil(len(compound["Product_of"]) / chunks)

    # Generate the InsertOne requests
    p_of_chunks = utils.Chunks(compound["Product_of"], chunk_size, return_list=True)

    requests = []
    for i, p_of_chunk in enumerate(p_of_chunks):
        requests.append(
            {
                "_id": f"{compound['_id']}_{i}",
                "c_id": compound["_id"],
                "Product_of": p_of_chunk,
            }
        )

    return requests


# Core Compounds
def write_core_compounds(
    compounds: List[dict], db: MINE, mine: str, chunk_size: int = 10000, processes=1
) -> None:
    """Write core compounds to the core compound database.

    Calculates and formats compounds into appropriate form to insert into the
    core compound database in the mongo instance. Core compounds are attempted
    to be inserted and collisions are detected on the database. The list of
    MINEs a given compound is found in is updated as well.

    Parameters
    ----------
    compounds : dict
        List of compound dictionaries to write.
    db : MINE
        MINE object to write core compounds with.
    mine : str
        Name of the MINE.
    chunk_size : int, optional
        Size of chunks to break compounds into when writing, by default 10000.
    processes : int, optional
        The number of processors to use, by default 1.
    """
    n_cpds = len(compounds)
    if processes == 1:
        pool = None
    else:
        pool = multiprocessing.Pool(processes)

    for i, cpd_chunk in enumerate(utils.Chunks(compounds, chunk_size)):
        if i % 20 == 0:
            print(f"Writing Compounds: Chunk {i} of {int(n_cpds/chunk_size) + 1}")

        # Capture annoying RDKit output

        cpd_chunk = [deepcopy(cpd) for cpd in cpd_chunk if cpd["_id"].startswith("C")]
        if pool:
            core_requests = [req for req in pool.map(_get_core_cpd_insert, cpd_chunk)]
        else:
            core_requests = [req for req in map(_get_core_cpd_insert, cpd_chunk)]

        core_update_requests = [
            _get_core_cpd_update(cpd_dict, mine) for cpd_dict in cpd_chunk
        ]

        # Need to write update (i.e. what keeps track of which mines the
        # compound has been in)
        # first to ensure cpd exists
        db.core_compounds.bulk_write(core_requests)
        db.core_compounds.bulk_write(core_update_requests)
    if pool:
        pool.close()
        pool.join()


def _get_core_cpd_update(cpd_dict: dict, mine: str) -> pymongo.UpdateOne:
    return pymongo.UpdateOne({"_id": cpd_dict["_id"]}, {"$addToSet": {"MINES": mine}})


def _get_core_cpd_insert(cpd_dict: dict) -> pymongo.UpdateOne:
    """Generate core compound to be inserted"""
    core_keys = ["_id", "SMILES", "Inchi", "InchiKey", "Mass", "Formula"]
    core_dict = {
        key: cpd_dict.get(key) for key in core_keys if cpd_dict.get(key) != None
    }

    mol_object = AllChem.MolFromSmiles(core_dict["SMILES"])
    rdk_fp = [
        i
        for i, val in enumerate(list(AllChem.RDKFingerprint(mol_object, fpSize=512)))
        if val
    ]

    # Store all different representations of the molecule (SMILES, Formula,
    #  InChI key, etc.) as well as its properties in a dictionary
    if not "SMILES" in core_dict:
        core_dict["SMILES"] = AllChem.MolToSmiles(mol_object, True)
    if not "Inchi" in core_dict:
        core_dict["Inchi"] = AllChem.MolToInchi(mol_object)
    if not "Inchikey" in core_dict:
        core_dict["Inchikey"] = AllChem.InchiToInchiKey(core_dict["Inchi"])

    core_dict["Mass"] = AllChem.CalcExactMolWt(mol_object)
    core_dict["Charge"] = AllChem.GetFormalCharge(mol_object)
    core_dict["Formula"] = AllChem.CalcMolFormula(mol_object)
    core_dict["logP"] = AllChem.CalcCrippenDescriptors(mol_object)[0]
    core_dict["RDKit_fp"] = rdk_fp
    core_dict["len_RDKit_fp"] = len(rdk_fp)
    # core_dict['NP_likeness'] = nps.scoreMol(mol_object, nps_model)
    core_dict["Spectra"] = {}
    # Record which expansion it's coming from
    core_dict["MINES"] = []

    return pymongo.UpdateOne(
        {"_id": core_dict["_id"]}, {"$setOnInsert": core_dict}, upsert=True
    )


# Target Compounds
def write_targets_to_mine(
    targets: List[dict], db: MINE, chunk_size: int = 10000
) -> None:
    """Write target compounds to target collection of MINE.

    Parameters
    ----------
    targets : List[dict]
        Listt of target dictionaries to write.
    db : MINE
        MINE object to write targets with.
    chunk_size : int, optional
        Size of chunks to break compounds into when writing, by default 10000.
    """

    def _get_cpd_insert(cpd_dict: dict):
        output_keys = ["_id", "ID", "SMILES", "InChI_key"]
        return pymongo.InsertOne(
            {key: cpd_dict.get(key) for key in output_keys if cpd_dict.get(key) != None}
        )

    n_cpds = len(targets)
    for i, target_chunk in enumerate(utils.Chunks(targets, chunk_size)):
        if i % 20 == 0:
            print(f"Writing Targets: Chunk {i} of {int(n_cpds/chunk_size) + 1}")
        cpd_requests = [_get_cpd_insert(cpd_dict) for cpd_dict in target_chunk]
        db.target_compounds.bulk_write(cpd_requests, ordered=False)
