"""Databases.py: This file contains MINE database classes including database
loading and writing functions."""
import ast
import datetime
import multiprocessing
import os
import platform
import sys
from copy import deepcopy
from shutil import move
from subprocess import call
from typing import List

import pymongo
from pymongo.errors import ServerSelectionTimeoutError
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger

from minedatabase import utils

# from minedatabase import utils
# from minedatabase.NP_Score import npscorer as nps

# from minedatabase.NP_Score import npscorer as nps

# nps_model = nps.readNPModel()

lg = logger()
lg.setLevel(4)

def establish_db_client(uri: str=None)->pymongo.MongoClient:
    """ Establish a connection to a mongo database given a URI

    Uses the provided URI to connect to a mongoDB. If none is given
    the default URI is used when using pymongo.

    Parameters
    ----------
    uri : [type], optional
        URI to connect to mongo DB, by default None

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
        raise IOError("Failed to load database client. Please verify that "
                      "mongod is running")
    return client


class MINE:
    """
    This class basically exposes the underlying mongo database to manipulation
    but also defines expected database structure.
    """
    def __init__(self, name: str, uri: str='mongodb://localhost:27017/'):
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
        # self.nps_model = nps.readNPModel()
        self._mass_cache = {}  # for rapid calculation of reaction mass change

    # def add_rxn_pointers(self):
    #     """Add links to the reactions that each compound participates in
    #     allowing for users to follow paths in the network"""
    #     reactions_count = self.reactions.count()
    #     print(f"Linking compounds to {reactions_count} reactions")
    #     for reaction in self.reactions.find().batch_size(500):
    #         # Update pointers for compounds that are reactants
    #         for compound in reaction['Reactants']:
    #             self.compounds.update({'_id': compound['c_id']},
    #                                   {'$push': {
    #                                       'Reactant_in': reaction['_id']}})
    #         # Update pointers for compounds that are products
    #         for compound in reaction['Products']:
    #             self.compounds.update({'_id': compound['c_id']},
    #                                   {'$push': {
    #                                       'Product_of': reaction['_id']}})
    #     # Write to log file
    #     self.meta_data.insert({'Timestamp': datetime.datetime.now(), 'Action':
    #                            "Add Reaction Pointers"})

    # def add_compound_sources(self, rxn_key_type='_id'):
        # """Adds a field detailing the compounds and reactions from which a
        # compound is produced. This enables filtering search results to only
        # include compounds which are relevant to a specified metabolic context.
        # This is different from add_rxn_pointers in that this provides the
        # precursor compounds rather than the reactions that the compound
        # participates in."""
        # for compound in self.compounds.find({'Sources': {'$exists': 0}}):
        #     compound['Sources'] = []

        #     for reaction in self.reactions.find(
        #             {'Products.c_id': compound[rxn_key_type]}):
        #         compound['Sources'].append(
        #             {'Compounds': [x['c_id'] for x in reaction['Reactants']],
        #              'Operators': reaction['Operators']})
        #     # If there are sources, then save them and make sure there aren't
        #     #  too many for a single compound.
        #     if compound['Sources']:
        #         try:
        #             self.compounds.save(compound)
        #         except pymongo.errors.DocumentTooLarge:
        #             print(f"Too Many Sources for {compound['SMILES']}")
        # # Allow for quick searching within database
        # self.compounds.ensure_index([('Sources.Compound', pymongo.ASCENDING),
        #                              ('Sources.Operators', pymongo.ASCENDING)])
        # # Write to log file
        # self.meta_data.insert({'Timestamp': datetime.datetime.now(), 'Action':
        #                        "Add Compound Source field"})

    def add_reaction_mass_change(self, reaction: str=None) -> None:
        """Calculate the change in mass between reactant and product compounds.

        This is useful for discovering compounds in molecular networking.

        Parameters
        ----------
        reaction : str, optional
            Reaction ID to calculate the mass change for, by default None

        Returns
        -------
        float
            Mass change of specified reaction

        """
        def _get_mass(_id):
            if _id not in self._mass_cache:
                try:
                    self._mass_cache['_id'] = self.compounds.find_one(
                        {'_id': _id}, {'Mass': 1})['Mass']
                except (TypeError, KeyError):
                    raise ValueError(f'A mass value for {_id} was not found in '
                                     'compounds collection')
            return self._mass_cache['_id']

        def _get_mass_change(rxn):
            if isinstance(rxn['Reactants'][0], dict):
                key = 'c_id'
            else:
                key = 1
            start_mass = _get_mass(rxn['Reactants'][0][key])
            return [_get_mass(product[key]) - start_mass for product
                    in rxn['Products'] if product[key][0] == 'C']

        if reaction:
            return _get_mass_change(reaction)

        bulk = self.reactions.initialize_unordered_bulk_op()
        reactions = self.reactions.find({'Mass_Change': {'$exists': 0}})
        for rxn in reactions:
            try:
                bulk.find({'_id': rxn['_id']}).update({'$set': {
                    'Mass_Change': _get_mass_change(rxn)}})
            except ValueError as e:
                print(e.args[0])
                print(f"Failed to parse rxn {rxn['_id']}")
        bulk.execute()


    def generate_image_files(self, path: str, query: dict=None, dir_depth: int=0,
                             img_type: str='svg:-a,nosource,w500,h500',
                             convert_r: bool=False) -> None:
        """Generates image files for compounds in database using ChemAxon's
        MolConvert.

        Parameters
        ----------
        path : str
            Target directory for image file
        query : dict, optional
            Query to limit number of files generated, by default None
        dir_depth : int, optional
            The number of directory levels to split the compounds
            into for files system efficiency. Ranges from 0 (all in top
            level directory) to the length of the file name (40 for MINE hashes), by default 0
        img_type : str, optional
            Type of image file to be generated. See molconvert
            documentation for valid options, by default 'svg:-a,nosource,w500,h500'
        convert_r : bool, optional
            Convert R in the smiles to *, by default False
        """
        ids = []
        extension = img_type.split(':')[0]
        structure_file = os.path.join(path, 'tmp.smiles')
        if not query:
            query = {}

        if not os.path.exists(path):
            if sys.platform == 'linux':
                os.system(f"sudo mkdir {path} -m 777")
            else:
                os.mkdir(path)
        with open(structure_file, 'w') as outfile:
            for comp in self.compounds.find(query, {'SMILES': 1}):
                if convert_r:
                    outfile.write(f"{comp['SMILES'].replace('R', '*')}\n")
                else:
                    outfile.write(f"{comp['SMILES']}\n")
                ids.append(comp['_id'])
        if platform.system() == 'Windows':
            rc = call([f"molconvert -mo '{path}/.{extension}' {img_type} '{structure_file}'"],
                      shell=True)
        else:
            rc = call([f"molconvert -mo '{path}/.{extension}' {img_type} '{structure_file}'"],
                      shell=True)
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
            new = os.path.join(new, _id + '.' + extension)
            if os.path.isfile(old):
                move(old, new)

    def build_indexes(self) -> None:
        """Build indexes for efficient querying of the database"""
        self.core_compounds.drop_indexes()
        self.core_compounds.create_index([('Mass', pymongo.ASCENDING)])
        self.core_compounds.create_index('Inchikey')
        self.compounds.create_index('Inchikey')
        self.compounds.create_index('Inchi')
        self.compounds.create_index('SMILES')
        self.reactions.create_index('Reactants.c_id')
        self.reactions.create_index('Products.c_id')
        self.meta_data.insert_one({'Timestamp': datetime.datetime.now(),
                                   'Action': "Database indexes built"})

    # def link_to_external_database(self, external_database, compound=None,
    #                               match_field='Inchikey', fields_to_copy=None):
    #     """This function looks for matching compounds in other databases (i.e.
    #     PubChem) and adds links where found.

    #     :param external_database: The name of the database to search for
    #         matching compounds
    #     :type external_database: str
    #     :param compound: The compound to search for external links. If none,
    #         link all compounds in the database.
    #     :type compound: dict
    #     :param match_field: The field to search on for matching compounds
    #     :type match_field: str
    #     :param fields_to_copy: Data to copy into the mine database. The first
    #         field is the field name in the external database. The second field
    #         is the field name in the MINE database where the data will be
    #         copied.
    #     :type fields_to_copy: list(tuple)
    #     """
    #     if compound:
    #         ext = MINE(external_database)
    #         projection = dict([('_id', 0,)] + [(x[0], 1,) for x
    #                                            in fields_to_copy])
    #         # Find compounds that have same name in another database
    #         for ext_comp in ext.compounds.find(
    #                 {match_field: compound[match_field]}, projection):
    #             for field in fields_to_copy:
    #                 if field[0] in ext_comp:
    #                     # dict_merge merges two dictionaries using sets to
    #                     # avoid duplicate values
    #                     utils.dict_merge(compound, utils.save_dotted_field(
    #                         field[1], utils.get_dotted_field(ext_comp,
    #                                                          field[0])))
    #         return utils.convert_sets_to_lists(compound)

    #     # If compound is None, link all compounds in database
    #     else:
    #         for comp in self.compounds.find():
    #             self.compounds.save(self.link_to_external_database(
    #                 external_database, compound=comp, match_field=match_field,
    #                 fields_to_copy=fields_to_copy))

    # def map_reactions(self, ext_db, match_field='_id'):
    #     """Update operators by adding the reactions they participate in.

    #     :param ext_db: name of MongoDB
    #     :param match_field: name of reaction field to match based on
    #     """
    #     # Update operators by adding the reactions they participate in
    #     lit_db = MINE(ext_db)
    #     for lit_rxn in lit_db.reactions.find():
    #         if match_field == 'InChI_hash':
    #             mine_rxn = self.reactions.find_one(
    #                 {match_field: {'$regex': lit_rxn[match_field].split('-')[0]
    #                                }}, {'Operators': 1})
    #         else:
    #             mine_rxn = self.reactions.find_one(
    #                 {match_field: lit_rxn[match_field]}, {'Operators': 1})
    #         if mine_rxn:
    #             for op in mine_rxn['Operators']:
    #                 self.operators.update(
    #                     {'_id': op},
    #                     {'$addToSet': {'Mapped_Rxns': lit_rxn['_id'],
    #                                    'References': {
    #                                        '$each': lit_rxn['References']}}})
    #             lit_db.reactions.update(
    #                 {'_id': lit_rxn['_id']},
    #                 {'$set': {'Mapped_Rules': mine_rxn['Operators']}})


def save_document(collection: pymongo.collection.Collection, doc) -> None:
    """Save document to given MongoDB collection.

    If document _id already exists, replaces it. Else, inserts it.

    Parameters
    ----------
    collection : mongodb.collection
        Collection to save document to.
    doc : dict
        Document to save.
    """
    if '_id' in doc and collection.find_one({'_id': doc['_id']}):
        collection.replace_one({'_id': doc['_id']}, doc)
    else:
        collection.insert_one(doc)
# Functions to write data to MINE

# Reactions
def write_reactions_to_mine(reactions: List[dict], db: MINE, chunk_size: int = 10000) -> None:
    """Write reactions to reaction collection of MINE.

    Parameters
    ----------
    reactions : List[dict]
        Dictionary of reactions to write.
    db : MINE
        MINE object to write reactions with.
    chunk_size : int, optional
        Size of chunks to break reactions into when writing, by default 10000
    """
    n_rxns = len(reactions)
    for i, rxn_chunk in enumerate(utils.Chunks(reactions, chunk_size)):
        if i%20 == 0:
            print(f"Writing Reactions: Chunk {i} of {int(n_rxns/chunk_size) + 1}")
        rxn_requests = [pymongo.InsertOne(utils.convert_sets_to_lists(rxn_dict)) for rxn_dict in rxn_chunk]

        db.reactions.bulk_write(rxn_requests, ordered=False)


# Compounds
def write_compounds_to_mine(compounds: List[dict], db: MINE, chunk_size: int = 10000) -> None:
    """Write compounds to reaction collection of MINE.

    Parameters
    ----------
    compounds : List[dict]
        Dictionary of compounds to write.
    db : MINE
        MINE object to write compounds with.
    chunk_size : int, optional
        Size of chunks to break compounds into when writing, by default 10000
    """
    def _get_cpd_insert(cpd_dict: dict):
        output_keys = ["_id", "ID", "SMILES", "InChi_key", "Type", "Generation", "Reactant_in",
                       "Product_of", "Expand", "Matched_Peak_IDs", "Matched_Adducts"]
        return pymongo.InsertOne({key: cpd_dict.get(key) for key in output_keys if cpd_dict.get(key) != None})

    n_cpds = len(compounds)
    for i, cpd_chunk in enumerate(utils.Chunks(compounds, chunk_size)):
        if i%20 == 0:
            print(f"Writing Compounds: Chunk {i} of {int(n_cpds/chunk_size) + 1}")
        cpd_requests = [_get_cpd_insert(cpd_dict) for cpd_dict in cpd_chunk]
        db.compounds.bulk_write(cpd_requests, ordered=False)


# Core Compounds
def write_core_compounds(compounds: List[dict], db: MINE, mine: str, chunk_size: int = 10000, processes = 1) -> None:
    """Write core compounds to the core compound database

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
        Size of chunks to break compounds into when writing, by default 10000
    processes : int, optional
        The number of processors to use, by default 1
    """
    n_cpds = len(compounds)
    pool = multiprocessing.Pool(processes)

    for i, cpd_chunk in enumerate(utils.Chunks(compounds, chunk_size)):
        if i%20 == 0:
            print(f"Writing Compounds: Chunk {i} of {int(n_cpds/chunk_size) + 1}")

        # Capture annoying RDKit output

        cpd_chunk = [deepcopy(cpd) for cpd in cpd_chunk if cpd["_id"].startswith("C")]
        core_requests = [req for req in pool.map(_get_core_cpd_insert, cpd_chunk)]
        core_update_requests = [_get_core_cpd_update(cpd_dict, mine) for cpd_dict in cpd_chunk]

        # Need to write update (i.e. what keeps track of which mines the compound has been in)
        # first to ensure cpd exists
        db.core_compounds.bulk_write(core_requests)
        db.core_compounds.bulk_write(core_update_requests)

    pool.close()


def _get_core_cpd_update(cpd_dict: dict, mine: str) -> pymongo.UpdateOne:
    return pymongo.UpdateOne({'_id': cpd_dict['_id']}, {'$addToSet': {'MINES': mine}})


def _get_core_cpd_insert(cpd_dict: dict) -> pymongo.InsertOne:

    # cpd_dict = deepcopy(cpd_dict)
    core_keys = ["_id", "SMILES", "Inchi", "InchiKey", "Mass", "Formula"]
    core_dict = {key: cpd_dict.get(key) for key in core_keys if cpd_dict.get(key) != None}

    mol_object = AllChem.MolFromSmiles(core_dict['SMILES'])

    # Store all different representations of the molecule (SMILES, Formula,
    #  InChI key, etc.) as well as its properties in a dictionary
    if not 'SMILES' in core_dict:
        core_dict['SMILES'] = AllChem.MolToSmiles(mol_object, True)
    if not 'Inchi' in core_dict:
        core_dict['Inchi'] = AllChem.MolToInchi(mol_object)
    if not 'Inchikey' in core_dict:
        core_dict['Inchikey'] = AllChem.InchiToInchiKey(core_dict['Inchi'])

    core_dict['Mass'] = AllChem.CalcExactMolWt(mol_object)
    core_dict['Formula'] = AllChem.CalcMolFormula(mol_object)
    core_dict['logP'] = AllChem.CalcCrippenDescriptors(mol_object)[0]
    # core_dict['NP_likeness'] = nps.scoreMol(mol_object, nps_model)
    core_dict['Spectra'] = {}
    # Record which expansion it's coming from
    core_dict['MINES'] = []

    return pymongo.UpdateOne({'_id': core_dict["_id"]}, {'$setOnInsert': core_dict}, upsert=True)

# Target Compounds
def write_targets_to_mine(targets: List[dict], db: MINE, chunk_size: int = 10000) -> None:
    """Write target compounds to target collection of MINE.

    Parameters
    ----------
    targets : List[dict]
        Listt of target dictionaries to write.
    db : MINE
        MINE object to write targets with.
    chunk_size : int, optional
        [description], by default 10000
    """
    def _get_cpd_insert(cpd_dict: dict):
        output_keys = ["_id", "ID", "SMILES", "InChi_key"]
        return pymongo.InsertOne({key: cpd_dict.get(key) for key in output_keys if cpd_dict.get(key) != None})

    n_cpds = len(targets)
    for i, target_chunk in enumerate(utils.Chunks(targets, chunk_size)):
        if i%20 == 0:
            print(f"Writing Targets: Chunk {i} of {int(n_cpds/chunk_size) + 1}")
        cpd_requests = [_get_cpd_insert(cpd_dict) for cpd_dict in target_chunk]
        db.target_compounds.bulk_write(cpd_requests, ordered=False)
