"""Databases.py: This file contains MINE database classes including database
loading functions."""
import ast
import datetime
import os
import sys
import platform
from copy import copy
from shutil import move
from subprocess import call

import pymongo
from pymongo.errors import ServerSelectionTimeoutError
from rdkit.Chem import AllChem

from minedatabase import utils
from minedatabase.NP_Score import npscorer as nps

nps_model = nps.readNPModel()

def establish_db_client(con_string=None):
    """This establishes a mongo database client in various environments"""
    # If a connection string is given use that, otherwise go to local
    try:
        if con_string:
            client = pymongo.MongoClient(con_string, ServerSelectionTimeoutMS=5)
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
    def __init__(self, name, con_string='mongodb://localhost:27017/'):
        self.client = establish_db_client(con_string)
        self.con_string = con_string
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
        self.nps_model = nps.readNPModel()
        self._mass_cache = {}  # for rapid calculation of reaction mass change

    def add_rxn_pointers(self):
        """Add links to the reactions that each compound participates in
        allowing for users to follow paths in the network"""
        reactions_count = self.reactions.count()
        print(f"Linking compounds to {reactions_count} reactions")
        for reaction in self.reactions.find().batch_size(500):
            # Update pointers for compounds that are reactants
            for compound in reaction['Reactants']:
                self.compounds.update({'_id': compound['c_id']},
                                      {'$push': {
                                          'Reactant_in': reaction['_id']}})
            # Update pointers for compounds that are products
            for compound in reaction['Products']:
                self.compounds.update({'_id': compound['c_id']},
                                      {'$push': {
                                          'Product_of': reaction['_id']}})
        # Write to log file
        self.meta_data.insert({'Timestamp': datetime.datetime.now(), 'Action':
                               "Add Reaction Pointers"})

    def add_compound_sources(self, rxn_key_type='_id'):
        """Adds a field detailing the compounds and reactions from which a
        compound is produced. This enables filtering search results to only
        include compounds which are relevant to a specified metabolic context.
        This is different from add_rxn_pointers in that this provides the
        precursor compounds rather than the reactions that the compound
        participates in."""
        for compound in self.compounds.find({'Sources': {'$exists': 0}}):
            compound['Sources'] = []

            for reaction in self.reactions.find(
                    {'Products.c_id': compound[rxn_key_type]}):
                compound['Sources'].append(
                    {'Compounds': [x['c_id'] for x in reaction['Reactants']],
                     'Operators': reaction['Operators']})
            # If there are sources, then save them and make sure there aren't
            #  too many for a single compound.
            if compound['Sources']:
                try:
                    self.compounds.save(compound)
                except pymongo.errors.DocumentTooLarge:
                    print(f"Too Many Sources for {compound['SMILES']}")
        # Allow for quick searching within database
        self.compounds.ensure_index([('Sources.Compound', pymongo.ASCENDING),
                                     ('Sources.Operators', pymongo.ASCENDING)])
        # Write to log file
        self.meta_data.insert({'Timestamp': datetime.datetime.now(), 'Action':
                               "Add Compound Source field"})

    def add_reaction_mass_change(self, reaction=None):
        """Calculate the change in mass between reactant and product compounds.
            This is useful for discovering compounds in molecular networking
        :param reaction: If specified, the function returns the mass change for
            the supplied reaction. If false, it calculates the mass change for
            all reactions missing this value in the database.
        :type reaction: str
        :return: If reaction specified, the mass change
        :rtype: float
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

    def generate_image_files(self, path, query=None, dir_depth=0,
                             img_type='svg:-a,nosource,w500,h500',
                             convert_r=False):
        """Generates image files for compounds in database using ChemAxon's
        MolConvert.

        :param path: Target directory for image files
        :type path: str
        :param query: A query to limit the number of files generated
        :type query: dict
        :param dir_depth: The number of directory levels to split the
            compounds into for files system efficiency. Ranges from 0 (all in
            top level directory to the length of the file name (40 for MINE
            hashes)
        :type dir_depth: int
        :param img_type: The type of image file to be generated. See molconvert
            documentation for valid options
        :type img_type: str

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

    def build_indexes(self):
        """Builds indexes for efficient querying of MINE databases"""
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

    def link_to_external_database(self, external_database, compound=None,
                                  match_field='Inchikey', fields_to_copy=None):
        """This function looks for matching compounds in other databases (i.e.
        PubChem) and adds links where found.

        :param external_database: The name of the database to search for
            matching compounds
        :type external_database: str
        :param compound: The compound to search for external links. If none,
            link all compounds in the database.
        :type compound: dict
        :param match_field: The field to search on for matching compounds
        :type match_field: str
        :param fields_to_copy: Data to copy into the mine database. The first
            field is the field name in the external database. The second field
            is the field name in the MINE database where the data will be
            copied.
        :type fields_to_copy: list(tuple)
        """
        if compound:
            ext = MINE(external_database)
            projection = dict([('_id', 0,)] + [(x[0], 1,) for x
                                               in fields_to_copy])
            # Find compounds that have same name in another database
            for ext_comp in ext.compounds.find(
                    {match_field: compound[match_field]}, projection):
                for field in fields_to_copy:
                    if field[0] in ext_comp:
                        # dict_merge merges two dictionaries using sets to
                        # avoid duplicate values
                        utils.dict_merge(compound, utils.save_dotted_field(
                            field[1], utils.get_dotted_field(ext_comp,
                                                             field[0])))
            return utils.convert_sets_to_lists(compound)

        # If compound is None, link all compounds in database
        else:
            for comp in self.compounds.find():
                self.compounds.save(self.link_to_external_database(
                    external_database, compound=comp, match_field=match_field,
                    fields_to_copy=fields_to_copy))


    def insert_core_compound(self, compound_dict, requests=None):
        """This method generates a mongo request to save a compound into the core database.
        The necessary fields for the API are calculated.
        If a list of requests are given the request is appended for later bulk writing.
        Otherwise a single entry is made. If a compound is already in the core database
        nothing is written.

        :param compound_dict: Compound Dictionary
        :type compound_dict: dict
        :param requests: List of requests for bulk insert
        :type requests: None
        """
        core_dict = copy(compound_dict)
        cpd_id = core_dict['_id']
        mol_object = AllChem.MolFromSmiles(core_dict['SMILES'])

        if 'Generation' in core_dict:
            del(core_dict['Generation'])
        if 'Expand' in core_dict:
            del(core_dict['Expand'])
        if 'Type' in core_dict:
            del(core_dict['Type'])
        if 'Product_of' in core_dict:
            del(core_dict['Product_of'])
        if 'Reactant_in' in core_dict:
            del(core_dict['Reactant_in'])
        # Store all different representations of the molecule (SMILES, Formula,
        #  InChI key, etc.) as well as its properties in a dictionary
        if not 'SMILES' in core_dict:
            core_dict['SMILES'] = AllChem.MolToSmiles(mol_object, True)
        if not 'Inchi' in core_dict:
            core_dict['Inchi'] = AllChem.MolToInchi(mol_object)
        if not 'Inchikey' in core_dict:
            core_dict['Inchikey'] = AllChem.InchiToInchiKey(
                core_dict['Inchi'])
        core_dict['Mass'] = AllChem.CalcExactMolWt(mol_object)
        core_dict['Formula'] = AllChem.CalcMolFormula(mol_object)
        core_dict['logP'] = AllChem.CalcCrippenDescriptors(mol_object)[0]
        core_dict['NP_likeness'] = nps.scoreMol(mol_object, self.nps_model)
        core_dict['Spectra'] = {}
        # Record which expansion it's coming from
        core_dict['MINES'] = []

        if requests != None:
            requests.append(pymongo.UpdateOne({'_id': cpd_id}, {'$setOnInsert': core_dict}, upsert=True))
        else:
            self.core_compounds.update_one({'_id': cpd_id}, {'$setOnInsert': core_dict}, upsert=True)
            
        return None
    
    def update_core_compound_MINES(self, compound_dict, requests=None):
        """This method adds the current MINE to the core compound's MINES set.

        :param compound_dict: RDKit object of the compound
        :type compound_dict: dict
        :param requests: List of requests for bulk insert
        :type requests: None
        """
        cpd_id = compound_dict['_id']
        if requests != None:
            requests.append(pymongo.UpdateOne({'_id': cpd_id}, {'$addToSet': {'MINES': self.name}}))
        else:
            self.core_compounds.update_one({'_id': cpd_id}, {'$addToSet': {'MINES': self.name}})

        return None

    def insert_mine_compound(self, compound_dict=None, requests=None):
        """This method saves a RDKit Molecule as a compound entry in the MINE.
        Calculates necessary fields for API and includes additional
        information passed in the compound dict. Overwrites preexisting
        compounds in MINE on _id collision.

        :param compound_dict: Information of compound to insert to database
        :type compound_dict: dict
        :param requests: A list of requests for pymongo bulk write
        : type requests: list
        """

        if compound_dict is None:
            return None

        # Store all different representations of the molecule (SMILES, Formula,
        #  InChI key, etc.) as well as its properties in a dictionary
        if '_atom_count' in compound_dict:
            del compound_dict['_atom_count']

        if 'Inchikey' in compound_dict:
            del compound_dict['Inchikey']

        if 'ID' in compound_dict:
            del compound_dict['ID']
        # If the compound is a reactant, then make sure the reactant name is
        # in a correct format.
        if 'Reactant_in' in compound_dict and isinstance(
                compound_dict['Reactant_in'], str) \
                and compound_dict['Reactant_in']:
            compound_dict['Reactant_in'] = ast.literal_eval(
                compound_dict['Reactant_in'])
        # If the compound is a product, then make sure the reactant name is
        # in a correct format.
        if 'Product_of' in compound_dict \
                and isinstance(compound_dict['Product_of'], str) \
                and compound_dict['Product_of']:
            compound_dict['Product_of'] = ast.literal_eval(
                compound_dict['Product_of'])

        # If bulk insertion, upsert (insert and update) the database
        if requests != None:
            requests.append(pymongo.ReplaceOne({'_id':compound_dict['_id']}, compound_dict, upsert=True))
        else:
            self.compounds.replace_one({'_id':compound_dict['_id']}, compound_dict, upsert=True)

        return None

    def insert_reaction(self, reaction_dict, requests=None):
        """Inserts a reaction into the MINE database and returns _id of the
         reaction in the mine database.

        :param reaction_dict: A dictionary containing 'Reactants' and
         'Products' lists of StoichTuples
        :type reaction_dict: dict        
        :return: The hashed _id of the reaction
        :rtype: str
        """

        # By converting to a dict, mongo stores the data as objects not
        # arrays allowing for queries by compound hash
        # if isinstance(reaction_dict['Reactants'][0], utils.StoichTuple):
        #     reaction_dict['Reactants'] = [x._asdict() for x
        #                                   in reaction_dict['Reactants']]
        #     reaction_dict['Products'] = [x._asdict() for x
        #                                  in reaction_dict['Products']]

        reaction_dict = utils.convert_sets_to_lists(reaction_dict)
        # If bulk insertion, upsert (insert and update) the database
        if requests != None:
            requests.append(pymongo.ReplaceOne({'_id':{'$eq':reaction_dict['_id']}}, reaction_dict, upsert=True))
        else:
            save_document(self.reactions, reaction_dict)
        return reaction_dict['_id']

    def map_reactions(self, ext_db, match_field='_id'):
        """Update operators by adding the reactions they participate in.

        :param ext_db: name of MongoDB
        :param match_field: name of reaction field to match based on
        """
        # Update operators by adding the reactions they participate in
        lit_db = MINE(ext_db)
        for lit_rxn in lit_db.reactions.find():
            if match_field == 'InChI_hash':
                mine_rxn = self.reactions.find_one(
                    {match_field: {'$regex': lit_rxn[match_field].split('-')[0]
                                   }}, {'Operators': 1})
            else:
                mine_rxn = self.reactions.find_one(
                    {match_field: lit_rxn[match_field]}, {'Operators': 1})
            if mine_rxn:
                for op in mine_rxn['Operators']:
                    self.operators.update(
                        {'_id': op},
                        {'$addToSet': {'Mapped_Rxns': lit_rxn['_id'],
                                       'References': {
                                           '$each': lit_rxn['References']}}})
                lit_db.reactions.update(
                    {'_id': lit_rxn['_id']},
                    {'$set': {'Mapped_Rules': mine_rxn['Operators']}})


def save_document(collection, doc):
    """Saves document to given MongoDB collection.

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

def insert_mine_compound(compound_dict=None):
        """This method saves a RDKit Molecule as a compound entry in the MINE.
        Calculates necessary fields for API and includes additional
        information passed in the compound dict.

        :param compound_dict: Compound info to add to database
        :type compound_dict: dict
        :return: The mongo request
        :rtype: pymongo.InsertOne
        """

        # Store all different representations of the molecule (SMILES, Formula,
        #  InChI key, etc.) as well as its properties in a dictionary
        if '_atom_count' in compound_dict:
            del compound_dict['_atom_count']

        if 'Inchikey' in compound_dict:
            del compound_dict['Inchikey']

        if 'ID' in compound_dict:
            del compound_dict['ID']
        # If the compound is a reactant, then make sure the reactant name is
        # in a correct format.
        if 'Reactant_in' in compound_dict and isinstance(
                compound_dict['Reactant_in'], str) \
                and compound_dict['Reactant_in']:
            compound_dict['Reactant_in'] = ast.literal_eval(
                compound_dict['Reactant_in'])
        # If the compound is a product, then make sure the reactant name is
        # in a correct format.
        if 'Product_of' in compound_dict \
                and isinstance(compound_dict['Product_of'], str) \
                and compound_dict['Product_of']:
            compound_dict['Product_of'] = ast.literal_eval(
                compound_dict['Product_of'])

        return pymongo.InsertOne(compound_dict)

def insert_core_compound(compound_dict):
        """This method generates a mongo request to save a compound into the core database.
        The necessary fields for the API are calculated.
        If a list of requests are given the request is appended for later bulk writing.
        Otherwise a single entry is made. If a compound is already in the core database
        nothing is written.

        :param compound_dict: Compound Dictionary
        :type compound_dict: dict
        :return requests: List of requests for bulk insert
        :rtype requests: pymongo.UpdateOne
        """
        core_dict = copy(compound_dict)
        cpd_id = core_dict['_id']
        mol_object = AllChem.MolFromSmiles(core_dict['SMILES'])

        if 'Generation' in core_dict:
            del(core_dict['Generation'])
        if 'Expand' in core_dict:
            del(core_dict['Expand'])
        if 'Type' in core_dict:
            del(core_dict['Type'])
        if 'Product_of' in core_dict:
            del(core_dict['Product_of'])
        if 'Reactant_in' in core_dict:
            del(core_dict['Reactant_in'])
                     
        # Store all different representations of the molecule (SMILES, Formula,
        #  InChI key, etc.) as well as its properties in a dictionary
        if not 'SMILES' in core_dict:
            core_dict['SMILES'] = AllChem.MolToSmiles(mol_object, True)
        if not 'Inchi' in core_dict:
            core_dict['Inchi'] = AllChem.MolToInchi(mol_object)
        if not 'Inchikey' in core_dict:
            core_dict['Inchikey'] = AllChem.InchiToInchiKey(
                core_dict['Inchi'])
        core_dict['Mass'] = AllChem.CalcExactMolWt(mol_object)
        core_dict['Formula'] = AllChem.CalcMolFormula(mol_object)
        core_dict['logP'] = AllChem.CalcCrippenDescriptors(mol_object)[0]
        core_dict['NP_likeness'] = nps.scoreMol(mol_object, nps_model)
        core_dict['Spectra'] = {}
        # Record which expansion it's coming from
        core_dict['MINES'] = []   
        
        return pymongo.UpdateOne({'_id': cpd_id}, {'$setOnInsert': core_dict}, upsert=True)

def update_core_compound_MINES(compound_dict, MINE):
        """This method adds the current MINE to the core compound's MINES set.

        :param compound_dict: RDKit object of the compound
        :type compound_dict: dict
        :return requests: Request for bulk insert
        :rtype requests: pymongo.UpdateOne
        """
        cpd_id = compound_dict['_id']
        return pymongo.UpdateOne({'_id': cpd_id}, {'$addToSet': {'MINES': MINE}})

def insert_reaction(reaction_dict):
    """Inserts a reaction into the MINE database and returns _id of the
        reaction in the mine database.

    :param reaction_dict: A dictionary containing 'Reactants' and
        'Products' lists of StoichTuples
    :type reaction_dict: dict    
    :return: Request for bulk insert
    :rtype: pymongo.InsertOne
    """
    reaction_dict['_id'] = utils.rxn2hash(reaction_dict['Reactants'],
                                            reaction_dict['Products'])

    # By converting to a dict, mongo stores the data as objects not
    # arrays allowing for queries by compound hash
    if isinstance(reaction_dict['Reactants'][0], utils.StoichTuple):
        reaction_dict['Reactants'] = [x._asdict() for x
                                        in reaction_dict['Reactants']]
        reaction_dict['Products'] = [x._asdict() for x
                                        in reaction_dict['Products']]

    reaction_dict = utils.convert_sets_to_lists(reaction_dict)

    return pymongo.InsertOne(reaction_dict)