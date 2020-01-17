"""Databases.py: This file contains MINE database classes including database
loading functions."""
import ast
import datetime
import os
import platform
from shutil import move
from subprocess import call

import pymongo
from pymongo.errors import ServerSelectionTimeoutError
from rdkit.Chem import AllChem

from minedatabase import utils
from minedatabase.NP_Score import npscorer as np


def establish_db_client():
    """This establishes a mongo database client in various environments"""
    try:
        # Special case for working on a sapphire node
        if 'node' in platform.node():
            client = pymongo.MongoClient(host='master',
                                         serverSelectionTimeoutMS=5)
        # Special case for working on a SEED cluster
        elif 'bio' in platform.node() or platform.node() == 'twig' \
                or platform.node() == 'branch':
            client = pymongo.MongoClient(host='branch',
                                         serverSelectionTimeoutMS=5)
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        # Local database
        else:
            client = pymongo.MongoClient(serverSelectionTimeoutMS=5)
    except ServerSelectionTimeoutError:
        raise IOError("Failed to load database client. Please verify that "
                      "mongod is running")
    return client


class MINE:
    """
    This class basically exposes the underlying mongo database to manipulation
    but also defines expected database structure.
    """
    def __init__(self, name):
        client = establish_db_client()
        db = client[name]
        self._db = db
        self.name = name
        self.meta_data = db.meta_data
        self.compounds = db.compounds
        self.reactions = db.reactions
        self.operators = db.operators
        self.models = db.models
        self.np_model = None
        self.id_db = client['UniversalMINE']
        self._mass_cache = {}  # for rapid calculation of reaction mass change

    def add_rxn_pointers(self):
        """Add links to the reactions that each compound participates in
        allowing for users to follow paths in the network"""
        reactions_count = self.reactions.count()
        print("Linking compounds to %s reactions" % reactions_count)
        for reaction in self.reactions.find().batch_size(500):
            # Update pointers for compounds that are reactants
            for compound in reaction['Reactants']:
                self.compounds.update({"_id": compound["c_id"]},
                                      {'$push': {
                                          "Reactant_in": reaction['_id']}})
            # Update pointers for compounds that are products
            for compound in reaction['Products']:
                self.compounds.update({"_id": compound["c_id"]},
                                      {'$push': {
                                          "Product_of": reaction['_id']}})
        # Write to log file
        self.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action":
                               "Add Reaction Pointers"})

    def add_compound_sources(self, rxn_key_type="_id"):
        """Adds a field detailing the compounds and reactions from which a
        compound is produced. This enables filtering search results to only
        include compounds which are relevant to a specified metabolic context.
        This is different from add_rxn_pointers in that this provides the
        precursor compounds rather than the reactions that the compound
        participates in."""
        for compound in self.compounds.find({"Sources": {"$exists": 0}}):
            compound['Sources'] = []

            for reaction in self.reactions.find(
                    {"Products.c_id": compound[rxn_key_type]}):
                compound['Sources'].append(
                    {"Compounds": [x['c_id'] for x in reaction['Reactants']],
                     "Operators": reaction["Operators"]})
            # If there are sources, then save them and make sure there aren't
            #  too many for a single compound.
            if compound['Sources']:
                try:
                    self.compounds.save(compound)
                except pymongo.errors.DocumentTooLarge:
                    print("Too Many Sources for %s" % compound['SMILES'])
        # Allow for quick searching within database
        self.compounds.ensure_index([("Sources.Compound", pymongo.ASCENDING),
                                     ("Sources.Operators", pymongo.ASCENDING)])
        # Write to log file
        self.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action":
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
                        {"_id": _id}, {'Mass': 1})['Mass']
                except (TypeError, KeyError):
                    raise ValueError('A mass value for %s was not found in '
                                     'compounds collection' % _id)
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
        reactions = self.reactions.find({"Mass_Change": {'$exists': 0}})
        for rxn in reactions:
            try:
                bulk.find({'_id': rxn['_id']}).update({'$set': {
                    'Mass_Change': _get_mass_change(rxn)}})
            except ValueError as e:
                print(e.args[0])
                print('Failed to parse rxn %s' % rxn['_id'])
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
        extension = img_type.split(":")[0]
        structure_file = os.path.join(path, 'tmp.smiles')
        if not query:
            query = {}

        if not os.path.exists(path):
            os.mkdir(path)
        with open(structure_file, 'w') as outfile:
            for comp in self.compounds.find(query, {'SMILES': 1}):
                if convert_r:
                    outfile.write("%s\n" % comp['SMILES'].replace('R', "*"))
                else:
                    outfile.write("%s\n" % comp['SMILES'])
                ids.append(comp['_id'])
        if platform.system() == 'Windows':
            rc = call("molconvert -mo \"%s/.%s\" %s \"%s\""
                      % (path, extension, img_type, structure_file),
                      shell=True)
        else:
            rc = call(["molconvert -mo %s/.%s %s %s"
                       % (path, extension, img_type, structure_file)],
                      shell=True)
        if rc:
            raise RuntimeError("molconvert returned %s" % rc)
        os.remove(structure_file)

        for i, _id in enumerate(ids):
            old = os.path.join(path, "%s.%s" % ((i + 1), extension))
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
        self.compounds.create_index([('Mass', pymongo.ASCENDING),
                                     ('Charge', pymongo.ASCENDING),
                                     ('DB_links.Model_SEED',
                                      pymongo.ASCENDING)])
        self.compounds.create_index([('Names', 'text'), ('Pathways', 'text')])
        self.compounds.create_index("DB_links.Model_SEED")
        self.compounds.create_index("DB_links.KEGG")
        self.compounds.create_index("MACCS")
        self.compounds.create_index("len_MACCS")
        self.compounds.create_index("RDKit")
        self.compounds.create_index("len_RDKit")
        self.compounds.create_index("Inchikey")
        self.compounds.create_index("MINE_id")
        self.compounds.create_index("Names")
        self.reactions.create_index("Reactants.c_id")
        self.reactions.create_index("Products.c_id")
        self.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                   "Action": "Database indexes built"})

    def link_to_external_database(self, external_database, compound=None,
                                  match_field="Inchikey", fields_to_copy=None):
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
            projection = dict([("_id", 0,)] + [(x[0], 1,) for x
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

    def insert_compound(self, mol_object, compound_dict=None, bulk=None,
                        kegg_db="KEGG", pubchem_db='PubChem-8-28-2015',
                        modelseed_db='ModelSEED'):
        """This class saves a RDKit Molecule as a compound entry in the MINE.
        Calculates necessary fields for API and includes additional
        information passed in the compound dict. Overwrites preexisting
        compounds in MINE on _id collision.

        :param mol_object: The compound to be stored
        :type mol_object: RDKit Mol object
        :param compound_dict: Additional information about the compound to be
            stored. Overwritten by calculated values.
        :type compound_dict: dict
        :param bulk: A pymongo bulk operation object. If None, reaction is
         immediately inserted in the database
        :param kegg_db: The ID of the KEGG Mongo database
        :type kegg_db: str
        :param pubchem_db: The ID of the PubChem Mongo database
        :type pubchem_db: str
        :param modelseed_db: The ID of the ModelSEED Mongo database
        :type modelseed_db: str
        :return: The hashed _id of the compound
        :rtype: str
        """

        if compound_dict is None:
            compound_dict = {}

        # Store all different representations of the molecule (SMILES, Formula,
        #  InChI key, etc.) as well as its properties in a dictionary
        compound_dict['SMILES'] = AllChem.MolToSmiles(mol_object, True)
        compound_dict['Inchi'] = AllChem.MolToInchi(mol_object)
        compound_dict['Inchikey'] = AllChem.InchiToInchiKey(
            compound_dict['Inchi'])
        compound_dict['Mass'] = AllChem.CalcExactMolWt(mol_object)
        compound_dict['Formula'] = AllChem.CalcMolFormula(mol_object)
        compound_dict['Charge'] = AllChem.GetFormalCharge(mol_object)
        # Get indices where bits are 1
        compound_dict['MACCS'] = list(AllChem.GetMACCSKeysFingerprint(
            mol_object).GetOnBits())
        compound_dict['len_MACCS'] = len(compound_dict['MACCS'])
        # Get indices where bits are 1
        compound_dict['RDKit'] = list(AllChem.RDKFingerprint(
            mol_object).GetOnBits())
        compound_dict['len_RDKit'] = len(compound_dict['RDKit'])
        compound_dict['logP'] = AllChem.CalcCrippenDescriptors(mol_object)[0]
        compound_dict['_id'] = utils.compound_hash(
            compound_dict['SMILES'],
            ('Type' in compound_dict
             and compound_dict['Type'] == 'Coreactant'))
        if '_atom_count' in compound_dict:
            del compound_dict['_atom_count']
        # Caching this for rapid reaction mass change calculation
        self._mass_cache[compound_dict['_id']] = compound_dict['Mass']

        # If the compound is a reactant, then make sure the reactant name is
        # in a correct format.
        if "Reactant_in" in compound_dict and isinstance(
                compound_dict['Reactant_in'], str) \
                and compound_dict['Reactant_in']:
            compound_dict['Reactant_in'] = ast.literal_eval(
                compound_dict['Reactant_in'])
        # If the compound is a product, then make sure the reactant name is
        # in a correct format.
        if "Product_of" in compound_dict \
                and isinstance(compound_dict['Product_of'], str) \
                and compound_dict['Product_of']:
            compound_dict['Product_of'] = ast.literal_eval(
                compound_dict['Product_of'])

        # Store links to external databases where compound is present
        if compound_dict['Inchikey']:
            if kegg_db:
                compound_dict = self.link_to_external_database(
                    kegg_db, compound=compound_dict,
                    fields_to_copy=[('Pathways', 'Pathways'),
                                    ('Names', 'Names'),
                                    ('DB_links', 'DB_links'),
                                    ('Enzymes', 'Enzymes')])

            if pubchem_db:
                compound_dict = self.link_to_external_database(
                    pubchem_db, compound=compound_dict,
                    fields_to_copy=[('COMPOUND_CID', 'DB_links.PubChem')])

            if modelseed_db:
                compound_dict = self.link_to_external_database(
                    modelseed_db, compound=compound_dict,
                    fields_to_copy=[('DB_links', 'DB_links')])

        # Calculate natural product likeness score and store in dict
        if not self.np_model:
            self.np_model = np.readNPModel()
        compound_dict["NP_likeness"] = np.scoreMol(mol_object, self.np_model)

        compound_dict = utils.convert_sets_to_lists(compound_dict)
        # Assign an id to the compound
        if self.id_db:
            mine_comp = self.id_db.compounds.find_one(
                {"Inchikey": compound_dict['Inchikey']},
                {'MINE_id': 1, "Pos_CFM_spectra": 1, "Neg_CFM_spectra": 1})
            # If compound already exists in MINE, store its MINE id in the dict
            if mine_comp:
                compound_dict['MINE_id'] = mine_comp['MINE_id']
                if 'Pos_CFM_spectra' in mine_comp:
                    compound_dict['Pos_CFM_spectra'] = \
                        mine_comp['Pos_CFM_spectra']
                if 'Neg_CFM_spectra' in mine_comp:
                    compound_dict['Neg_CFM_spectra'] = \
                        mine_comp['Neg_CFM_spectra']
            # If compound does not exist, create new id based on number of
            # current ids in the MINE
            else:
                compound_dict['MINE_id'] = self.id_db.compounds.count()
                self.id_db.compounds.save(compound_dict)

        # If bulk insertion, upsert (insert and update) the database
        if bulk:
            bulk.find({'_id': compound_dict['_id']}).upsert().\
                replace_one(compound_dict)
        elif self.compounds.find_one({'_id': compound_dict['_id']}):           
            self.compounds.replace_one({'_id': compound_dict['_id']},
                                       compound_dict)
        else:
            self.compounds.insert_one(compound_dict)
        return compound_dict['_id']

    def insert_reaction(self, reaction_dict, bulk=None):
        """Inserts a reaction into the MINE database and returns _id of the
         reaction in the mine database.

        :param reaction_dict: A dictionary containing "Reactants" and
         "Products" lists of StoichTuples
        :type reaction_dict: dict
        :param bulk: A pymongo bulk operation object. If None, reaction is
         immediately inserted in the database
        :return: The hashed _id of the reaction
        :rtype: str
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
        # If bulk insertion, upsert (insert and update) the database
        if bulk:
            bulk.find({'_id': reaction_dict['_id']}).upsert().\
                replace_one(reaction_dict)
        else:
            if '_id' in reaction_dict:
                self.reactions.replace_one({'_id': reaction_dict['_id']},
                                           reaction_dict)
            else:
                self.reactions.insert_one(reaction_dict)
        return reaction_dict['_id']

    def map_reactions(self, ext_db, match_field='_id'):
        """Update operators by adding the reactions they participate in.

        :param ext_db: name of MongoDB
        :param match_field: name of reaction field to match based on
        """
        # Update operators by adding the reactions they participate in
        lit_db = MINE(ext_db)
        for lit_rxn in lit_db.reactions.find():
            if match_field == "InChI_hash":
                mine_rxn = self.reactions.find_one(
                    {match_field: {"$regex": lit_rxn[match_field].split('-')[0]
                                   }}, {'Operators': 1})
            else:
                mine_rxn = self.reactions.find_one(
                    {match_field: lit_rxn[match_field]}, {'Operators': 1})
            if mine_rxn:
                for op in mine_rxn['Operators']:
                    self.operators.update(
                        {'_id': op},
                        {"$addToSet": {"Mapped_Rxns": lit_rxn["_id"],
                                       "References": {
                                           "$each": lit_rxn['References']}}})
                lit_db.reactions.update(
                    {"_id": lit_rxn["_id"]},
                    {"$set": {"Mapped_Rules": mine_rxn['Operators']}})
