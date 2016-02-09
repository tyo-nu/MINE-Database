#!/usr/bin/env python

"""Databases.py: This file contains MINE database classes including database loading functions."""

__author__ = 'JGJeffryes'

import pymongo
import platform
import hashlib
import utils
import datetime
from rdkit.Chem import AllChem

def establish_db_client():
    """This establishes a mongo database client in various environments"""
    try:
        #special case for working on a sapphire node
        if 'node' in platform.node():
            client = pymongo.MongoClient(host='master')
        #special case for working on a SEED cluster
        elif 'bio' in platform.node() or 'twig' == platform.node() or 'branch' == platform.node():
            client = pymongo.MongoClient(host='branch')
            admin = client['admin']
            admin.authenticate('worker', 'bnice14bot')
        #local database
        else:
            client = pymongo.MongoClient()
    except:
        raise IOError("Failed to load database client. Please verify that mongod is running")
    return client


class MINE:
    """
    This class basically exposes the underlying mongo database to manipulation but also defines expected database
    structure.
    """
    def __init__(self, name):
        self.client = establish_db_client()
        db = self.client[name]
        self._db = db
        self.name = name
        self.meta_data = db.meta_data
        self.compounds = db.compounds
        self.reactions = db.reactions
        self.operators = db.operators
        self.models = db.models
        self.id_db = self.client['UniversalMINE']

    def add_rxn_pointers(self):
        """Add links to the reactions that each compound participates in allowing for users to follow paths in the
         network"""
        reactions_count = self.reactions.count()
        print("Linking compounds to %s reactions" % reactions_count)
        for reaction in self.reactions.find().batch_size(500):
            for compound in reaction['Reactants']:
                self.compounds.update({"_id": compound[1]}, {'$push': {"Reactant_in": reaction['_id']}})
            for compound in reaction['Products']:
                self.compounds.update({"_id": compound[1]}, {'$push': {"Product_of": reaction['_id']}})
        self.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Add Reaction Pointers"})

    def fix_rxn_pointers(self, new_id, comp_dict):
        if self.reactions.count() and new_id != comp_dict['_id']:
            try:
                for reaction in comp_dict['Product_of']:
                    rxn = self.reactions.find_one({'_id': str(reaction)}, {'Products': 1})
                    for i, product in enumerate(rxn['Products']):
                        if product[1] == comp_dict['_id']:
                            rxn['Products'][i][1] = new_id
                    self.reactions.update({'_id': str(reaction)}, {'$set': {'Products': rxn['Products']}})
            except KeyError:
                pass

            try:
                for reaction in comp_dict['Reactant_in']:
                    rxn = self.reactions.find_one({'_id': str(reaction)}, {'Reactants': 1})
                    for i, reactant in enumerate(rxn['Reactants']):
                        if reactant[1] == comp_dict['_id']:
                            rxn['Reactants'][i][1] = new_id
                    self.reactions.update({'_id': str(reaction)}, {'$set': {'Reactants': rxn['Reactants']}})
            except KeyError:
                pass

        comp_dict['_id'] = new_id

        return comp_dict

    def add_compound_sources(self, rxn_key_type="_id"):
        for compound in self.compounds.find({"Sources": {"$exists": 0}}):
            compound['Sources'] = []
            for reaction in self.reactions.find({"Products.c_id": compound[rxn_key_type]}):
                compound['Sources'].append({"Compounds": [x['c_id'] for x in reaction['Reactants']], "Operators": reaction["Operators"]})
            if compound['Sources']:
                try:
                    self.compounds.save(compound)
                except pymongo.errors.DocumentTooLarge:
                    print("Too Many Sources for %s" % compound['SMILES'])
        self.compounds.ensure_index([("Sources.Compound", pymongo.ASCENDING), ("Sources.Operators", pymongo.ASCENDING)])
        self.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Add Compound Source field"})

    def link_to_external_database(self, external_database, match_field="Inchikey", fields_to_copy=None):
        """
        This function looks for matching compounds in other databases (i.e. PubChem) and adds links where found

        :param external_database: String, the name of the database to search for matching compounds
        :param match_field: String, The field to search on for matching compunds
        :param fields_to_copy: List of tuples, data to copy into the mine database. The first field is the field name in
        the external database. The second field is the field name in the MINE database where the data will be copied.
        :return:
        """
        pass

    def insert_compound(self, mol_object, compound_dict={}, kegg_db="KEGG", pubchem_db='PubChem-8-28-2015',
                        modelseed_db='ModelSEED'):
        """
        This class saves a RDKit Molecule as a compound entry in the MINE. Calculates necessary fields for API and
        includes additional information passed in the compound dict. Overwrites preexisting compounds in MINE on _id
        collision.
        :param mol_object: The compound to be stored
        :type mol_object: RDKit Mol object
        :param compound_dict: Additional information about the compound to be stored. Overwritten by calculated values.
        :type compound_dict: dict
        :return:
        :rtype:
        """

        compound_dict['SMILES'] = AllChem.MolToSmiles(mol_object, True)
        compound_dict['Inchi'] = AllChem.MolToInchi(mol_object)
        compound_dict['Inchikey'] = AllChem.InchiToInchiKey(compound_dict['Inchi'])
        compound_dict['Mass'] = AllChem.CalcExactMolWt(mol_object)
        compound_dict['Formula'] = AllChem.CalcMolFormula(mol_object)
        compound_dict['Charge'] = AllChem.GetFormalCharge(mol_object)
        compound_dict['MACCS'] = [i for i, bit in enumerate(AllChem.GetMACCSKeysFingerprint(mol_object)) if bit]
        compound_dict['len_MACCS'] = len(compound_dict['MACCS'])
        compound_dict['RDKit'] = [i for i, bit in enumerate(AllChem.RDKFingerprint(mol_object)) if bit]
        compound_dict['len_RDKit'] = len(compound_dict['RDKit'])
        comphash = hashlib.sha1(compound_dict['SMILES'].encode('utf-8')).hexdigest()
        if '_id' in compound_dict:
            if "X" in compound_dict['_id']:
                    compound_dict = self.fix_rxn_pointers('X' + comphash, compound_dict)
            else:
                compound_dict = self.fix_rxn_pointers('C' + comphash, compound_dict)
        else:
            compound_dict['_id'] = 'C' + comphash

        compound_dict['DB_links'] = {}
        if compound_dict['Inchikey']:
            if kegg_db:
                kegg_comps = self.client[kegg_db].compounds.find({"Inchikey": compound_dict['Inchikey']},
                                                    {'Pathways': 1, 'Names': 1, 'DB_links': 1, 'Enzymes': 1, '_id': 0})
                for kcomp in kegg_comps:
                    utils.dict_merge(compound_dict, kcomp)

            if pubchem_db:
                pubchem_comps = self.client[pubchem_db].compounds.find({"Inchikey": compound_dict['Inchikey']},
                                                                       {'COMPOUND_CID': 1, "IUPAC_NAME": 1})
                if pubchem_comps.count() and "PubChem" not in compound_dict['DB_links']:
                    compound_dict['DB_links']['PubChem'] = [x['COMPOUND_CID'] for x in pubchem_comps]

            if modelseed_db:
                for seed_comp in self.client[modelseed_db].compounds.find({"Inchikey": compound_dict['Inchikey']},
                                                                          {'DB_links.Model_SEED': 1}):
                    if 'Model_SEED' not in compound_dict['DB_links']:
                        compound_dict['DB_links']['Model_SEED'] = []
                    compound_dict['DB_links']['Model_SEED'].extend(seed_comp['DB_links']['Model_SEED'])

        if self.id_db:
            mine_comp = self.id_db.compounds.find_one({"Inchikey": compound_dict['Inchikey']})
            if mine_comp:
                compound_dict['MINE_id'] = mine_comp['MINE_id']
            else:
                compound_dict['MINE_id'] = self.id_db.compounds.count()
                self.id_db.compounds.save(utils.convert_sets_to_lists(compound_dict))
        self.compounds.save(utils.convert_sets_to_lists(compound_dict))

    def insert_reaction(self, reaction_dict):
        pass
