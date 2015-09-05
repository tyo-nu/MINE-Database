#!/usr/bin/env python

"""Databases.py: This file contains MINE database classes including database loading functions."""

__author__ = 'JGJeffryes'

import pymongo
import platform

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
        client = establish_db_client()
        db = client[name]
        self.db = db
        self.name = name
        self.meta_data = self.meta_data
        self.compounds = db.compounds
        self.reactions = db.reactions
        self.models = db.models

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
