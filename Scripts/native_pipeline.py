"""For testing getting comps that are product of model compounds"""


import pymongo

uri = 'mongodb://username:password@minedatabase.ci.northwestern.edu:27017'
db_name = 'kegg_lte600_500mcy'


cpd_ids = ['C000000588c5adcceec7bc0033b2404db527f006f',  # just 1 rxn
           'C000000bceed4d371ac1164f77db9170e0944693e',  # just 1 rxn
           'C3f8881998951e22f122db498b619af1eea58d2c1',  # has 2 product_of vals bc has 1,000s of rxns
           'C32128af794b3adefa47a88712cc1967af41098de']  # last one has 3 rxns

model_ids = ['C12051', 'C21976', 'C99999', 'C00024', 'C00125', 'C13144', 'C14574']


client = pymongo.MongoClient(uri, ServerSelectionTimeoutMS=5)

db = client[db_name]

pipeline = [
    {
        '$match': {
            '_id': {'$in': cpd_ids}
        }
    },
    {
        '$project': {
            'Product_of': 1,
        }
    },
    {
        '$unwind': '$Product_of'
    },
    {
        '$lookup': {
            'from': 'product_of',
            'localField': 'Product_of',
            'foreignField': '_id',
            'as': 'product_of_doc'
        }
    },
    {
        '$unwind': '$product_of_doc'
    },
    {
        '$project': {
            '_id': 1,
            'producing_rxns': '$product_of_doc.Product_of',
        }
    },
    {
        '$unwind': '$producing_rxns'
    },
    {
        '$group': {
            '_id': '$_id',
            'producing_rxns': {'$addToSet': '$producing_rxns'},
            'n_rxns': {'$sum': 1}
        }
    },
    {
        '$unwind': '$producing_rxns'
    },
    {
        '$lookup': {
            'from': 'reactions',
            'localField': 'producing_rxns',
            'foreignField': '_id',
            'as': 'reaction_doc'
        }
    },
    {
        '$project': {
            'reactants': '$reaction_doc.Reactants',
        }
    },
    {
        '$project': {
            'reaction_doc': 0,
            'producing_rxns': 0
        }
    },
    {
        '$unwind': '$reactants'
    },
    {
        '$unwind': '$reactants'
    },
    {
        '$project': {
            'reactant': {'$arrayElemAt': ['$reactants', 1]}
        }
    },
    {
        '$match': {
            'reactant': {'$regex': '^C.*'}
        }
    },
    {
        '$lookup': {
            'from': 'compounds',
            'localField': 'reactant',
            'foreignField': '_id',
            'as': 'reactant_doc'
        }
    },
    {
        '$project': {
            '_id': 1,
            'KEGG_id': {'$arrayElemAt': ['$reactant_doc.ID', 0]}
        }
    },
    {
        '$group': {
            '_id': '$_id',
            'KEGG_ids': {'$addToSet': '$KEGG_id'}
        }
    },
    {
        '$match': {
            'KEGG_ids': {'$elemMatch': {'$in': model_ids}}
        }
    },
    {
        '$project': {
            '_id': 1,
        }
    }
]

compounds = db.compounds.aggregate(pipeline)

compounds = list(compounds)

print(compounds)
