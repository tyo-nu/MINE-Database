__author__ = 'JGJeffryes'

import pylab as P
from minedatabase import databases

db = databases.MINE('cdmine_iaf')
seed_comps = set()


def is_seed(rxn):
    for tup in rxn['Reactants']+rxn['Products']:
        _id = tup['c_id']
        if _id not in seed_comps:
            if not db.compounds.count(
                    {'_id': _id, 'DB_links.Model_SEED': {'$exists': 1}}):
                return 0
            seed_comps.add(_id)
    return 1
for rxn in db.reactions.find():
    rxn['AllSeed'] = is_seed(rxn)
    db.reactions.save(rxn)
print(db.reactions.count({'AllSeed': 1}))

"""
db_list = ['YMDBexp2', 'MetaCyc']
results = []
for name in db_list:
    db = databases.MINE(name)

    res = db.compounds.aggregate([{"$match": {"Generation": 0}}, {'$group': {'_id': {'$substr': ['$Inchikey', 0, 14]}, 'Isomers': {'$sum': 1}}},
                              {'$group': {"_id": '$Isomers', "Groups": {'$sum': 1}}}, {"$sort": {"_id": 1}}
                              ])
    results.append([x["Isomers"] for x in res])
P.figure()  # dpi=600, figsize=(3.35, 3.35))
#P.hist(results, range=(-3, 3), bins=200, normed=True, histtype='stepfilled', alpha= 0.5, label=["PubChem", "KEGG", 'KEGG MINE'])
P.hist(results, range=(1, 10), bins=10, normed=True, histtype='stepfilled', alpha=0.5, label=db_list)
P.legend(prop={'size': 14})
P.rcParams.update({'font.size': 14})
P.ylabel("Density")
P.xlabel("Natural Product Likeness Score")
P.show()"""