__author__ = 'JGJeffryes'

import pylab as P
from minedatabase import databases

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
P.show()