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
