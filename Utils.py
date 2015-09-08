__author__ = 'JGJeffryes'

def convert_sets_to_lists(obj):
    """Recursively converts dictionaries that contain sets to lists"""
    if isinstance(obj, set):
        try:
            obj = sorted(list(obj), key=lambda x: len(x))  # this brings short names to the top of the list
        except TypeError:
            obj = list(obj)
    elif isinstance(obj, dict):
        for key in obj:
            obj[key] = convert_sets_to_lists(obj[key])
    return obj


def fix_rxn_pointers(db, new_id, comp_dict):
    if new_id != comp_dict['_id']:
        try:
            for reaction in comp_dict['Product_of']:
                rxn = db.reactions.find_one({'_id': str(reaction)}, {'Products': 1})
                for i, product in enumerate(rxn['Products']):
                    if product[1] == comp_dict['_id']:
                        rxn['Products'][i][1] = new_id
                db.reactions.update({'_id': str(reaction)}, {'$set': {'Products': rxn['Products']}})
        except KeyError:
            pass

        try:
            for reaction in comp_dict['Reactant_in']:
                rxn = db.reactions.find_one({'_id': str(reaction)}, {'Reactants': 1})
                for i, reactant in enumerate(rxn['Reactants']):
                    if reactant[1] == comp_dict['_id']:
                        rxn['Reactants'][i][1] = new_id
                db.reactions.update({'_id': str(reaction)}, {'$set': {'Reactants': rxn['Reactants']}})
        except KeyError:
            pass

        comp_dict['_id'] = new_id

    return comp_dict
