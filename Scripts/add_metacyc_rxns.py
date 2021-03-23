import sys
import re
import collections

import csv
from minedatabase import utils
from minedatabase.databases import MINE
from rdkit.Chem import AllChem


def make_hash_dict(db, key_field):
    hash_dict = {}
    for comp in db.compounds.find({key_field: {'$exists': 1}}, {key_field: 1}):
        for name in utils.get_dotted_field(comp, key_field):
            hash_dict[name] = comp['_id']
    return hash_dict

def dict_from_sdf(path):
    suppl = AllChem.SDMolSupplier(path)
    comp_dict = {}
    for mol in suppl:
        if mol:
            comp_dict[mol.GetProp('FRAME-ID')] = mol
    return comp_dict


def add_metacyc_rxns(mine_db, csv_path, metacyc2hash):
    inserted = set()
    def parse_comps(field):
        atoms = collections.Counter()
        compounds = collections.Counter(field.split(' // '))
        half_rxn = []
        for comp, stoich in compounds.items():
            if comp in metacyc2hash:
                mol = metacyc2hash[comp]
                for pair in re.findall('([A-Z][a-z]*)(\d*)', AllChem.CalcMolFormula(mol)):
                    if pair[1]:
                        atoms[pair[0]] += int(pair[1]) * stoich
                    else:
                        atoms[pair[0]] += 1 * stoich
                if comp not in inserted:
                    mine_db.insert_compound(mol, {'Generation': 0})
                    inserted.add(comp)
                half_rxn.append(utils.stoich_tuple(stoich, utils.get_compound_hash(mol)))
            else:
                raise ValueError('Undefined Compound: %s' % comp)
        return half_rxn, atoms

    with open(csv_path) as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            rxn = {'Metabolite': '', 'Type': '', 'MetaCyc ID': row['MetaCyc ID']}
            if isinstance(row['Citations'], str):
                rxn['References'] = [x.strip('"[]"') for x in row['Citations'].split(' // ')]
            else:
                rxn['References'] = [str(row['Citations'])]
            rxn['References'].append("MetaCyc: %s" % row['MetaCyc ID'])
            try:
                rxn['Reactants'], r_atoms = parse_comps(row['Reactants of reaction'])
                rxn['Products'], p_atoms = parse_comps(row['Products of reaction'])
                if r_atoms - p_atoms or p_atoms - r_atoms:
                    print(r_atoms, p_atoms)
                    raise ValueError('Unbalanced Reaction: %s' % rxn['MetaCyc ID'])
                if sorted(rxn['Reactants']) == sorted(rxn['Products']):
                    raise ValueError('No Change: %s' % rxn['MetaCyc ID'])

            except ValueError as e:
                print(e)
                continue
            mine_db.insert_reaction(rxn)
            """reactions = pd.read_csv(csv_path, sep='\t', error_bad_lines=False).fillna("")
                for i, row in reactions.iterrows():
                    rxn = row[['MetaCyc ID']].to_dict()
                    rxn['Metabolite'], rxn['Type'] = "", "" """


def add_metacyc_comps(metacyc_db, mine_db):
    c_ids = set(mine_db.reactions.distinct("Reactants.c_id"))
    c_ids |= set(mine_db.reactions.distinct("Products.c_id"))
    for _id in c_ids:
        if not mine_db.compounds.count({"_id": _id}):
            comp = metacyc_db.compounds.find_one({"_id": _id})
            mine_db.compounds.insert(comp)


if __name__ == '__main__':
    AllChem.WrapLogs()
    db = MINE(sys.argv[1])
    hash_dict = dict_from_sdf(sys.argv[2])
    add_metacyc_rxns(db, sys.argv[3], hash_dict)
    add_metacyc_comps(db, MINE(sys.argv[1]))