__author__ = 'JGJeffryes'

from rdkit.Chem import AllChem
from pymongo import ASCENDING
from databases import MINE
from argparse import ArgumentParser
import itertools
import collections
import time
import utils
from copy import deepcopy
import datetime
import hashlib
import csv
import multiprocessing

stoich_tuple = collections.namedtuple("stoich_tuple", 'compound,stoich')

class Pickaxe:
    """
        This class generates new compounds from user-specified starting compounds using a set of SMARTS-based reaction
        rules. It may be initialized with a text file containing the reaction rules and cofactors or this may be
        done on an ad hock basis.
    """
    def __init__(self, rule_list=None, cofactor_list=None, explicit_h=True, kekulize=True, errors=True,
                 raceimze=False, split_stereoisomers=True, mine=None):
        self.rxn_rules = {}
        self.cofactors = {}
        self._raw_compounds = {}
        self.compounds = collections.OrderedDict()
        self.reactions = collections.OrderedDict()
        self.mine = mine
        self.generation = -1
        self.explicit_h = explicit_h
        self.split_stereoisomers = split_stereoisomers
        self.kekulize = kekulize
        self.raceimize = raceimze
        self.stoich_tuple = stoich_tuple
        self.errors = errors
        from rdkit import RDLogger
        lg = RDLogger.logger()
        if not errors:
            lg.setLevel(4)
        if cofactor_list:
            with open(cofactor_list) as infile:
                for cofactor in infile:
                    self._load_cofactor(cofactor)
        self.cofactors['set'] = set()
        if rule_list:
            with open(rule_list) as infile:
                for rule in infile:
                    if rule[0] == "#":
                        continue
                    self.load_rxn_rule(rule)

    def _load_cofactor(self, cofactor_text):
        """
        Loads a cofactor into the cofactor dictionary from a tab-delimited string
        :param cofactor_text: tab-delimited string with the compound name and SMILES
        :type cofactor_text: basestring
        :return:
        :rtype:
        """
        split_text = cofactor_text.strip().split('\t')
        try:
            mol = AllChem.MolFromSmiles(split_text[1])
            if not mol:
                raise ValueError
            smi = AllChem.MolToSmiles(mol, True)
        except (IndexError, ValueError):
            raise ValueError("Unable to load cofactor: %s" % cofactor_text)
        self._add_compound(split_text[0], smi, mol=mol)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        self.cofactors[split_text[0]] = mol
        self._raw_compounds[AllChem.MolToSmiles(mol, True)] = smi

    def load_rxn_rule(self, rule_text):
        """
        Loads a reaction rule into the rxn_rule dictionary from a tab-delimited string
        :param rule_text: tab-delimited string with the rule name, cofactor names, and rule as SMARTS
        :type rule_text: basestring
        :return:
        :rtype:
        """
        if rule_text[0] == "#":
            return
        split_text = rule_text.strip().split('\t')
        if len(split_text) < 3:
            raise ValueError("Unable to parse reaction rule: %s" % rule_text)
        reactant_names = split_text[1].split(';')
        for name in reactant_names:
            if name == "Any":
                continue
            if name not in self.cofactors:  # try to proceed as if name is SMILES
                self._load_cofactor(name+"\t"+name)
        rxn = AllChem.ReactionFromSmarts(split_text[2])
        if rxn.GetNumReactantTemplates() != len(reactant_names):
            raise ValueError("Number of cofactors does not match supplied reaction rule: %s" % rule_text)
        self.rxn_rules[split_text[0]] = (reactant_names, rxn)

    def load_compound_set(self, compound_file=None, structure_field='structure', id_field='id'):
        """
        If a compound file is provided, this function loads the compounds into it's internal dictionary. If not, it
        attempts to find the compounds in it's associated MINE database.
        :param compound_file: Path to a file containing compounds as tsv
        :type compound_file: basestring
        :param structure_field: the name of the column containing the structure incarnation as Inchi or SMILES (Default:
        'structure')
        :type structure_field: str
        :param id_field: the name of the column containing the desired compound ID (Default: 'id)
        :type id_field: str
        :return: compound SMILES
        :rtype: list
        """
        compound_smiles = []
        self.generation = 0
        if compound_file:
            with open(compound_file) as infile:
                reader = csv.DictReader(infile.readlines(), delimiter="\t")
                for line in reader:
                    if "InChI=" in line[structure_field]:
                        mol = AllChem.MolFromInchi(line[structure_field])
                    else:
                        mol = AllChem.MolFromSmiles(line[structure_field])
                    if not mol:
                        if self.errors:
                            print("Unable to Parse %s" % line[structure_field])
                        continue
                    smi = AllChem.MolToSmiles(mol, True)
                    id = line[id_field]
                    self._add_compound(id, smi, mol=mol)
                    compound_smiles.append(smi)
        else:
            if not self.mine:
                raise ValueError('MINE database not specified')
            db = MINE(self.mine)
            for compound in db.compounds.find():
                id = compound['_id']
                smi = compound['SMILES']
                self._add_compound(id, smi)
                compound_smiles.append(smi)
        print("%s compounds loaded" % len(compound_smiles))
        return compound_smiles

    def _add_compound(self, id, smi, mol=None):
        if not mol:
            mol = AllChem.MolFromSmiles(smi)
        try:
            i_key = AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))
            self.compounds[smi] = {'ID': id, '_id': smi, "SMILES": smi, 'Inchikey': i_key, 'Generation': self.generation}
        except ValueError:
            self.compounds[smi] = {'ID': id, '_id': smi, "SMILES": smi, 'Inchikey': "None", 'Generation': self.generation}
        self._raw_compounds[smi] = smi

    def transform_compound(self, compound_SMILES, rules=None):
        """
        Perform transformations to a compound returning the products and the predicted reactions
        :param compound_SMILES: The compound on which to operate represented as SMILES
        :type compound_SMILES: string
        :param rules: The names of the reaction rules to apply. If none, all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as tuple of id abd SMILES string and reactions as a tuple of reactants & products
        :rtype: tuple of lists
        """
        pred_rxns = set()
        pred_compounds = set()
        if not rules:
            rules = self.rxn_rules.keys()
        mol = AllChem.MolFromSmiles(compound_SMILES)
        if not mol:
            if self.errors:
                raise ValueError('Unable to parse: %s' % compound_SMILES)
            else:
                print('Unable to parse: %s' % compound_SMILES)
                return
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        self.cofactors['Any'] = mol
        for rule_name in rules:
            rule = self.rxn_rules[rule_name]
            reactant_mols = tuple([self.cofactors[x] for x in rule[0]])
            try:
                product_sets = rule[1].RunReactants(reactant_mols)
            except RuntimeError:  # I need to do more to untangle the causes of this error
                print("Runtime ERROR!"+rule_name)
                print(compound_SMILES)
                continue
            reactants = next(self._make_compound_tups(reactant_mols, rule_name))  # no enumeration for reactants
            for product_mols in product_sets:
                try:
                    for stereo_prods in self._make_compound_tups(product_mols, rule_name, split_stereoisomers=self.split_stereoisomers):
                        pred_compounds.update(x.compound for x in stereo_prods)
                        rhash = self._calculate_rxn_hash(reactants, stereo_prods)
                        text_rxn = ' + '.join(['(%s) %s' % (x.stoich, x.compound) for x in reactants]) + ' => ' + \
                           ' + '.join(['(%s) %s' % (x.stoich, x.compound) for x in stereo_prods])
                        pred_rxns.add(text_rxn)
                        if rhash not in self.reactions:
                            reaction_data = {"_id": rhash, "Reactants": reactants, "Products": stereo_prods,
                                             "Operators": {rule_name}, "SMILES_rxn": text_rxn}
                            self.reactions[rhash] = reaction_data
                        else:
                            self.reactions[rhash]['Operators'].add(rule_name)
                except ValueError:
                    continue
        return pred_compounds, pred_rxns

    def _make_compound_tups(self, mols, rule_name, split_stereoisomers=False):
        """Takes a list of mol objects and returns an generator for (compound, stoich) tuples"""
        comps = []
        products = []
        for m in mols:
            r_smiles = AllChem.MolToSmiles(m, True)
            try:
                comps.append(self._calculate_compound_information(r_smiles, m))
            except ValueError:
                print("Product ERROR!: %s %s" % (rule_name, r_smiles))
                raise ValueError
        if split_stereoisomers:
            for subrxn in itertools.product(*comps):
                products.append(collections.Counter(subrxn))
        else:
            products = [collections.Counter(comps)]
        for rxn in products:
            yield [self.stoich_tuple(x, y) if len(x) > 1 else self.stoich_tuple(x[0], y) for x, y in rxn.items()]

    def _calculate_compound_information(self, raw, mol_obj):
        """Calculate the standard data for a compound & return a tuple with compound_ids. Memoized with _raw_compound
        dict"""
        if raw not in self._raw_compounds:
            if self.explicit_h:
                mol_obj = AllChem.RemoveHs(mol_obj)  # this step slows down the process quite a bit
            AllChem.SanitizeMol(mol_obj)
            if self.raceimize:
                mols = self._racemization(mol_obj)
            else:
                mols = [mol_obj]
            smiles = [AllChem.MolToSmiles(m, True) for m in mols]
            for s, m in zip(smiles, mols):
                if s not in self._raw_compounds:
                    self._add_compound(None, s, mol=m)
            self._raw_compounds[raw] = tuple([self._raw_compounds[s] for s in smiles])
        return self._raw_compounds[raw] if isinstance(self._raw_compounds[raw], tuple) else (self._raw_compounds[raw],)

    def _racemization(self, compound, max_centers=3, carbon_only=True):
        """
        Enumerates all possible stereoisomers for unassigned chiral centers.
        :param compound: A compound
        :type compound: rdMol object
        :param max_centers: The maximum number of unspecified stereocenters to enumerate. Sterioisomers grow 2^n_centers
        :type max_centers: int
        :param carbon_only: Only enumerate unspecified carbon centers. (other centers are often not tautomeric artifacts)
        :type carbon_only:  bool
        :return: list of stereoisomers
        :rtype: list of rdMol objects
        """
        new_comps = []
        unassigned_centers = [c[0] for c in AllChem.FindMolChiralCenters(compound, includeUnassigned=True) if c[1] == "?"]
        if carbon_only:
            unassigned_centers = list(filter(lambda x: compound.GetAtomWithIdx(x).GetAtomicNum() == 6, unassigned_centers))
        if not unassigned_centers or len(unassigned_centers) > max_centers:
            return [compound]
        for seq in itertools.product([1, 0], repeat=len(unassigned_centers)):
            for atomId, cw in zip(unassigned_centers, seq):
                if cw:
                    compound.GetAtomWithIdx(atomId).SetChiralTag(AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
                else:
                    compound.GetAtomWithIdx(atomId).SetChiralTag(AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            new_comps.append(deepcopy(compound))
        return new_comps

    def _calculate_rxn_hash(self, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is connectivity only, second block is sterio
        only"""
        def __get_blocks(tups):
            first_block, second_block = [], []
            for x in tups:
                if self.compounds[x.compound]["Inchikey"]:
                    split_inchikey = self.compounds[x.compound]["Inchikey"].split('-')
                    if len(split_inchikey) > 1:
                        first_block.append("%s,%s" % (x.stoich, split_inchikey[0]))
                        second_block.append("%s,%s" % (x.stoich, split_inchikey[1]))
                else:
                    print("No Inchikey for %s" % x.compound)
            return "+".join(first_block), "+".join(first_block)

        reactants.sort()
        products.sort()
        r_1, r_2 = __get_blocks(reactants)
        p_1, p_2 = __get_blocks(products)
        first_block = r_1+'<==>'+p_1
        second_block = r_2+'<==>'+p_2

        return hashlib.sha256(first_block.encode()).hexdigest()+"-"+hashlib.md5(second_block.encode()).hexdigest()

    def _assign_ids(self):
        self.compounds = dict(self.compounds)
        self.reactions = dict(self.reactions)
        i = 1
        for comp in self.compounds.values():
            if not comp['ID']:
                comp['ID'] = '_cpd'+str(i).zfill(7)
                i += 1
                self.compounds[comp['_id']] = comp
        for rxn in self.reactions.values():
            rxn['ID_rxn'] = ' + '.join(['(%s) %s[c0]' % (x.stoich, self.compounds[x.compound]["ID"]) for x in rxn["Reactants"]]) \
                            + ' => ' + ' + '.join(['(%s) %s[c0]' % (x.stoich, self.compounds[x.compound]["ID"]) for x in rxn["Products"]])
            self.reactions[rxn['_id']] = rxn

    def transform_all(self, num_workers=1, max_generations=1):
        while self.generation < max_generations:
            self.generation += 1
            ti = time.time()
            n_comps = len(self.compounds)
            n_rxns = len(self.reactions)
            compound_smiles = [c['SMILES'] for c in self.compounds.values() if c['Generation'] == self.generation - 1]
            print_on = max(round(.05 * len(compound_smiles)), 1)
            if compound_smiles:
                if num_workers > 1:
                    pool = multiprocessing.Pool(processes=num_workers)
                    manager = multiprocessing.Manager()
                    self.compounds = manager.dict(self.compounds)
                    self.reactions = manager.dict(self.reactions)
                    for i, res in enumerate(pool.imap_unordered(self.transform_compound, compound_smiles)):
                        if not (i+1) % print_on:
                            print("Generation %s: %s percent complete" % (self.generation, round((i+1) * 100 / len(compound_smiles))))
                else:
                    for i, smi in enumerate(compound_smiles):
                        self.transform_compound(smi)
                        if not (i+1) % print_on:
                            print("Generation %s: %s percent complete" % (self.generation, round((i+1) * 100 / len(compound_smiles))))
            print("Generation %s produced %s new compounds and %s new reactions in %s sec" % (self.generation,
                len(self.compounds)-n_comps, len(self.reactions) - n_rxns, time.time()-ti))
        self._assign_ids()

    def write_compound_output_file(self, path, delimiter='\t'):
        """
        Writes all compound data to the specified path.
        :param path: path to output
        :type path: basestring
        :param delimiter: the character with which to separate data entries
        :type delimiter: basestring
        :return:
        :rtype:
        """
        path = utils.prevent_overwrite(path)
        with open(path, 'w') as outfile:
            outfile.write('ID\tInChIKey\tSMILES\n')
            for c in sorted(self.compounds.values(), key=lambda x: x['ID']):
                outfile.write(delimiter.join([c['ID'], str(c['Inchikey']), c['SMILES']])+'\n')

    def write_reaction_output_file(self, path, delimiter='\t'):
        """
        Writes all reaction data to the specified path.
        :param path: path to output
        :type path: basestring
        :param delimiter: the character with which to separate data entries
        :type delimiter: basestring
        :return:
        :rtype:
        """
        path = utils.prevent_overwrite(path)
        with open(path, 'w') as outfile:
            outfile.write('ID\tName\tID Equation\tSMILES Equation\tRxn Hash\tOperators\n')
            for i, rxn in enumerate(self.reactions.values()):
                outfile.write(delimiter.join(['_rxn'+str(i).zfill(7), '', rxn['ID_rxn'], rxn["SMILES_rxn"], rxn['_id'],
                                              ';'.join(rxn['Operators'])])+'\n')

    def save_to_MINE(self, db_id):
        """
        Save compounds to a MINE database.
        :param db_id: The name of the target database
        :type db_id: basestring
        :return:
        :rtype:
        """
        db = MINE(db_id)
        for comp_dict in self.compounds.values():
            db.insert_compound(AllChem.MolFromSmiles(comp_dict['SMILES']), comp_dict)
        db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Compounds Inserted"})
        for rxn in self.reactions.values():
            rxn['Reactants'] = [tup._asdict() for tup in rxn['Reactants']]
            rxn['Products'] = [tup._asdict() for tup in rxn['Products']]
            rxn = utils.convert_sets_to_lists(rxn)
            db.reactions.save(rxn)
        db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Reactions Inserted"})

        db.compounds.ensure_index([('Mass', ASCENDING), ('Charge', ASCENDING), ('DB_links.Model_SEED', ASCENDING)])
        db.compounds.ensure_index([('Names', 'text'), ('Enzymes', 'text'), ('Pathways', 'text')])
        db.compounds.ensure_index("DB_links.Model_SEED")
        db.compounds.ensure_index("DB_links.KEGG")
        db.compounds.ensure_index("MACCS")
        db.compounds.ensure_index("len_MACCS")
        db.compounds.ensure_index("Inchikey")
        db.compounds.ensure_index("MINE_id")

        db.reactions.ensure_index("Reactants.compound")
        db.reactions.ensure_index("Products.compound")


if __name__ == "__main__":
    t1 = time.time()
    parser = ArgumentParser()
    parser.add_argument('-C', '--cofactor_list', default="Tests/Cofactor_SMILES.tsv", help="Specify a list of cofactors"
                                                                                           " as a tab-separated file")
    parser.add_argument('-r', '--rule_list', default="Tests/operators_smarts.tsv", help="Specify a list of reaction "
                                                                                        "rules as a tab-separated file")
    parser.add_argument('-c', '--compound_file', default="Tests/test_compounds.tsv", help="Specify a list of starting "
                        "compounds as a tab-separated file")
    parser.add_argument('-s', '--smiles', default=None, help="Specify a starting compound as SMILES.")
    parser.add_argument('-o', '--output_dir', default=".", help="The directory in which to place files")
    parser.add_argument('-d', '--database', default=None, help="The name of the database in which to store output. If "
                                                               "not specified, data is still written as tsv files")
    parser.add_argument('-R', '--raceimize', action='store_true', default=False, help="Enumerate the possible chiral "
                        "forms for all unassigned steriocenters in compounds & reactions")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="Display RDKit errors & warnings")
    parser.add_argument('--bnice', action='store_true', default=False, help="Set several options to enable "
                                                                            "compatability with bnice operators.")
    parser.add_argument('-m', '--max_workers', default=1, type=int, help="Set the nax number of processes to spawn to "
                                                                         "perform calculations.")
    parser.add_argument('-g', '--generations', default=1, type=int, help="Set the numbers of time to apply the reaction"
                                                                         " rules to the compound set.")
    options = parser.parse_args()
    pk = Pickaxe(cofactor_list=options.cofactor_list, rule_list=options.rule_list, raceimze=options.raceimize,
                 errors=options.verbose, explicit_h=options.bnice, kekulize=options.bnice)
    if options.smiles:
        pk._add_compound("", "Start", options.smiles)
    else:
        pk.load_compound_set(compound_file=options.compound_file)
    pk.transform_all(max_generations=options.generations, num_workers=options.max_workers)
    pk.write_compound_output_file(options.output_dir+'/compounds.tsv')
    pk.write_reaction_output_file(options.output_dir+'/reactions.tsv')
    if options.database:
        print("Saving results to %s" % options.database)
        pk.save_to_MINE(options.database)
    print("Execution took %s seconds." % (time.time()-t1))
