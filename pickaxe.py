__author__ = 'JGJeffryes'

from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToFile
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
import os

stoich_tuple = collections.namedtuple("stoich_tuple", 'compound,stoich')

class Pickaxe:
    def __init__(self, rule_list=None, cofactor_list=None, explicit_h=True, kekulize=True, errors=True,
                 raceimze=False, split_stereoisomers=True, database=None, image_dir=None):
        """
        This class generates new compounds from user-specified starting compounds using a set of SMARTS-based reaction
        rules. It may be initialized with a text file containing the reaction rules and cofactors or this may be
        done on an ad hock basis.

        :param rule_list: Path to a list of reaction rules in TSV form
        :type rule_list: str
        :param cofactor_list: Path to list of cofactors in TSV form
        :type cofactor_list: str
        :param explicit_h: Explicitly represent bound hydrogen atoms
        :type explicit_h: bool
        :param kekulize: Kekulize structures before applying reaction rules
        :type kekulize: bool
        :param errors: Print underlying RDKit warnings and halt on error
        :type errors: bool
        :param raceimze: Enumerate all possible chiral forms of a molecule if unspecified stereocenters exist
        :type raceimze: bool
        :param split_stereoisomers: Split each possible stereoisomers into a separate reaction
        :type split_stereoisomers: bool
        :param database: Name of desired Mongo Database
        :type database: str
        :param image_dir: Path to desired image folder
        :type image_dir: str
        """
        self.rxn_rules = {}
        self.cofactors = {}
        self._raw_compounds = {}
        self.compounds = {}
        self.reactions = {}
        self.mine = database
        self.generation = -1
        self.explicit_h = explicit_h
        self.split_stereoisomers = split_stereoisomers
        self.kekulize = kekulize
        self.raceimize = raceimze
        self.image_dir = image_dir
        self.errors = errors
        # TODO: Test database and warn on overwrite

        from rdkit import RDLogger
        lg = RDLogger.logger()
        if not errors:
            lg.setLevel(4)

        if cofactor_list:
            with open(cofactor_list) as infile:
                for cofactor in infile:
                    self._load_cofactor(cofactor)

        if rule_list:
            with open(rule_list) as infile:
                for rule in infile:
                    self.load_rxn_rule(rule)

    def _load_cofactor(self, cofactor_text):
        """
        Loads a cofactor into the cofactor dictionary from a tab-delimited string
        :param cofactor_text: tab-delimited string with the compound name and SMILES
        :type cofactor_text: basestring
        :return:
        :rtype:
        """
        #TODO: replace with csv function
        if cofactor_text[0] == "#":
            return
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

    def load_rxn_rule(self, rule_text):
        """
        Loads a reaction rule into the rxn_rule dictionary from a tab-delimited string
        :param rule_text: tab-delimited string with the rule name, cofactor names, and rule as SMARTS
        :type rule_text: basestring
        :return:
        :rtype:
        """
        # TODO: replace with csv function
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
        self.rxn_rules[split_text[0]] = (rxn, {'Reactants': reactant_names, 'SMARTS': split_text[2], 'Name': split_text[0],
                                         "Reactions_predicted": 0, "_id": hashlib.sha256(split_text[0].encode()).hexdigest()})

    def load_compound_set(self, compound_file=None, structure_field='structure', id_field='id', fragmented_mols=False):
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
        :param fragmented_mols: Permit the loading of disconnected molecules
        :type fragmented_mols: bool
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
                    if not fragmented_mols and len(AllChem.GetMolFrags(mol)) > 1:
                        continue
                    smi = AllChem.MolToSmiles(mol, True)
                    id = line[id_field]
                    self._add_compound(id, smi, mol=mol)
                    compound_smiles.append(smi)
        elif self.mine:
            db = MINE(self.mine)
            for compound in db.compounds.find():
                id = compound['_id']
                smi = compound['SMILES']
                self._add_compound(id, smi)
                compound_smiles.append(smi)
        else:
            raise ValueError('No input file or database specified for starting compounds')
        print("%s compounds loaded" % len(compound_smiles))
        return compound_smiles

    def _add_compound(self, id, smi, mol=None):
        """Adds a compound to the internal compound dictionary"""
        if not mol:
            mol = AllChem.MolFromSmiles(smi)
        i_key = AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))
        _id = 'C' + hashlib.sha1(smi.encode('utf-8')).hexdigest()
        self.compounds[smi] = {'ID': id, '_id': _id, "SMILES": smi, 'Inchikey': i_key, 'Generation': self.generation,
                               'Reactant_in': [], 'Product_of': [], "Sources": []}
        # if we are building a mine and generating images, do so here
        if self.image_dir and self.mine:
            try:
                MolToFile(mol, os.path.join(self.image_dir, _id + '.svg'), kekulize=False)
            except OSError:
                print("Unable to generate image for %s" % smi)
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
        # these sets collect predictions if this function is being used in isolation
        pred_rxns = set()
        pred_compounds = set()
        if not rules:
            rules = self.rxn_rules.keys()
        mol = AllChem.MolFromSmiles(compound_SMILES)
        mol = AllChem.RemoveHs(mol)
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
            reactant_mols = tuple([self.cofactors[x] for x in rule[1]['Reactants']])
            try:
                product_sets = rule[0].RunReactants(reactant_mols)
            except RuntimeError:  # This error should be addressed in a new version of RDKit
                print("Runtime ERROR!"+rule_name)
                print(compound_SMILES)
                continue
            reactants = next(self._make_compound_tups(reactant_mols, rule_name))  # no enumeration for reactants
            reactants.sort()  # By sorting the reactant (and later products) we ensure that reactions are unique
            for product_mols in product_sets:
                try:
                    for stereo_prods in self._make_compound_tups(product_mols, rule_name, split_stereoisomers=self.split_stereoisomers):
                        pred_compounds.update(x.compound for x in stereo_prods)
                        stereo_prods.sort()
                        # TODO: extract this into _add_reaction method
                        text_rxn = ' + '.join(['(%s) %s' % (x.stoich, x.compound) for x in reactants]) + ' => ' + \
                           ' + '.join(['(%s) %s' % (x.stoich, x.compound) for x in stereo_prods])
                        #this hash function is less informative that the one that appears later, but much faster.
                        rhash = hashlib.sha256(text_rxn.encode()).hexdigest()
                        pred_rxns.add(text_rxn)
                        if rhash not in self.reactions:
                            reaction_data = {"_id": rhash, "Reactants": reactants, "Products": stereo_prods,
                                             "Operators": {rule_name}, "SMILES_rxn": text_rxn,
                                             "Generation": self.generation}
                            self.reactions[rhash] = reaction_data
                        else:
                            self.reactions[rhash]['Operators'].add(rule_name)
                except ValueError:
                    continue
        return pred_compounds, pred_rxns

    def _make_compound_tups(self, mols, rule_name, split_stereoisomers=False):
        """Takes a list of mol objects and returns an generator for (compound, stoich) tuples"""
        # TODO: comment this function & next
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
            yield [stoich_tuple(x, y) if len(x) > 1 else stoich_tuple(x[0], y) for x, y in rxn.items()]

    def _calculate_compound_information(self, raw, mol_obj):
        """Calculate the standard data for a compound & return a tuple with compound_ids. Memoized with _raw_compound
        dict"""
        if raw not in self._raw_compounds:
            if self.explicit_h:
                try:
                    mol_obj = AllChem.RemoveHs(mol_obj)  # this step slows down the process quite a bit
                except:
                    raise ValueError
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
        so this cutoff prevents lag
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
        """Calculates a unique reaction hash using inchikeys. First block is connectivity only, second block is stereo
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
            return "+".join(first_block), "+".join(second_block)

        reactants.sort()
        products.sort()
        r_1, r_2 = __get_blocks(reactants)
        p_1, p_2 = __get_blocks(products)
        first_block = r_1+'<==>'+p_1
        second_block = r_2+'<==>'+p_2

        return hashlib.sha256(first_block.encode()).hexdigest()+"-"+hashlib.md5(second_block.encode()).hexdigest()

    def assign_ids(self):
        """Assigns a numerical ID to compounds for ease of reference. Unique only to the CURRENT run."""
        # If we were running a multiprocess expansion, this removes the dicts from Manger control
        self.compounds = dict(self.compounds)
        self.reactions = dict(self.reactions)
        i = 1
        for comp in sorted(self.compounds.values(), key=lambda x: (x['Generation'], x['_id'])):
            if not comp['ID']:
                comp['ID'] = 'pk_cpd'+str(i).zfill(7)
                i += 1
                self.compounds[comp['SMILES']] = comp
                # if we are not loading into the mine, we generate the image here
                if self.image_dir and not self.mine:
                    mol = AllChem.MolFromSmiles(comp['SMILES'])
                    try:
                        MolToFile(mol, os.path.join(self.image_dir, comp['ID'] + '.png'), fitImage=True, kekulize=False)
                    except OSError:
                        print("Unable to generate image for %s" % comp['SMILES'])
        i = 1
        for rxn in sorted(self.reactions.values(), key=lambda x: (x['Generation'], x['_id'])):
            rxn['ID_rxn'] = ' + '.join(['(%s) %s[c0]' % (x.stoich, self.compounds[x.compound]["ID"]) for x in rxn["Reactants"]]) \
                            + ' => ' + ' + '.join(['(%s) %s[c0]' % (x.stoich, self.compounds[x.compound]["ID"]) for x in rxn["Products"]])
            rxn['ID'] = 'pk_rxn'+str(i).zfill(7)
            i += 1
            self.reactions[rxn['_id']] = rxn

    def transform_all(self, num_workers=1, max_generations=1):
        """
        This function applies all of the reaction rules to all the compounds until the generation cap is reached.

        :param num_workers: The number of CPUs to for the expansion process.
        :type num_workers: int
        :param max_generations: The maximum number of times an reaction rule may be applied
        :type max_generations: int
        :return:
        :rtype:
        """
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
                    # The compound and reaction processes will be shared across processes and thus need a Manager to
                    # handle requests. This can be a bottleneck if running on a server with a large number of CPUs but
                    # that's not the predominate use case
                    manager = multiprocessing.Manager()
                    self.compounds = manager.dict(self.compounds)
                    self.reactions = manager.dict(self.reactions)
                    for i, res in enumerate(pool.imap_unordered(self.transform_compound, compound_smiles)):
                        if not (i+1) % print_on:
                            print("Generation %s: %s percent complete" % (self.generation, round((i+1) / len(compound_smiles)) * 100))
                else:
                    for i, smi in enumerate(compound_smiles):
                        self.transform_compound(smi)
                        if not (i+1) % print_on:
                            print("Generation %s: %s percent complete" % (self.generation, round((i+1) / len(compound_smiles)) * 100))
            print("Generation %s produced %s new compounds and %s new reactions in %s sec" %
                  (self.generation, len(self.compounds)-n_comps, len(self.reactions) - n_rxns, time.time()-ti))

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
            for rxn in sorted(self.reactions.values(), key=lambda x: x['ID']):
                outfile.write(delimiter.join([rxn['ID'], '', rxn['ID_rxn'], rxn["SMILES_rxn"], rxn['_id'],
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
        bulk_c = db.compounds.initialize_unordered_bulk_op()
        bulk_r = db.reactions.initialize_unordered_bulk_op()

        # This loop performs 4 functions to reactions 1. convert stoich_tuples to dicts with hashes 2. add reaction
        # links to comps 3. add source information to compounds 4. iterate the reactions predicted for each relevant
        # operator
        for rxn in self.reactions.values():
            _tmpr, _tmpc = [], []  # having temp variables for the lists avoids pointer issues
            for i, x in enumerate(rxn['Reactants']):
                self.compounds[x.compound]['Reactant_in'].append(rxn['_id'])
                _tmpc.append({"stoich": x.stoich, "c_id": "C"+hashlib.sha1(x.compound.encode('utf-8')).hexdigest()})
            rxn['Reactants'] = _tmpc
            for i, x in enumerate(rxn['Products']):
                self.compounds[x.compound]['Product_of'].append(rxn['_id'])
                self.compounds[x.compound]['Sources'].append(
                    {"Compounds": [x['c_id'] for x in rxn['Reactants']], "Operators": list(rxn["Operators"])})
                _tmpr.append({"stoich": x.stoich, "c_id": "C" + hashlib.sha1(x.compound.encode('utf-8')).hexdigest()})
            rxn["Products"] = _tmpc
            # iterate the number of reactions predicted
            for op in rxn['Operators']:
                self.rxn_rules[op][1]['Reactions_predicted'] += 1
            rxn = utils.convert_sets_to_lists(rxn)
            bulk_r.find({'_id': rxn['_id']}).upsert().replace_one(rxn)
        bulk_r.execute()
        db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Reactions Inserted"})

        for comp_dict in self.compounds.values():
            db.insert_compound(AllChem.MolFromSmiles(comp_dict['SMILES']), comp_dict, bulk=bulk_c)
        bulk_c.execute()
        db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Compounds Inserted"})

        for x in self.rxn_rules.values():
            db.operators.save(x[1])  # there are fewer operators so bulk operations are not really faster
        db.build_indexes()


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
                                                                            "compatibility with bnice operators.")
    parser.add_argument('-m', '--max_workers', default=1, type=int, help="Set the nax number of processes to spawn to "
                                                                         "perform calculations.")
    parser.add_argument('-g', '--generations', default=1, type=int, help="Set the numbers of time to apply the reaction"
                                                                         " rules to the compound set.")
    parser.add_argument('-i', '--image_dir', default=None, help="Specify a directory to store images of all created "
                                                                "compounds")
    options = parser.parse_args()
    pk = Pickaxe(cofactor_list=options.cofactor_list, rule_list=options.rule_list, raceimze=options.raceimize,
                 errors=options.verbose, explicit_h=options.bnice, kekulize=options.bnice, image_dir=options.image_dir,
                 database=options.database)
    if options.image_dir and not os.path.exists(options.image_dir):
        os.mkdir(options.image_dir)
    if options.smiles:
        pk.generation = 0
        pk._add_compound("Start", options.smiles)
    else:
        pk.load_compound_set(compound_file=options.compound_file)
    pk.transform_all(max_generations=options.generations, num_workers=options.max_workers)
    if options.database:
        print("Saving results to %s" % options.database)
        pk.save_to_MINE(options.database)
    else:
        pk.assign_ids()
        pk.write_compound_output_file(options.output_dir+'/compounds.tsv')
        pk.write_reaction_output_file(options.output_dir+'/reactions.tsv')

    print("Execution took %s seconds." % (time.time()-t1))
