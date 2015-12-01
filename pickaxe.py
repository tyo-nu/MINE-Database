__author__ = 'JGJeffryes'

from rdkit.Chem import AllChem
from databases import MINE
from argparse import ArgumentParser
import itertools
import collections
import time
import utils
from copy import deepcopy
import datetime
import hashlib

class Pickaxe:
    """
        This class generates new compounds from user-specified starting compounds using a set of SMARTS-based reaction
        rules. It may be initialized with a text file containing the reaction rules and cofactors or this may be
        done on an ad hock basis.
    """
    def __init__(self, rule_list=None, cofactor_list=None, explicit_h=True, strip_cof=True, kekulize=True, errors=True,
                 raceimze=False, mine=None):
        self.rxn_rules = {}
        self.cofactors = {}
        self._raw_compounds = {}
        self.compounds = collections.OrderedDict()
        self.reactions = collections.OrderedDict()
        self.mine = mine
        self.generation = 1
        self.explicit_h = explicit_h
        self.strip_cof = strip_cof
        self.kekulize = kekulize
        self.raceimize = raceimze
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
        from rdkit import RDLogger
        lg = RDLogger.logger()
        if not errors:
            lg.setLevel(4)

    def _load_cofactor(self, cofactor_text):
        """
        Loads a cofactor into the cofactor dictionary from a tab-delimited string
        :param cofactor_text: tab-delimited string with the compound name and SMILES
        :type cofactor_text: basestring
        :return:
        :rtype:
        """
        split_text = cofactor_text.strip().split('\t')
        # TODO: add input validation
        mol = AllChem.MolFromSmiles(split_text[1])
        self.compounds[split_text[0]] = {'_id': split_text[0], 'SMILES': split_text[1], "Inchikey": 'None', 'Generation': -1}
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        if self.kekulize:
            AllChem.Kekulize(mol)
            # also need to unset the aromatic flags in case the ring is hydrolysed (else will throw errors)
            for atom in mol.GetAtoms():
                atom.SetIsAromatic(False)
        self.cofactors[split_text[0]] = mol
        self._raw_compounds[split_text[1]] = split_text[0]

    def load_rxn_rule(self, rule_text):
        """
        Loads a reaction rule into the rxn_rule dictionary from a tab-delimited string
        :param rule_text: tab-delimited string with the rule name, cofactor names, and rule as SMARTS
        :type rule_text: basestring
        :return:
        :rtype:
        """
        # TODO: add input validation
        split_text = rule_text.strip().split('\t')
        reactant_names = split_text[1].split(';')
        for name in reactant_names:
            if name == "Any":
                continue
            if name not in self.cofactors:  # try to proceed as if name is SMILES
                self._load_cofactor(name+"\t"+name)
        rxn = AllChem.ReactionFromSmarts(split_text[2])
        if rxn.GetNumReactantTemplates() != len(reactant_names):
            raise ValueError("Number of cofactors does not match supplied reaction rule")
        self.rxn_rules[split_text[0]] = (reactant_names, rxn)

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
        if self.kekulize:
            AllChem.Kekulize(mol)
            # also need to unset the aromatic flags in case the ring is hydrolysed (else will throw errors)
            for atom in mol.GetAtoms():
                atom.SetIsAromatic(False)
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
                continue
            reactants = next(self._make_compound_tups(reactant_mols, rule_name))  # no enumeration for reactants
            for product_mols in product_sets:
                try:
                    for stereo_prods in self._make_compound_tups(product_mols, rule_name, split_stereoisomers=True):
                        pred_compounds.update(x.compound for x in stereo_prods)
                        rid = self._calculate_rxn_hash(reactants, stereo_prods)
                        text_rxn = ' + '.join(['%s "%s"' % (x.stoich, x.compound) for x in reactants]) + ' --> ' + \
                           ' + '.join(['%s "%s"' % (x.stoich, x.compound) for x in stereo_prods])
                        pred_rxns.add(text_rxn)
                        if rid not in self.reactions:
                            reaction_data = {"_id": rid, "Reactants": reactants, "Products": stereo_prods,
                                         "Text_rxn": text_rxn, "Operators": {rule_name}}
                            self.reactions[rid] = reaction_data
                        else:
                            self.reactions[rid]['Operators'].add(rule_name)
                except ValueError:
                    continue
        return pred_compounds, pred_rxns

    def _make_compound_tups(self, mols, rule_name, split_stereoisomers=False):
        """Takes a list of mol objects and returns an generator for (compound, stoich) tuples"""
        comp_tuple = collections.namedtuple("stoich_tuple", 'compound,stoich')
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
            yield [comp_tuple(x, y) if len(x) > 1 else comp_tuple(x[0], y) for x, y in rxn.items()]

    def _calculate_compound_information(self, raw, mol_obj):
        """Calculate the standard data for a compound & return a tuple with compound_ids. Memoized with _raw_compound
        dict"""
        if raw not in self._raw_compounds:
            if self.explicit_h:
                mol_obj = AllChem.RemoveHs(mol_obj) # this step slows down the process quite a bit
            if self.raceimize:
                mols = self._racemization(mol_obj)
            else:
                mols = [mol_obj]
            smiles = [AllChem.MolToSmiles(m, True) for m in mols]
            for s, m in zip(smiles, mols):
                if s not in self._raw_compounds:
                    cid = str(len(self.compounds)+1).zfill(7)
                    try:
                        i_key = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                        self.compounds[cid] = {'_id': cid, "SMILES": s, 'Inchikey': i_key, 'Generation': self.generation}
                    except ValueError:
                        self.compounds[cid] = {'_id': cid, "SMILES": s, 'Inchikey': "None", 'Generation': self.generation}
                    self._raw_compounds[s] = cid
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
                split_inchikey = self.compounds[x.compound]["Inchikey"].split('-')
                if len(split_inchikey) > 1:
                    first_block.append("%s,%s" % (x.stoich, split_inchikey[0]))
                    second_block.append("%s,%s" % (x.stoich, split_inchikey[1]))
            return "+".join(first_block), "+".join(first_block)

        reactants.sort()
        products.sort()
        r_1, r_2 = __get_blocks(reactants)
        p_1, p_2 = __get_blocks(products)
        first_block = r_1+'<==>'+p_1
        second_block = r_2+'<==>'+p_2

        return hashlib.sha256(first_block.encode()).hexdigest()+"-"+hashlib.md5(second_block.encode()).hexdigest()

    def

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
            outfile.write('_id\tSMILES\n')
            for c in sorted(self.compounds.values(), key=lambda x: x['_id']):
                outfile.write(delimiter.join([c['_id'], c['Inchikey'], c['SMILES']])+'\n')

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
            outfile.write('_id\tText Rxn\tOperator\n')
            for rxn in self.reactions.values():
                outfile.write(delimiter.join([rxn['_id'], rxn["Text_rxn"], ';'.join(rxn['Operators'])])+'\n')

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
            rxn = utils.convert_sets_to_lists(rxn)
            db.reactions.save(rxn)
        db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action": "Reactions Inserted"})


if __name__ == "__main__":
    t1 = time.time()
    parser = ArgumentParser()
    parser.add_argument('-c', '--cofactor_list', default="Tests/Cofactor_SMILES.tsv", help="Specify a list of cofactors"
                                                                                           " as a tab-separated file")
    parser.add_argument('-r', '--rule_list', default="Tests/operators_smarts.tsv", help="Specify a list of reaction "
                                                                                        "rules as a tab-separated file")
    parser.add_argument('-o', '--output_dir', default=".", help="The directory in which to place files")
    parser.add_argument('-d', '--database', default=None, help="The name of the database in which to store output. If "
                                                               "not specified, data is still written as tsv files")
    parser.add_argument('-R', '--raceimize', action='store_true', default=False, help="Enumerate the possible chiral "
                        "forms for all unassigned steriocenters in compounds & reactions")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="Display RDKit errors & warnings")
    parser.add_argument('--bnice', action='store_true', default=False, help="Set several options to enable "
                                                                            "compatability with bnice operators.")
    options = parser.parse_args()

    pk = Pickaxe(cofactor_list=options.cofactor_list, rule_list=options.rule_list, raceimze=options.raceimize,
                 errors=options.verbose, explicit_h=options.bnice, kekulize=options.bnice)
    seed_comps = []
    with open('iAF1260.tsv') as infile:
        for i, line in enumerate(infile):
            mol = AllChem.MolFromInchi(line.split('\t')[6])
            smi = AllChem.MolToSmiles(mol, True)
            id = line.split('\t')[0]
            i_key = AllChem.InchiToInchiKey(line.split('\t')[6])
            pk.compounds[id] = {'_id': id, "SMILES": smi, 'Inchikey':i_key, 'Genration': 0}
            pk._raw_compounds[smi] = id
            seed_comps.append(smi)
    for i, smi in enumerate(seed_comps):
        print(i)
        prod, rxns = pk.transform_compound(smi, rules=[])
    pk.write_compound_output_file(options.output_dir+'/compounds.tsv')
    pk.write_reaction_output_file(options.output_dir+'/reactions.tsv')
    if options.database:
        pk.save_to_MINE(options.database)
    print("Exicution took %s seconds." % (time.time()-t1))
