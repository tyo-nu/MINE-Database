__author__ = 'JGJeffryes'

from rdkit.Chem import AllChem, Draw
from itertools import chain, combinations
from collections import defaultdict
from utils import memoize

from line_profiler import LineProfiler
def do_profile(func):
        def profiled_func(*args, **kwargs):
            try:
                profiler = LineProfiler()
                profiler.add_function(func)
                profiler.enable_by_count()
                return func(*args, **kwargs)
            finally:
                profiler.print_stats()
        return profiled_func

class Pickaxe:
    """
        This class generates new compounds from user-specified starting compounds using a set of SMARTS-based reaction
        rules. it may be initialized with a text file containing the reaction rules and cofactors or this may be
        done on an ad hock basis.
    """
    def __init__(self, rule_list=None, cofactor_list=None, explicit_h=True, strip_cof=True, kekulize=True, errors=True):
        self.rxn_rules = {}
        self.cofactors = {}
        self.cof_set = {'[H+]'}
        self.explicit_h = explicit_h
        self.strip_cof = strip_cof
        self.kekulize = kekulize
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
        # add input validation
        mol = AllChem.MolFromSmiles(split_text[1])
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        if self.kekulize:
            AllChem.Kekulize(mol)
            for atom in mol.GetAtoms():
                atom.SetIsAromatic(False)
        self.cofactors[split_text[0]] = mol
        self.cof_set.add(split_text[1])

    def load_rxn_rule(self, rule_text):
        """
        Loads a reaction rule into the rxn_rule dictionary from a tab-delimited string
        :param rule_text: tab-delimited string with the rule name, cofactor names, and rule as SMARTS
        :type rule_text: basestring
        :return:
        :rtype:
        """
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

    def transform_compound(self, compound_SMILES, rules=None, reactions=False):
        """
        Perform transformations to a compound returning the products and the roles predicting the conversion
        :param compound_SMILES: The compound on which to operate represented as SMILES
        :type compound_SMILES: string
        :param rules: The names of the reaction rules to apply. If none, all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as SMILES string and predicting rules as a list of strings
        :rtype: list of tuples
        """
        n_compounds = 0
        if not rules:
            rules = self.rxn_rules.keys()
        mol = AllChem.MolFromSmiles(compound_SMILES)
        if self.kekulize:
            AllChem.Kekulize(mol)
            for atom in mol.GetAtoms():
                atom.SetIsAromatic(False)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        self.cofactors['Any'] = mol
        products = defaultdict(set)
        rxns = []
        for rule_name in rules:
            rule = self.rxn_rules[rule_name]
            comps = tuple([self.cofactors[x] for x in rule[0]])
            try:
                ps = rule[1].RunReactants(comps)
            except RuntimeError:
                print("Runtime ERROR!"+rule_name)
                continue
            for compound in chain.from_iterable(ps):
                n_compounds += 1
                try:
                    compound = AllChem.RemoveHs(compound)
                    smiles = AllChem.MolToSmiles(compound, True)
                except ValueError:
                    smiles = AllChem.MolToSmiles(compound, True)
                    print("Product ERROR!:%s %s" % (rule_name, smiles))
                if not self.strip_cof or smiles not in self.cof_set:
                    products[smiles].add(rule_name)
            if reactions:
                reactants = tuple(self._make_compound_tups(comps))
                rxns.extend([(reactants, tuple(self._make_compound_tups(prods))) for prods in ps])
        print("%s compounds" % n_compounds)
        if reactions:
            return [x for x in products.items()], rxns
        else:
            return [x for x in products.items()]

    def _make_compound_tups(self, mols):
        comps = defaultdict(int)
        for m in mols:
            try:
                AllChem.RemoveHs(m)
                inchikey = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                comps[inchikey] += 1
            except (RuntimeError, ValueError):
                return ()
        return [(y, x) for x, y in comps.items()]

if __name__ == "__main__":
    pk = Pickaxe(cofactor_list="Cofactor_SMILES.tsv", rule_list="operators_smarts.tsv")
    operators = defaultdict(int)
    predicted_rxns = set()
    compounds = defaultdict(set)
    with open('iAF1260.tsv') as infile:
        for i,line in enumerate(infile):
            mol = AllChem.MolFromInchi(line.split('\t')[6])
            smi = AllChem.MolToSmiles(mol, True)
            print(i)
            prod = pk.transform_compound(smi)
            for c in prod:
                compounds[c[0]].update(c[1])
                for op in c[1]:
                    operators[op] += 1
            #for r in rxns:
                #predicted_rxns.add(r)
    with open('compounds', 'w') as outfile:
        for compound in compounds:
            outfile.write('%s\t"%s"\n' %(compound, compounds[compound]))
    for x in sorted(operators.keys()):
        print(x, operators[x])
    #with open('flat_reactions', 'w') as outfile:
        #for rxn in predicted_rxns:
            #outfile.write('+'.join([" ".join(x) for x in rxn[0]])+'-->'+'+'.join([" ".join(x) for x in rxn[1]])+'\n')
