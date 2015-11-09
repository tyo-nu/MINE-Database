__author__ = 'JGJeffryes'

from rdkit.Chem import AllChem, Draw
from itertools import chain, combinations
from collections import defaultdict
from utils import memoize

class Pickaxe:
    """
        This class generates new compounds from user-specified starting compounds using a set of SMARTS-based reaction
        rules. it may be initialized with a text file containing the reaction rules and cofactors or this may be
        done on an ad hock basis.
    """
    def __init__(self, rule_list=None, cofactor_list=None, explicit_h=True, strip_cof=False, kekulize=True):
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

    def transform_compound(self, compound_SMILES, rules=None):
        """
        Perform transformations to a compound returning the products and the roles predicting the conversion
        :param compound_SMILES: The compound on which to operate represented as SMILES
        :type compound_SMILES: string
        :param rules: The names of the reaction rules to apply. If none, all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as SMILES string and predicting rules as a list of strings
        :rtype: list of tuples
        """
        n_rxns = 0
        if not rules:
            rules = self.rxn_rules.keys()
        mol = AllChem.MolFromSmiles(compound_SMILES)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        if self.kekulize:
            AllChem.Kekulize(mol)
        print(AllChem.MolToSmiles(mol, kekuleSmiles=True))
        self.cofactors['Any'] = mol
        products = defaultdict(set)
        for rule_name in rules:
            rule = self.rxn_rules[rule_name]
            comps = tuple([self.cofactors[x] for x in rule[0]])
            ps = rule[1].RunReactants(comps)
            for compound in chain.from_iterable(ps):
                try:
                    smiles = AllChem.MolToSmiles(AllChem.RemoveHs(compound))
                except (RuntimeError, ValueError):
                    smiles = AllChem.MolToSmiles(compound)
                if not self.strip_cof or smiles not in self.cof_set:
                    products[smiles].add(rule_name)
        print("%s reactions" % n_rxns)
        return [x for x in products.items()]

    def _make_compound_tups(self, mols):
        comps = defaultdict(int)
        for m in mols:
            inchikey = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
            comps[inchikey] += 1
        return [(y, x) for x, y in comps.items()]

if __name__ == "__main__":
    pk = Pickaxe(cofactor_list="Cofactor_SMILES.tsv", rule_list="operators_smarts.tsv")
    operators = set()
    for prod in pk.transform_compound('OC(=O)c1ccccc1'):
        operators.update(prod[1])
    print(sorted(operators))
