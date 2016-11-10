from rdkit.Chem import AllChem
import hashlib
import collections
from os import path

"""Utils.py: contains basic functions reused in various contexts in other modules"""

stoich_tuple = collections.namedtuple("stoich_tuple", 'stoich,compound')


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


def get_dotted_field(input_dict, accessor_string):
    """Gets data from a dictionary using a dotted accessor-string"""
    current_data = input_dict
    for chunk in accessor_string.split('.'):
        current_data = current_data.get(chunk, {})
    return current_data


def save_dotted_field(accessor_string, data):
    """Gets data from a dictionary using a dotted accessor-string"""
    for chunk in accessor_string.split('.')[::-1]:
        data = {chunk: data}
    return data


def memoize(f):
    """ Memoization decorator for a function taking one or more arguments. """
    class memodict(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            ret = self[key] = f(*key)
            return ret
    return memodict().__getitem__


def prevent_overwrite(write_path):
    """
    Prevents overwrite of existing output files by appending "_new" when needed
    :param write_path: potential write path
    :type write_path: string
    :return:
    :rtype:
    """
    while path.exists(write_path):
        sp = write_path.split('.')
        if len(sp) > 1:
            sp[-2] += '_new'
            write_path = '.'.join(sp)
        else:
            write_path += '_new'
    return write_path


def approximate_matches(list1, list2, epsilon=0.01):
    """
    Takes two list of tuples and searches for matches of tuples first value within the supplied epsilon. Emits tuples
    with the tuples second values where found. if a value in one dist does not match the other list, it is emitted alone.
    :param list1: first list of tuples
    :type list1: list
    :param list2: second list of tuples
    :type list2: list
    :param epsilon: maximum difference in
    :type epsilon: float
    :return: second values of tuples
    :rtype: generator
    """
    list1.sort()
    list2.sort()
    list1_index = 0
    list2_index = 0

    while list1_index < len(list1) or list2_index < len(list2):
        if list1_index == len(list1):
            yield (0, list2[list2_index][1])
            list2_index += 1
            continue
        if list2_index == len(list2):
            yield (list1[list1_index][1], 0)
            list1_index += 1
            continue

        list1_element = list1[list1_index][0]
        list2_element = list2[list2_index][0]

        difference = abs(list1_element - list2_element)

        if difference < epsilon:
            yield (list1[list1_index][1], list2[list2_index][1])
            list1_index += 1
            list2_index += 1
        elif list1_element < list2_element:
            yield (list1[list1_index][1], 0)
            list1_index += 1
        elif list2_element < list1_element:
            yield (0, list2[list2_index][1])
            list2_index += 1
        else:
            raise AssertionError('Unexpected else taken')


def dict_merge(finaldict, sourcedict):
    """Merges two dictionaries using sets to avoid duplication of values"""
    for key, val in sourcedict.items():
        if (key in finaldict) and isinstance(finaldict[key], list):
            finaldict[key] = set(finaldict[key])
        if isinstance(val, list):
            if key in finaldict:
                finaldict[key].update(val)
            else:
                finaldict[key] = set(val)
        elif isinstance(val, str):
            if key in finaldict:
                finaldict[key].update(val)
            else:
                finaldict[key] = set(val)
                finaldict[key].update(val)
        elif isinstance(val, float):
            if key not in finaldict:
                finaldict[key] = val
        elif isinstance(val, dict):
            if not key in finaldict:
                finaldict[key] = {}
            dict_merge(finaldict[key], val)


def rxn2hash(reactants, products, return_text=False):
    """Hashes reactant and product lists"""
    def to_str(half_rxn):
        return ['(%s) %s' % (x[0], x[1]) if (len(x) == 2 and not isinstance(x, str)) else '(1) %s' % x
                for x in sorted(half_rxn)]

    text_rxn = ' + '.join(to_str(reactants)) + ' => ' + ' + '.join(to_str(products))
    rhash = hashlib.sha256(text_rxn.encode()).hexdigest()
    if return_text:
        return rhash, text_rxn
    else:
        return rhash


def parse_text_rxn(rxn, rp_del, cp_del, translation_dict=None):
    """Makes a list of product and reactant stoich_tuples"""

    def parse_half(half_rxn, td):
        if translation_dict:
            return [stoich_tuple(1, td[x.strip()]) if len(x.split()) == 1
                    else stoich_tuple(int(x.split()[0].strip('()')), td[x.split()[1].strip()])
                    for x in half_rxn.split(cp_del)]
        else:
            return [stoich_tuple(1, x.strip()) if len(x.split()) == 1
                    else stoich_tuple(int(x.split()[0].strip('()')), x.split()[1].strip())
                    for x in half_rxn.split(cp_del)]

    return [parse_half(x, translation_dict) for x in rxn.split(rp_del)]


_reactions = None
def neutralise_charges(mol, reactions=None):
    def _InitialiseNeutralisationReactions():
        patts = (
            # Imidazoles
            ('[n+;H]', 'n'),
            # Amines
            ('[N+;!H0]', 'N'),
            # Carboxylic acids and alcohols
            ('[$([O-]);!$([O-][#7])]', 'O'),
            # Thiols
            ('[S-;X1]', 'S'),
            # Sulfonamides
            ('[$([N-;X2]S(=O)=O)]', 'N'),
            # Enamines
            ('[$([N-;X2][C,N]=C)]', 'N'),
            # Tetrazoles
            ('[n-]', '[nH]'),
            # Sulfoxides
            ('[$([S-]=O)]', 'S'),
            # Amides
            ('[$([N-]C=O)]', 'N'),

        )
        return [(AllChem.MolFromSmarts(x), AllChem.MolFromSmiles(y, False)) for x, y in patts]

    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol

def do_profile(func):
    from line_profiler import LineProfiler

    def profiled_func(*args, **kwargs):
        try:
            profiler = LineProfiler()
            profiler.add_function(func)
            profiler.enable_by_count()
            return func(*args, **kwargs)
        finally:
            profiler.print_stats()
    return profiled_func