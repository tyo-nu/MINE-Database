"""Utils.py: contains basic functions reused in various contexts in other
modules"""
from rdkit.Chem import AllChem
import hashlib
import collections
from os import path

stoich_tuple = collections.namedtuple("stoich_tuple", 'stoich,c_id')


def compound_hash(compound, cofactor=False):
    """Creates hash string for given compound
    
    :param compound: The compound to be hashed
    :type compound: str or Mol Object
    :param cofactor: is the compound a cofactor
    :type cofactor: bool
    :return: A hashed compound _id
    :rtype: str
    """
    # Check to see if compound is a Mol object. If true, convert that Mol object
    # to a SMILES string
    if isinstance(compound, AllChem.Mol):
        compound = AllChem.MolToSmiles(compound, True)
    # Create hash string using hashlib module
    chash = hashlib.sha1(compound.encode('utf-8')).hexdigest()
    # Mark cofactors with an X at the beginning, all else with a C
    if cofactor:
        return "X"+chash
    else:
        return "C"+chash


def convert_sets_to_lists(obj):
    """Recursively converts dictionaries that contain sets to lists"""
    if isinstance(obj, set):
        # This brings short names to the top of the list
        try:
            obj = sorted(list(obj), key=len)
        except TypeError:
            obj = list(obj)
    elif isinstance(obj, dict):
        for key in obj:
            obj[key] = convert_sets_to_lists(obj[key])
    return obj


def get_dotted_field(input_dict, accessor_string):
    """Gets data from a dictionary using a dotted accessor-string
    
    :param input_dict: A nested dictionary
    :type input_dict: dict
    :param accessor_string: A dotted path description i.e. "DBLinks.KEGG"
    :type accessor_string: str
    :return: The value in the nested dict
    """
    current_data = input_dict
    for chunk in accessor_string.split('.'):
        current_data = current_data.get(chunk, {})
    return current_data


def save_dotted_field(accessor_string, data):
    """Saves data to a dictionary using a dotted accessor-string
    
    :param accessor_string: A dotted path description i.e. "DBLinks.KEGG"
    :type accessor_string: str
    :param data: The value to be stored
    :return: A nested dictionary
    :rtype: dict
    """
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
    """Prevents overwrite of existing output files by appending "_new" when
    needed
    
    :param write_path: potential write path
    :type write_path: string
    :return: new write path
    :rtype: str
    """
    while path.exists(write_path):
        sp = write_path.split('.')
        # Make sure that files without an extension are still valid (otherwise,
        # split would create a list of one string which would give an index
        # error when sp[-2] is called)
        if len(sp) > 1:
            sp[-2] += '_new'
            write_path = '.'.join(sp)
        else:
            write_path += '_new'
    return write_path


def approximate_matches(list1, list2, epsilon=0.01):
    """Takes two list of tuples and searches for matches of tuples first value
    within the supplied epsilon. Emits tuples with the tuples second values
    where found. if a value in one dist does not match the other list, it is
    emitted alone.
    
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
    # Get text reaction to be hashed
    def to_str(half_rxn):
        return ['(%s) %s' % (x[0], x[1]) if (len(x) == 2
                                             and not isinstance(x, str))
                else '(1) %s' % x for x in sorted(half_rxn)]

    text_rxn = ' + '.join(to_str(reactants)) + ' => ' + \
               ' + '.join(to_str(products))
    # Hash text reaction
    rhash = hashlib.sha256(text_rxn.encode()).hexdigest()
    if return_text:
        return rhash, text_rxn
    else:
        return rhash


def _calculate_rxn_hash(db, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is
        connectivity only, second block is stereo only"""

        def __get_blocks(tups):
            first_block, second_block = [], []
            for x in tups:
                comp = db.compounds.find_one({"_id": x.c_id})
                if comp and comp["Inchikey"]:
                    split_inchikey = comp["Inchikey"].split('-')
                    if len(split_inchikey) > 1:
                        first_block.append("%s,%s" %
                                           (x.stoich, split_inchikey[0]))
                        second_block.append("%s,%s" %
                                            (x.stoich, split_inchikey[1]))
                else:
                    print("No Inchikey for %s" % x.c_id)
            return "+".join(first_block), "+".join(second_block)

        reactants.sort()
        products.sort()
        r_1, r_2 = __get_blocks(reactants)
        p_1, p_2 = __get_blocks(products)
        first_block = r_1 + '<==>' + p_1
        second_block = r_2 + '<==>' + p_2
        return hashlib.sha256(first_block.encode()).hexdigest() + "-" + \
            hashlib.md5(second_block.encode()).hexdigest()


def parse_text_rxn(rxn, rp_del, cp_del, translation_dict=None):
    """Makes a list of product and reactant stoich_tuples"""

    def parse_half(half_rxn, td):
        if translation_dict:
            return [stoich_tuple(1, td[x.strip()]) if len(x.split()) == 1
                    else stoich_tuple(int(x.split()[0].strip('()')),
                                      td[x.split()[1].strip()])
                    for x in half_rxn.split(cp_del)]
        else:
            return [stoich_tuple(1, x.strip()) if len(x.split()) == 1
                    else stoich_tuple(int(x.split()[0].strip('()')),
                                      x.split()[1].strip())
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
        return [(AllChem.MolFromSmarts(x), AllChem.MolFromSmiles(y, False))
                for x, y in patts]

    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions = _InitialiseNeutralisationReactions()
        reactions = _reactions
    for (reactant, product) in reactions:
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol
