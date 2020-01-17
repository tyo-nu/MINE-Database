"""Tests for queries.py using pytest."""
# pylint: disable=redefined-outer-name

import os
from collections import OrderedDict
import minedatabase.utils as utils

DATA_DIR = os.path.dirname(__file__) + '/data/'


def test_file_to_dict_list():
    """
    GIVEN compound files of different types (.tsv, .csv, and .json)
    WHEN the file contents are converted to a list of compound dicts
    THEN check that the list of compound dicts is produced as expected
    """
    res_7 = OrderedDict(
        [('id', 'cpd01211'), ('abbreviation', 'tcynt'),
         ('name', 'Thiocyanate'), ('formula', 'CNS'), ('mass', '58'),
         ('source', 'ModelSEED'), ('structure', 'InChI=1S/CHNS/c2-1-3/h3H'),
         ('charge', '-1'), ('is_core', '1'), ('is_obsolete', '0'),
         ('linked_compound', 'null'), ('is_cofactor', '0'), ('deltag', '22.2'),
         ('deltagerr', '5.68687'), ('pka', '3:0.5'), ('pkb', ''),
         ('abstract_compound', 'null'), ('comprised_of', 'null'),
         ('aliases', 'null')])

    filenames = ['test_compounds.tsv', 'test_compounds.csv',
                 'test_compounds.json']

    for file in filenames:
        res = utils.file_to_dict_list(DATA_DIR + file)
        assert len(res) == 15
        assert res[7] == res_7
