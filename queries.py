"""Queries.py: Contains functions which power the API queries"""

__author__ = 'JGJeffryes'

from ast import literal_eval
from rdkit.Chem import AllChem

default_projection = {'SMILES': 1, 'Formula': 1, 'MINE_id': 1, 'Names': 1, 'Inchikey': 1, 'Mass': 1, 'Sources': 1,
                      'Generation': 1, 'NP_likeness': 1}

def quick_search(db, query, search_projection=default_projection):
    """
    This function takes user provided compound identifiers and attempts to find a related database ID

    :param db: A Mongo Database, DB to search
    :param query: String, a MINE id, KEGG code, ModelSEED id, Inchikey or Name
    :param search_projection: Dictionary, The fields which should be returned
    :return:
    """
    results = []
    #check if query already is a _id
    if (len(query) == 41) and (query[0] == 'C'):
        query_field = '_id'
    elif (len(query) == 6) and (query[0] == 'C'):
        query_field = 'DB_links.KEGG'
    elif (len(query) == 8) and (query[0:2] == 'cpd'):
        query_field = 'DB_links.Model_SEED'
    elif len(query.split('-')[0]) == 14 and query.isupper():
        query_field = 'Inchikey'
        query = query.split('-')[0]
    elif query.isdigit():
        query_field = "MINE_id"
        query = int(query)
    else:
        query_field = 'Names'

    if query_field == 'Inchikey':
        results = [x for x in db.compounds.find({query_field: {'$regex': '^'+query}}, search_projection).limit(500)
                   if x['_id'][0] == "C"]
    elif query_field == 'Names':
        results = [x for x in db.compounds.find({"Names": {'$regex': '^'+query+'$', '$options': 'i'}}, search_projection) if x['_id'][0] == "C"]
        if not results:
            cursor = db.compounds.find({"$text": {"$search": query}}, {"score": {"$meta": "textScore"}, 'Formula': 1,
                                                                           'MINE_id': 1, 'Names': 1, 'Inchikey': 1,
                                                                           'SMILES': 1, 'Mass': 1})
            results.extend(x for x in cursor.sort([("score", {"$meta": "textScore"})]).limit(500) if x['_id'][0] == "C")
    else:
        results = [x for x in db.compounds.find({query_field: query}, search_projection).limit(500)
                   if x['_id'][0] == "C"]

    return results

def advanced_search(db, mongo_query, search_projection=default_projection):
    """
    Returns compounds in the indicated database which match the provided mongo query

    :param db: A Mongo Database, DB to search
    :param mongo_query: String, A valid mongo query
    :param search_projection: Dictionary, The fields which should be returned
    :return:
    """
    if db != 'admin':  # we don't want users poking around here
        query_dict = literal_eval(mongo_query)  # this transforms the string into a dictionary
        return [x for x in db.compounds.find(query_dict, search_projection)]
    else:
        return ['Illegal query']

def similarity_search(db, comp_structure, min_tc, fp_type, limit, search_projection=default_projection):
    """
    Returns compounds in the indicated database which have structural similarity to the provided compound

    :param db: A Mongo Database, DB to search
    :param comp_structure: String, A molecule in Molfile or or SMILES format
    :param min_tc: Float, Minimum Tannimoto score
    :param fp_type: String, Fingerprint type. Currently accepts MACCS or RDKit
    :param limit: Integer, the maximum number of compounds to return
    :param search_projection: Dictionary, The fields which should be returned
    :return:
    """
    similarity_search_results = []
    fp_type = str(fp_type)
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    else:
        mol = AllChem.MolFromSmiles(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    if fp_type == 'MACCS':
        query_fp = set([i for i, bit in enumerate(AllChem.GetMACCSKeysFingerprint(mol)) if bit])
    elif fp_type == 'RDKit':
        query_fp = set([i for i, bit in enumerate(AllChem.RDKFingerprint(mol)) if bit])
    else:
        raise ValueError("Invalid FP_type")

    len_fp = len(query_fp)
    search_projection[fp_type] = 1
    for x in db.compounds.find({"$and": [{"len_"+fp_type: {"$gte": min_tc*len_fp}},
                               {"len_"+fp_type: {"$lte": len_fp/min_tc}}]}, search_projection):
        test_fp = set(x[fp_type])
        tc = len(query_fp & test_fp)/float(len(query_fp | test_fp))
        if tc >= min_tc:
            del x[fp_type]
            similarity_search_results.append(x)
            if len(similarity_search_results) == limit:
                break

    return similarity_search_results

def structure_search(db, comp_structure, search_projection=default_projection):
    """
    Returns compounds in the indicated database which are exact matches to the provided structure

    :param db: A Mongo Database, DB to search
    :param comp_structure: String, A molecule in Molfile or or SMILES format
    :param search_projection: Dictionary, The fields which should be returned
    :return:
    """
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    else:
        mol = AllChem.MolFromSmiles(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    inchi = AllChem.MolToInchi(mol)
    inchi_key = AllChem.InchiToInchiKey(inchi)
    # sure, we could look for a matching SMILES but this is faster
    return quick_search(db, inchi_key, search_projection)

def substructure_search(db, comp_structure, limit, search_projection=default_projection):
    """
    Returns compounds in the indicated database which contain the provided structure

    :param db: A Mongo Database, DB to search
    :param comp_structure: String, A molecule in Molfile or or SMILES format
    :param limit: Integer, the maximum number of compounds to return
    :param search_projection: Dictionary, The fields which should be returned
    :return:
    """
    substructure_search_results = []
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    else:
        mol = AllChem.MolFromSmarts(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")
    query_fp = [i for i, bit in enumerate(AllChem.RDKFingerprint(mol)) if bit]
    for x in db.compounds.find({"RDKit": {"$all": query_fp}}, search_projection):
        comp = AllChem.MolFromSmiles(x['SMILES'])
        if comp.HasSubstructMatch(mol):
            substructure_search_results.append(x)
            if len(substructure_search_results) == limit:
                break
    return substructure_search_results
