"""Queries.py: Contains functions which power the API queries"""
from ast import literal_eval
from rdkit.Chem import AllChem

default_projection = {'SMILES': 1, 'Formula': 1, 'MINE_id': 1, 'Names': 1,
                      'Inchikey': 1, 'Mass': 1, 'Sources': 1, 'Generation': 1,
                      'NP_likeness': 1}


def quick_search(db, query, search_projection=default_projection.copy()):
    """This function takes user provided compound identifiers and attempts to
     find a related database ID

    :param db: DB to search
    :type db: A Mongo Database
    :param query: A MINE id, KEGG code, ModelSEED id, Inchikey or Name
    :type query: str
    :param search_projection: The fields which should be returned
    :type search_projection: str
    :return: Query results
    :rtype: list
    """

    # Determine what kind of query was input (e.g. KEGG code, MINE id, etc.)
    # If it can't be determined to be an id or a key, assume it is the name
    # of a compound.
    if (len(query) == 41) and (query[0] == 'C'):
        query_field = '_id'
    elif (len(query) == 6) and (query[0] == 'C'):
        query_field = 'DB_links.KEGG'
    elif (len(query) == 8) and (query[0:2] == 'cpd'):
        query_field = 'DB_links.Model_SEED'
    elif len(query.split('-')[0]) == 14 and query.isupper():
        query_field = 'Inchikey'
    elif query.isdigit():
        query_field = "MINE_id"
        query = int(query)
    else:
        query_field = 'Names'

    if query_field == 'Names':
        # Return results for all compounds with specified name
        # Make sure that cofactors are not included
        results = [x for x in db.compounds.find(
                              {"Names": {'$regex': '^'+query+'$', '$options':
                               'i'}}, search_projection) if x['_id'][0] == "C"]
        if not results:
            cursor = db.compounds.find({"$text": {"$search": query}},
                                       {"score": {"$meta": "textScore"},
                                        'Formula': 1, 'MINE_id': 1, 'Names': 1,
                                        'Inchikey': 1, 'SMILES': 1, 'Mass': 1})
            results.extend(x for x in cursor.sort(
                [("score",
                 {"$meta": "textScore"})]).limit(500) if x['_id'][0] == "C")
    else:
        # If query isn't a compound name, then return compounds with
        # specified query
        results = [x for x in db.compounds.find({query_field: query},
                   search_projection).limit(500) if x['_id'][0] == "C"]

    return results


def advanced_search(db, mongo_query,
                    search_projection=default_projection.copy()):
    """Returns compounds in the indicated database which match the provided mongo
    query

    :param db: DB to search
    :type db: A Mongo Database
    :param query: A valid mongo query
    :type query: str
    :param search_projection: The fields which should be returned
    :type search_projection: str
    :return: Query results
    :rtype: list
    """
    # We don't want users poking around here
    if db.name == 'admin' or not mongo_query:
        raise ValueError('Illegal query')
    # This transforms the string into a dictionary
    query_dict = literal_eval(mongo_query)
    return [x for x in db.compounds.find(query_dict, search_projection)]


def similarity_search(db, comp_structure, min_tc, fp_type, limit,
                      search_projection=default_projection.copy()):
    """Returns compounds in the indicated database which have structural
     similarity to the provided compound

    :param db: DB to search
    :type db: A Mongo Database
    :param comp_structure: A molecule in Molfile or SMILES format
    :type comp_structure: str
    :param min_tc: Minimum Tanimoto score
    :type min_tc: float
    :param fp_type: Fingerprint type. Currently accepts MACCS or RDKit
    :type fp_type: str
    :param limit: The maximum number of compounds to return
    :type limit: int
    :param search_projection: The fields which should be returned
    :type search_projection: str
    :return: Query results
    :rtype: list
    """
    similarity_search_results = []
    fp_type = str(fp_type)
    # Create Mol object from Molfile (has newlines)
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    # Create Mol object from SMILES string (does not have newlines)
    else:
        mol = AllChem.MolFromSmiles(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    # Based on fingerprint type specified by user, get the finger print as an
    # explicit bit vector (series of 1s and 0s). Then, return a set of all
    # indices where a bit is 1 in the bit vector.
    if fp_type == 'MACCS':
        query_fp = set(AllChem.GetMACCSKeysFingerprint(mol).GetOnBits())
    elif fp_type == 'RDKit':
        query_fp = set(AllChem.RDKFingerprint(mol).GetOnBits())
    else:
        raise ValueError("Invalid FP_type")

    len_fp = len(query_fp)
    # Return only id and fingerprint vector
    search_projection[fp_type] = 1
    # Filter compounds that meet tanimoto coefficient size requirements
    for x in db.compounds.find({"$and": [{"len_"+fp_type: {"$gte": min_tc*len_fp
                                                           }},
                               {"len_"+fp_type: {"$lte": len_fp/min_tc}}]},
                               search_projection):
        # Put fingerprint in set for fast union (&) and intersection (|)
        # calculations
        test_fp = set(x[fp_type])
        # Calculate tanimoto coefficient
        tc = len(query_fp & test_fp)/float(len(query_fp | test_fp))
        # If a sufficient tanimoto coefficient is calculated, append the
        # compound to the search results (until the limit is reached)
        if tc >= min_tc:
            del x[fp_type]
            similarity_search_results.append(x)
            if len(similarity_search_results) == limit:
                break

    del search_projection[fp_type]
    return similarity_search_results


def structure_search(db, comp_structure, stereo=True,
                     search_projection=default_projection.copy()):
    """Returns compounds in the indicated database which are exact matches to the
    provided structure

    :param db: DB to search
    :type db: A Mongo Database
    :param comp_structure: A molecule in Molfile or SMILES format
    :type comp_structure: str
    :param stereo: If true, uses stereochemistry in finding exact match
    :type stereo: bool
    :param search_projection: The fields which should be returned
    :type search_projection: str
    :return: Query results
    :rtype: list    """

    # Create Mol object from Molfile (has newlines)
    if "\n" in comp_structure:
        mol = AllChem.MolFromMolBlock(str(comp_structure))
    # Create Mol object from SMILES string (does not have newlines)
    else:
        mol = AllChem.MolFromSmiles(str(comp_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")

    # Get InChI string from mol file (rdkit)
    inchi = AllChem.MolToInchi(mol)
    # Get InChI key from InChI string (rdkit)
    inchi_key = AllChem.InchiToInchiKey(inchi)
    # Sure, we could look for a matching SMILES but this is faster
    if stereo:
        return quick_search(db, inchi_key, search_projection)
    else:
        return [x for x in db.compounds.find(
            {"Inchikey": {'$regex': '^'+inchi_key.split('-')[0]}},
            search_projection)]


def substructure_search(db, sub_structure, limit,
                        search_projection=default_projection.copy()):
    """Returns compounds in the indicated database which contain the provided
    structure

    :param db: DB to search
    :type db: A Mongo Database
    :param comp_structure: A molecule in Molfile or SMILES format
    :type comp_structure: str
    :param limit: The maximum number of compounds to return
    :type limit: int
    :param search_projection: The fields which should be returned
    :type search_projection: str
    :return: Query results
    :rtype: list
    """

    substructure_search_results = []
    # Create Mol object from Molfile (has newlines)
    if "\n" in sub_structure:
        mol = AllChem.MolFromMolBlock(str(sub_structure))
    # Create Mol object from SMILES string (does not have newlines)
    else:
        mol = AllChem.MolFromSmiles(str(sub_structure))
    if not mol:
        raise ValueError("Unable to parse comp_structure")
    # Based on fingerprint type specified by user, get the finger print as an
    # explicit bit vector (series of 1s and 0s). Then, return a set of all
    # indices where a bit is 1 in the bit vector.
    query_fp = list(AllChem.RDKFingerprint(mol).GetOnBits())
    for x in db.compounds.find({"RDKit": {"$all": query_fp}},
                               search_projection):
        # Get Mol object from SMILES string (rdkit)
        comp = AllChem.MolFromSmiles(x['SMILES'])
        # Use HasSubstructMatch (rdkit) to determine if compound has a
        # specified substructure. If so, append it to the results (until limit).
        if comp and comp.HasSubstructMatch(mol):
            substructure_search_results.append(x)
            if len(substructure_search_results) == limit:
                break
    return substructure_search_results
