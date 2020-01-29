"""Queries.py: Contains functions which power the API queries"""
import re
from ast import literal_eval

from rdkit.Chem import AllChem

DEFAULT_PROJECTION = {'SMILES': 1, 'Formula': 1, 'MINE_id': 1, 'Names': 1,
                      'Inchikey': 1, 'Mass': 1, 'Sources': 1, 'Generation': 1,
                      'NP_likeness': 1}


def quick_search(db, query, search_projection=DEFAULT_PROJECTION.copy()):
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
    if re.match(r"C\w{40}", query):
        query_field = '_id'
    elif re.match(r"C\d{5}", query):
        query_field = 'DB_links.KEGG'
    elif re.match(r'cpd\d{5}', query):
        query_field = 'DB_links.Model_SEED'
    elif re.search(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", query):
        query_field = 'Inchikey'
        query = query.split("=", 1)[-1]
    elif query.isdigit():
        query_field = "MINE_id"
        query = int(query)
    else:
        query_field = 'Names'

    if query_field == 'Names':
        # Return results for all compounds with specified name
        # Make sure that cofactors are not included
        results = db.compounds.find({"Names": {'$regex': '^' + query + '$',
                                               '$options': 'i'}},
                                    search_projection)
        results = [x for x in results if x['_id'][0] == "C"]
        if not results:
            cursor = db.compounds.find({"$text": {"$search": query}},
                                       {"score": {"$meta": "textScore"},
                                        'Formula': 1, 'MINE_id': 1, 'Names': 1,
                                        'Inchikey': 1, 'SMILES': 1, 'Mass': 1})
            top_x = cursor.sort([("score", {"$meta": "textScore"})]).limit(500)
            results.extend(x for x in top_x if x['_id'][0] == "C")
    else:
        # If query isn't a compound name, then return compounds with
        # specified query
        results = [x for x in db.compounds.find({query_field: query},
                                                search_projection).limit(500)
                   if x['_id'][0] == "C"]

    return results


def advanced_search(db, mongo_query,
                    search_projection=DEFAULT_PROJECTION.copy()):
    """Returns compounds in the indicated database which match the provided
    mongo query

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


def similarity_search(db, comp_structure, min_tc, limit, fp_type='RDKit',
                      search_projection=DEFAULT_PROJECTION.copy()):
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
    for x in db.compounds.find(
            {"$and": [{"len_" + fp_type: {"$gte": min_tc * len_fp}},
                      {"len_" + fp_type: {"$lte": len_fp / min_tc}}]},
            search_projection):
        # Put fingerprint in set for fast union (&) and intersection (|)
        # calculations
        test_fp = set(x[fp_type])
        # Calculate tanimoto coefficient
        tmc = len(query_fp & test_fp) / float(len(query_fp | test_fp))
        # If a sufficient tanimoto coefficient is calculated, append the
        # compound to the search results (until the limit is reached)
        if tmc >= min_tc:
            del x[fp_type]
            similarity_search_results.append(x)
            if len(similarity_search_results) == limit:
                break

    del search_projection[fp_type]
    return similarity_search_results


def structure_search(db, comp_structure, stereo=True,
                     search_projection=DEFAULT_PROJECTION.copy()):
    """Returns compounds in the indicated database which are exact matches to
    the provided structure

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
            {"Inchikey": {'$regex': '^' + inchi_key.split('-')[0]}},
            search_projection)]


def substructure_search(db, sub_structure, limit,
                        search_projection=DEFAULT_PROJECTION.copy()):
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
        # specified substructure. If so, append it to the results (until
        # limit).
        if comp and comp.HasSubstructMatch(mol):
            substructure_search_results.append(x)
            if len(substructure_search_results) == limit:
                break
    return substructure_search_results


def get_ids(db, collection, query):
    """Returns ids for documents in database collection matching query.

    :param db: DB to search
    :type db: A Mongo Database
    :param collection: collection within DB to search
    :type collection: A Mongo Database Collection
    :param query: A valid Mongo query.
    :type query: str
    """

    if query:
        query = literal_eval(query)
    else:
        query = {}
    ids = [x['_id'] for x in db[collection].find(query)]

    return ids


def get_comps(db, id_list):
    """Returns compounds with associated IDs from a Mongo database.

    :param db: DB to search
    :type db: A Mongo Database
    :param id_list: IDs to get compounds for.
    :type id_list: list
    :return: Compound JSON documents
    :rtype: list
    """
    compounds = []
    for cpd_id in id_list:
        excluded_fields = {"len_FP2": 0, "FP2": 0, "len_FP4": 0, "FP4": 0}
        if isinstance(cpd_id, int):
            cpd = db.compounds.find_one({'MINE_id': cpd_id}, excluded_fields)
        else:
            cpd = db.compounds.find_one({'_id': cpd_id}, excluded_fields)
        # New MINEs won't have this precomputed
        if cpd and 'Reactant_in' not in cpd and 'Product_of' not in cpd:
            rxns_as_sub = db.reactions.find({'Reactants.c_id': cpd['_id']})
            cpd['Reactant_in'] = [x['_id'] for x in rxns_as_sub]
            rxns_as_prod = db.reactions.find({'Products.c_id': cpd['_id']})
            cpd['Product_of'] = [x['_id'] for x in rxns_as_prod]
        compounds.append(cpd)

    return compounds


def get_rxns(db, id_list):
    """Returns reactions with associated IDs from a Mongo database.

    :param db: DB to search
    :type db: A Mongo Database
    :param id_list: IDs to get compounds for.
    :type id_list: list
    :return: Reaction JSON documents
    :rtype: list
    """
    reactions = []
    for rxn_id in id_list:
        reactions.append(db.reactions.find_one({'_id': rxn_id}))

    return reactions


def get_ops(db, operator_ids):
    """Returns operators from a Mongo database.

    :param db: DB to search
    :type db: A Mongo Database
    :param operator_ids: Mongo _ids or operator names (e.g. 1.1.-1.h)
    :type operator_ids: list
    :return: Operator JSON documents
    :rtype: list
    """
    if not operator_ids:
        operators = [op for op in db.operators.find()]
    else:
        operators = []
        for op_id in operator_ids:
            op = db.operators.find_one({'$or': [{'_id': op_id},
                                                {"Name": op_id}]})
            operators.append(op)

    return operators


def get_op_w_rxns(db, operator_id):
    """Returns operator with all its associated reactions.

    :param db: DB to search
    :type db: A Mongo Database
    :param operator_id: Mongo _id or operator name (e.g. 1.1.-1.h)
    :type operator_id: str
    :return: Operator JSON document (with reactions)
    :rtype: list
    """
    operator = db.operators.find_one({'$or': [{'_id': operator_id},
                                              {"Name": operator_id}]})
    if operator:
        operator['Reaction_ids'] = \
            db.reactions.find({"Operators": operator_id}).distinct("_id")
    else:
        return None

    return operator
