"""Compound_io.py: Functions to load MINE databases from and dump compounds
into common cheminformatics formats"""
import csv
import datetime
import os
import sys
import collections

from minedatabase import utils
from minedatabase.databases import MINE
from rdkit.Chem import AllChem


def export_sdf(mine_db, dir_path, max_compounds=None):
    """Exports compounds from the database as an MDL SDF file
    
    :param mine_db: a MINE object
    :param dir_path: directory for files
    :param max_compounds: maximum number of compounds per file (defaults to
         unlimited)
    :return:
    """

    # Make sure that all compounds point to all their reactants
    if not mine_db.compounds.find_one({"Product_of": {'$exists': 1}}):
        mine_db.add_rxn_pointers()

    print("Exporting %s compounds from %s as an SDF file" %
          (mine_db.compounds.count(), mine_db.name))
    target = utils.prevent_overwrite(os.path.join(dir_path, mine_db.name)
                                     + "_1.sdf")
    # SDWriter (rdkit) writes Mol objects to SD files
    w = AllChem.SDWriter(target)
    w.SetKekulize(True)
    n_files = 1
    for compound in mine_db.compounds.find():
        # Convert SMILES string to Mol object, replacing 'CoA' and 'R' by '*'
        mol = AllChem.MolFromSmiles(compound['SMILES'], True,
                                    {'CoA': '*', 'R': "*"})
        # if Mol object successfully generated, annotate properties
        if mol:
            mol.SetProp('_id', compound['_id'])
            mol.SetProp('Generation', str(compound['Generation']))
            if 'Reactant_in' in compound:
                mol.SetProp('Reactant_in', str(compound['Reactant_in']))
            if 'Product_of' in compound:
                mol.SetProp('Product_of', str(compound['Product_of']))
            w.write(mol)
            # Start writing a new sdf file if the maximum (set by user) has
            # been reached for the current file
            if max_compounds and (w.NumMols() >= max_compounds):
                n_files += 1
                target = utils.prevent_overwrite(
                    os.path.join(dir_path, mine_db.name) + "_%s.sdf" % n_files)
                w = AllChem.SmilesWriter(target)
    w.close()


def export_smiles(mine_db, dir_path, max_compounds=None):
    """Exports compounds from the database as a SMILES file
    
    :param mine_db: a MINE object
    :param dir_path: directory for files
    :param max_compounds: maximum number of compounds per file (defaults to
        unlimited)
    :return:
    """
    header = ['SMILES', "_id", "Generation", 'Reactant_in', 'Product_of']
    # Make sure that all compounds point to all their reactants
    if not mine_db.compounds.find_one({"Product_of": {'$exists': 1}}):
        mine_db.add_rxn_pointers()

    print("Exporting %s compounds from %s as SMILES file" %
          (mine_db.compounds.count(), mine_db.name))
    target = open(utils.prevent_overwrite(os.path.join(dir_path, mine_db.name)
                                          + "_1.smiles"), 'w')

    # DictWriter allows for each key:value pair of a dictionary to be written
    # on its own row (by writerow)
    w = csv.DictWriter(target, fieldnames=header, dialect='excel-tab')
    n_files = 1
    i = 0
    for compound in mine_db.compounds.find({}, dict([(x, 1) for x in header])):
        w.writerow(compound)
        i += 1
        # If max compounds per file has been set by user and our number of
        # compounds that we have written so far is divisible by the max number,
        # then we start a new file
        if max_compounds and not (i % max_compounds):
            n_files += 1
            target = open(utils.prevent_overwrite(
                os.path.join(dir_path, mine_db.name) + "_%s.smiles" % n_files),
                'w')
            w = csv.DictWriter(target, fieldnames=header, dialect='excel-tab')


def export_mol(mine_db, target, name_field='_id'):
    """Exports compounds from the database as MDL molfiles
    
    :param mine_db: The database to export
    :type mine_db: a MINE object
    :param target: a directory in which to place the files
    :type target: str
    :param name_field: the field to provide names for the mol files. Must be
        unique & universal
    :type name_field: str
    :return:
    :rtype:
    """
    # Create the file if it doesn't yet exist
    if not os.path.exists(target):
        os.mkdir(target)

    # Let user know if an id does not exist for every compound in database
    if mine_db.compounds.find().count() != mine_db.compounds.find(
            {name_field: {'$exists': 1}}).count():
        raise ValueError('%s does not exist for every compound in the database'
                         % name_field)

    for compound in mine_db.compounds.find({'_id': {'$regex': '^C'}}):
        # Create Mol object from SMILES code for each compound using
        # MolFromSmiles (rdkit). Take stereochemistry into account (True),
        # and replace CoA and R with *.
        mol = AllChem.MolFromSmiles(compound['SMILES'], True,
                                    {'CoA': '*', 'R': "*"})
        if "." in name_field:
            compound[name_field] = utils.get_dotted_field(compound, name_field)
        # Make things more compact and look nicer
        if isinstance(compound[name_field], list):
            compound[name_field] = ','.join(compound[name_field])
        # Use MolToMolFile (rdkit) to create a mol file from the Mol object
        # with the file path specified.
        AllChem.MolToMolFile(mol, os.path.join(target,
                                               compound[name_field]+'.mol'))


def export_tsv(mine_db, target, compound_fields=('_id', 'Names', 'Model_SEED',
                                                 'Formula', 'Charge', 'Inchi'),
               reaction_fields=('_id', 'SMILES_rxn', 'C_id_rxn')):
    """Exports MINE compound and reaction data as tab-separated values files
    amenable to use in ModelSEED.
    
    :param mine_db: The database to export
    :type mine_db: a MINE object
    :param target: a directory in which to place the files
    :type target: string
    :param compound_fields: The fields to export in the compound table
    :type compound_fields: set
    :param reaction_fields: The fields to export in the reaction table
    :type reaction_fields: set
    """
    db_links = ('KEGG', 'Model_SEED', 'PubChem')
    print("Exporting %s compounds from %s to tsv" % (mine_db.compounds.count(),
                                                     mine_db.name))
    with open(utils.prevent_overwrite(os.path.join(target, mine_db.name)
                                      + "_compounds.tsv"), 'w') as out:
        w = csv.DictWriter(out, fieldnames=compound_fields, dialect='excel-tab')
        w.writeheader()
        for compound in mine_db.compounds.find(
                {}, dict([('SMILES', 1)]+[('DB_links.'+x, 1) if x in db_links
                                          else (x, 1) for x
                                          in compound_fields])):
            # This is a work around for supporting older MINEs which lack Inchi
            if 'Inchi' in compound_fields and 'Inchi' not in compound:
                compound['Inchi'] = AllChem.MolToInchi(AllChem.MolFromSmiles(
                    compound['SMILES']))
            if 'SMILES' not in compound_fields:
                del compound['SMILES']
            if 'DB_links' in compound:
                for k, v in compound['DB_links'].items():
                    compound[k] = ", ".join(v)
                del compound['DB_links']
            w.writerow(compound)

    print("Exporting %s reactions from %s to tsv" % (mine_db.reactions.count(),
                                                     mine_db.name))
    with open(utils.prevent_overwrite(os.path.join(target, mine_db.name)
                                      + "_reactions.tsv"), 'w') as out:
        w = csv.DictWriter(out, fieldnames=reaction_fields, dialect='excel-tab')
        w.writeheader()
        for rxn in mine_db.reactions.find(
                {}, dict([('Reactants', 1), ('Products', 1)] +
                         [(x, 1) for x in reaction_fields])):
            if 'C_id_rxn' in reaction_fields:
                def to_str(half_rxn):
                    return ['(%s) %s' % (x['stoich'], x['c_id'])
                            for x in half_rxn]
                rxn['C_id_rxn'] = ' + '.join(to_str(rxn['Reactants'])) + \
                                  ' => ' + ' + '.join(to_str(rxn['Products']))
            if 'Reactants' not in reaction_fields:
                del rxn['Reactants']
            if 'Products' not in reaction_fields:
                del rxn['Products']
            w.writerow(rxn)


def export_kbase(mine_db, target):
    """Exports MINE compound and reaction data as tab-separated values files
    amenable to use in ModelSEED.
    
    :param mine_db: The database to export
    :type mine_db: a MINE object
    :param target: a directory in which to place the files
    :type target: str
    """
    compound_fields = collections.OrderedDict([('id', "_id"), ('name', ""),
                                               ('formula','Formula'),
                                               ('charge', 'Charge'),
                                               ('aliases', "Names")])
    reaction_fields = collections.OrderedDict(
        [('id', "_id"), ('direction', ">"), ('compartment', "c0"), ('gpr', ''),
         ('name', ''), ('enzyme', ''), ('pathway', ''), ('reference', ''),
         ('equation', '')])
    print("Exporting %s compounds from %s to tsv" % (mine_db.compounds.count(),
                                                     mine_db.name))
    with open(utils.prevent_overwrite(os.path.join(target, mine_db.name)
                                      + "_compounds.tsv"), 'w') as out:
        w = csv.DictWriter(out, fieldnames=compound_fields,
                           dialect='excel-tab')
        w.writeheader()
        for compound in mine_db.compounds.find(
                {}, dict([("Names", 1), ('DB_links.Model_SEED',1)]
                         + [(x, 1) for x in compound_fields.values()])):
            if compound['_id'][0] == 'X':
                continue
            for k, v in compound_fields.items():
                if v in compound:
                    compound[k] = compound[v]
                    del compound[v]
            if 'name' in compound_fields and 'Names' in compound:
                compound['name'] = compound['Names'][0]
                del compound['Names']
            if 'aliases' in compound:
                compound['aliases'] = "|".join(compound['aliases'])
                if 'Model_SEED' in compound['DB_links']:
                    compound['aliases'] += "|" + "|".join(sorted(
                        compound['DB_links']['Model_SEED']))
            if 'DB_links' in compound:
                del compound['DB_links']
            w.writerow(compound)

    print("Exporting %s reactions from %s to tsv" % (mine_db.reactions.count(),
                                                     mine_db.name))
    with open(utils.prevent_overwrite(os.path.join(target, mine_db.name)
                                      + "_reactions.tsv"), 'w') as out:
        w = csv.DictWriter(out, fieldnames=reaction_fields, dialect='excel-tab')
        w.writeheader()
        for rxn in mine_db.reactions.find(
                {}, dict([('Reactants', 1), ('Products', 1)] +
                         [(x, 1) for x in reaction_fields.values()])):
            for k, v in reaction_fields.items():
                if v in rxn:
                    rxn[k] = rxn[v]
                    del rxn[v]
            if 'equation' in reaction_fields:
                def to_str(half_rxn):
                    return ['(%s) %s' % (x['stoich'],
                                         x['c_id'].replace('X', 'C'))
                            for x in half_rxn]
                rxn['equation'] = ' + '.join(to_str(rxn['Reactants'])) + \
                                  ' => ' + ' + '.join(to_str(rxn['Products']))
            if 'Reactants' not in reaction_fields:
                del rxn['Reactants']
            if 'Products' not in reaction_fields:
                del rxn['Products']
            w.writerow(rxn)


def export_inchi_rxns(mine_db, target, rxn_ids=None):
    reaction_fields = collections.OrderedDict(
        [("Reaction Rule", "Operators"), ('ID', "_id"), ('Equation', '')])
    comp_memo = {}

    def get_name_and_inchi(comp_id):
        if comp_id not in comp_memo:
            comp = mine_db.compounds.find_one({"_id": comp_id},
                                              {"Names": 1, "Inchi": 1, "MINE_id": 1})
            comp_memo[comp_id] = (comp.get('Names', [comp['MINE_id']])[0], comp.get('Inchi'))
        return comp_memo[comp_id]

    def to_str(half_rxn):
        lst = []
        for x in half_rxn:
            name, inchi = get_name_and_inchi(x['c_id'])
            lst.append('(%s) %s[%s]' % (x['stoich'], name, inchi))
        return lst

    with open(utils.prevent_overwrite(os.path.join(target, mine_db.name)
                                      + "_reactions.tsv"), 'w') as out:
        w = csv.DictWriter(out, fieldnames=reaction_fields, dialect='excel-tab')
        w.writeheader()
        if rxn_ids:
            query = {"_id": {"$in": rxn_ids}}
        else:
            query = {}
        for rxn in mine_db.reactions.find(
                query, dict([('Reactants', 1), ('Products', 1)] +
                            [(x, 1) for x in reaction_fields.values()])):
            for k, v in reaction_fields.items():
                if v in rxn:
                    if isinstance(rxn[v], list):
                        rxn[k] = ", ".join(rxn[v])
                    else:
                        rxn[k] = rxn[v]
                    del rxn[v]
            if 'Equation' in reaction_fields:
                rxn['Equation'] = ' + '.join(to_str(rxn['Reactants'])) + \
                                  ' => ' + ' + '.join(to_str(rxn['Products']))
            if 'Reactants' not in reaction_fields:
                del rxn['Reactants']
            if 'Products' not in reaction_fields:
                del rxn['Products']
            w.writerow(rxn)


def import_sdf(mine_db, target,):
    """Imports a SDF file as a MINE database
    
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the SDF file to be loaded
    :return:
    """
    # SDMolSupplier (rdkit) takes entries from sdf file and returns Mol objects
    sdf_gen = AllChem.SDMolSupplier(target)
    # Go through each generated Mol object and add each to MINE database
    for mol in sdf_gen:
        mine_db.insert_compound(mol, compound_dict=mol.GetPropsAsDict(),
                                pubchem_db=None, kegg_db=None,
                                modelseed_db=None)
    # Add to log file (metadata)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action":
                              "SDF Imported", "Filepath": target})


def import_smiles(mine_db, target,):
    """
    Imports a smiles file as a MINE database
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the SDF file to be loaded
    :return:
    """
    # SmilesMolSupplier (rdkit) generates Mol objects from smiles file (.smi)
    mols = AllChem.SmilesMolSupplier(target, delimiter='\t', nameColumn=0)
    # Go through each generated mol file and add molecule to MINE database
    # Stores compound properties in dict (GetPropsAsDict() from rdkit Mol class)
    for mol in mols:
        if mol:
            mine_db.insert_compound(mol, compound_dict=mol.GetPropsAsDict(),
                                    pubchem_db=None, kegg_db=None,
                                    modelseed_db=None)
    # Add to log file (metadata)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action":
                              "SDF Imported", "Filepath": target})


def import_mol_dir(mine_db, target, name_field="Name", overwrite=False):
    """Imports a directory of molfiles as a MINE database
    
    :param mine_db: a MINE object, the database to insert the compound into
    :param target: a path, the molfile directory to be loaded
    :param name_field: a string, the field for the compound name
    :param overwrite: a bool, if true, new compounds replace the old compounds
        in the database
    """
    # For each .mol file in the directory of the target folder (path):
    for file in os.listdir(target):
        if ".mol" in file:
            # MolFromMolFile (rdkit) generates Mol objects from .mol files
            mol = AllChem.MolFromMolFile(target+'/'+file)
            # Mol object name becomes name of mol file without .mol extension
            name = file.rstrip('.mol')
            # Check that Mol object is successfully generated
            if mol:
                # Create hashkey for the compound
                comphash = utils.compound_hash(mol)
                # If we don't want to overwrite, and the compound (comphash)
                # already exists, then add an extra comphash for that molecule
                if not overwrite and mine_db.compounds.count({"_id": comphash}):
                    mine_db.compounds.update({"_id": comphash},
                                             {"$addToSet": {name_field: name}})
                # If we don't care about overwriting, just insert the new
                # compound into the database
                else:
                    mine_db.insert_compound(mol,
                                            compound_dict={name_field: [name],
                                                           'Generation': 0},
                                            pubchem_db=None, kegg_db=None,
                                            modelseed_db=None)
    # Add to log file (metadata)
    mine_db.meta_data.insert({"Timestamp": datetime.datetime.now(), "Action":
                              "MolFiles Imported", "Filepath": target})

if __name__ == '__main__':
    # User inputs task as first argument (export-sdf, export-smi, export-mol,
    #  import-sdf, import-smi, or import-mol)
    task = sys.argv[1]
    # User inputs database name as second argument
    db_name = sys.argv[2]
    # User inputs file path as third argument
    path = sys.argv[3]
    database = MINE(db_name)
    if task == 'export-sdf':
        # If a maximum molecules per file is specified (fourth argument
        # entered by user), then pass that to the export function.
        if len(sys.argv) == 5:
            export_sdf(database, path, int(sys.argv[4]))
        # Otherwise, assume an unlimited number of molecules per file
        else:
            export_sdf(database, path)
    elif task == 'export-smi':
        # If a maximum molecules per file is specified (fourth argument
        # entered by user), then pass that to the export function.
        if len(sys.argv) == 5:
            export_smiles(database, path, int(sys.argv[4]))
        # Otherwise, assume an unlimited number of molecules per file
        else:
            export_smiles(database, path)
    elif task == 'export-mol':
        # If a maximum molecules per file is specified (fourth argument
        # entered by user), then pass that to the export function.
        if len(sys.argv) == 5:
            export_mol(database, path, sys.argv[4])
        # Otherwise, assume an unlimited number of molecules per file
        else:
            export_mol(database, path)
    elif task == 'export-tsv':
        export_tsv(database, path)
    elif task == 'export-kbase':
        export_kbase(database, path)
    elif task == 'import-sdf':
        import_sdf(database, path)
    elif task == 'import-smi':
        import_smiles(database, path)
    elif task == 'import-mol':
        import_mol_dir(database, path)
    else:
        print("ERROR: Unrecognised Task")
