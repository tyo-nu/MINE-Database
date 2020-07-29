"""Pickaxe.py: This module generates new compounds from user-specified starting
   compounds using a set of SMARTS-based reaction rules."""
import collections
import csv
import datetime
import hashlib
import itertools
import multiprocessing
import os
import re
import time
from sys import exit
from argparse import ArgumentParser
from copy import deepcopy

from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToFile, rdMolDraw2D

from minedatabase import utils
from minedatabase.databases import MINE
import minedatabase.databases as databases
from minedatabase.utils import rxn2hash, StoichTuple, get_fp


class Pickaxe:
    """This class generates new compounds from user-specified starting
    compounds using a set of SMARTS-based reaction rules. It may be initialized
    with a text file containing the reaction rules and coreactants or this may
    be done on an ad hoc basis."""
    def __init__(self, rule_list=None, coreactant_list=None, explicit_h=True,
                 kekulize=True, neutralise=True, errors=True,
                 racemize=False, database=None, database_overwrite=False,
                 con_string='mongodb://localhost:27017',
                 image_dir=None, quiet=False):
        """
        :param rule_list: Path to a list of reaction rules in TSV form
        :type rule_list: str
        :param coreactant_list: Path to list of coreactants in TSV form
        :type coreactant_list: str
        :param explicit_h: Explicitly represent bound hydrogen atoms
        :type explicit_h: bool
        :param kekulize: Kekulize structures before applying reaction rules
        :type kekulize: bool
        :param neutralise: Remove charges on structure before applying reaction
            rules
        :type neutralise: bool
        :param errors: Print underlying RDKit warnings and halt on error
        :type errors: bool
        :param racemize: Enumerate all possible chiral forms of a molecule if
            unspecified stereocenters exist
        :type racemize: bool
        :param database: Name of desired Mongo Database
        :type database: str
        :param image_dir: Path to desired image folder
        :type image_dir: str
        :param quiet: Silence unbalenced reaction warnings
        :type quiet: bool
        """
        self.operators = {}
        self.coreactants = {}
        self._raw_compounds = {}
        self.compounds = {}
        self.reactions = {}
        self.generation = 0
        self.explicit_h = explicit_h
        self.kekulize = kekulize
        self.racemize = racemize
        self.neutralise = neutralise
        self.image_dir = image_dir
        self.errors = errors
        self.quiet = quiet
        self.fragmented_mols = False
        self.radical_check = False
        self.structure_field = None
        self.con_string = con_string
        # For tanimoto filtering
        self.target_fps = []
        self.crit_tani = None
        self.tani_filter = False

        # If a database is specified make sure it does not exist already
        if database:
            db = MINE(database, con_string)
            if database in db.client.list_database_names():
                if database_overwrite:
                    print(f"Database {database} already exists. Deleting database and removing from core compound mines.")
                    db.core_compounds.update_many({}, {'$pull' : {'MINES' : database}})
                    db.client.drop_database(database)
                    self.mine = database
                else:
                    print(f"Warning! Database {database} already exists. Quitting to avoid overwriting")
                    exit("Exiting due to database name collision.")
                    self.mine = None
            else:
                self.mine = database
            del(db)

        else:
            self.mine = None        

        # Use RDLogger to catch errors in log file. SetLevel indicates mode (
        # 0 - debug, 1 - info, 2 - warning, 3 - critical). Default is no errors
        logger = RDLogger.logger()
        if not errors:
            logger.setLevel(4)

        # Load coreactants (if any) into Pickaxe object
        if coreactant_list:
            with open(coreactant_list) as infile:
                for coreactant in infile:
                    self._load_coreactant(coreactant)

        # Load rules (if any) into Pickaxe object
        if rule_list:
            self._load_operators(rule_list)    

    def load_target_compounds(self, target_compound_file=None, crit_tani=0, 
                                structure_field=None, id_field='id'):
        """
        Loads the target list into an list of fingerprints to later compare to compounds to determine 
        if those compounds should be expanded.

        :param target_compound_file: Path to a file containing compounds as tsv
        :type target_compound_file: basestring
        :param crit_tani: The critical tanimoto cutoff for expansion
        :type crit_tani: float
        :param structure_field: the name of the column containing the
            structure incarnation as Inchi or SMILES (Default:'structure')
        :type structure_field: str
        :param id_field: the name of the column containing the desired
            compound ID (Default: 'id)
        :type id_field: str
        :return: compound SMILES
        :rtype: list
        """
        # Specify options for tanimoto filtering
        self.tani_filter = True
        self.crit_tani = crit_tani

        # Set structure field to None otherwise load_compounds
        # could interferece and vice versa
        self.structure_field = None

        # Load target compounds
        target_smiles = []
        if target_compound_file:
            for line in utils.file_to_dict_list(target_compound_file):
                mol = self._mol_from_dict(line, structure_field)
                if not mol:
                    continue
                # Add compound to internal dictionary as a target
                # compound and store SMILES string to be returned
                smi = AllChem.MolToSmiles(mol, True)
                _id = line[id_field]
                # Do not operate on inorganic compounds
                if 'C' in smi or 'c' in smi:
                    AllChem.SanitizeMol(mol)
                    self._add_compound(_id, smi, mol=mol,
                                       cpd_type='Target Compound')
                    
                    target_smiles.append(smi)

                    # Generate fingerprints for tanimoto filtering
                    fp = AllChem.RDKFingerprint(mol)
                    self.target_fps.append(fp)
        
        else:
            raise ValueError("No input file or database specified for "
                             "target compounds")
        print(f"{len(target_smiles)} target compounds loaded")
        return target_smiles    

    

    def load_compound_set(self, compound_file=None, database=None, con_string='mongodb://localhost:27017',
                         structure_field=None, id_field='id'):
            """If a compound file is provided, this function loads the compounds
            into it's internal dictionary. If not, it attempts to find the
            compounds in its associated MINE database. 
            Compound_file supercedes database. Run once for each source.

            :param compound_file: Path to a file containing compounds as tsv
            :type compound_file: basestring
            :param database: Existing MINE to load compounds from
            :type database: basestring
            :param structure_field: the name of the column containing the
                structure incarnation as Inchi or SMILES (Default:'structure')
            :type structure_field: str
            :param id_field: the name of the column containing the desired
                compound ID (Default: 'id)
            :type id_field: str
            :return: compound SMILES
            :rtype: list
            """
            self.structure_field = None
            compound_smiles = []
            if compound_file:
                for line in utils.file_to_dict_list(compound_file):
                    mol = self._mol_from_dict(line, structure_field)
                    if not mol:
                        continue
                    # Add compound to internal dictionary as a starting
                    # compound and store SMILES string to be returned
                    smi = AllChem.MolToSmiles(mol, True)
                    _id = line[id_field]
                    # Do not operate on inorganic compounds
                    if 'C' in smi or 'c' in smi:
                        AllChem.SanitizeMol(mol)
                        self._add_compound(_id, smi, mol=mol,
                                        cpd_type='Starting Compound')
                        compound_smiles.append(smi)
            # If a MINE database is being used instead, search for compounds
            # annotated as starting compounds and return those as a list of
            # SMILES strings
            elif database:
                db = MINE(database, con_string)
                # Check to see if database has any compounds
                if db.compounds.find_one({'Type':'Starting Compound'}):
                    for compound in db.compounds.find({'Type':'Starting Compound'}):
                        _id = compound['_id']
                        smi = compound['SMILES']
                        # Assume unannotated compounds are starting compounds
                        if 'type' not in compound:
                            compound['Type'] = 'Starting Compound'
                        self._add_compound(_id, smi, cpd_type=compound['Type'])
                        compound_smiles.append(smi)
                else:
                    raise ValueError('Specified MINE contains no starting compounds.')
            else:
                raise ValueError('No input file or database specified for '
                                'starting compounds')
            print(f"{len(compound_smiles)} compounds loaded")
            return compound_smiles

    def transform_compound(self, compound_smiles, rules=None):
        """Perform transformations to a compound returning the products and the
        predicted reactions

        :param compound_smiles: The compound on which to operate represented
            as SMILES
        :type compound_smiles: string
        :param rules: The names of the reaction rules to apply. If none,
            all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as tuple of compound and reaction dicts
        :rtype: tuple
        """
        if not rules:
            rules = self.operators.keys()
        # Create Mol object from input SMILES string and remove hydrogens
        # (rdkit)
        mol = AllChem.MolFromSmiles(compound_smiles)
        mol = AllChem.RemoveHs(mol)
        if not mol:
            if self.errors:
                raise ValueError(f"Unable to parse: {compound_smiles}")
            else:
                print(f"Unable to parse: {compound_smiles}")
                return
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        for rule_name in rules:
            # Lookup rule in dictionary from rule name
            rule = self.operators[rule_name]
            # Get RDKit Mol objects for reactants
            reactant_mols = tuple([mol if x == 'Any'
                                   else self.coreactants[x][0]
                                   for x in rule[1]['Reactants']])
            # Perform chemical reaction on reactants for each rule
            try:
                product_sets = rule[0].RunReactants(reactant_mols)
            # This error should be addressed in a new version of RDKit
            except RuntimeError:
                print("Runtime ERROR!" + rule_name)
                print(compound_smiles)
                continue
            # No enumeration for reactant stereoisomers. _make_half_rxn
            # returns a generator, the first element of which is the reactants.
            reactants, reactant_atoms = next(
                self._make_half_rxn(reactant_mols, rule[1]['Reactants']))
            # By sorting the reactant (and later products) we ensure that
            # compound order is fixed.
            reactants.sort()
            # TODO: check for balanced rxn and remove
            for product_mols in product_sets:
                try:
                    for stereo_prods, product_atoms in self._make_half_rxn(
                            product_mols, rule[1]['Products'], self.racemize):
                        # Update predicted compounds list
                        stereo_prods.sort()
                        # Get reaction text (e.g. A + B <==> C + D)
                        text_rxn = self._add_reaction(reactants, rule_name,
                                                      stereo_prods)

                        # text_rxn is None if reaction has already been inserted
                        if text_rxn:
                            if not self.quiet:
                                self._check_atom_balance(product_atoms,
                                                        reactant_atoms, rule_name,
                                                        text_rxn)

                                # check this and remove unbalanced reactions
                except (ValueError, MemoryError) as e:
                    # TODO: silence
                    if not self.quiet:
                        print(e)
                        print("Error Processing Rule: " + rule_name)
                    continue
        return self.compounds, self.reactions        
    
    def transform_all(self, num_workers=1, max_generations=1):
        """This function applies all of the reaction rules to all the compounds
        until the generation cap is reached.

        :param num_workers: The number of CPUs to for the expansion process.
        :type num_workers: int
        :param max_generations: The maximum number of times an reaction rule
            may be applied
        :type max_generations: int
        """
        def print_progress(done, total):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.05 * total), 1)
            if not done % print_on:
                print(f"Generation {self.generation}: {round(done / total * 100)} percent complete")

        while self.generation < max_generations:
            if self.tani_filter == True:
                if not self.target_fps:
                    print(f'No targets to filter for. Can\'t expand.')
                    return None                

                # Flag compounds to be expanded
                self._filter_by_tani(num_workers=num_workers)

            self.generation += 1
            # Use to print out time per generation at end of loop
            time_init = time.time()
            n_comps = len(self.compounds)
            n_rxns = len(self.reactions)
            # Get all SMILES strings for compounds
            compound_smiles = [c['SMILES'] for c in self.compounds.values()
                            if c['Generation'] == self.generation - 1
                            and c['Type'] not in ['Coreactant', 'Target Compound']
                            and c['Expand'] == True]

            if not compound_smiles:
                continue
            if num_workers > 1:
                chunk_size = max(
                    [round(len(compound_smiles) / (num_workers * 10)), 1])
                print("Chunk Size:", chunk_size)
                new_comps = deepcopy(self.compounds)
                new_rxns = deepcopy(self.reactions)
                pool = multiprocessing.Pool(processes=num_workers)
                for i, res in enumerate(pool.imap_unordered(
                        self.transform_compound, compound_smiles, chunk_size)):
                    new_comps.update(res[0])
                    new_rxns.update(res[1])
                    print_progress((i + 1), len(compound_smiles))
                self.compounds = new_comps
                self.reactions = new_rxns

            else:
                for i, smi in enumerate(compound_smiles):
                    # Perform possible reactions on compound
                    try:
                        self.transform_compound(smi)
                        print_progress(i, len(compound_smiles))
                    except:
                        # TODO what error
                        continue

            print(f"Generation {self.generation} took {time.time()-time_init} sec and produced:")
            print(f"\t\t{len(self.compounds) - n_comps} new compounds")
            print(f"\t\t{len(self.reactions) - n_rxns} new reactions")
            print(f'----------------------------------------\n')


    def _load_coreactant(self, coreactant_text):
        """
        Loads a coreactant into the coreactant dictionary from a tab-delimited
            string
        :param coreactant_text: tab-delimited string with the compound name and
            SMILES
        """
        # If coreactant is commented out (with '#') then don't import
        if coreactant_text[0] == '#':
            return
        split_text = coreactant_text.strip().split('\t')
        # split_text[0] is compound name, split_text[1] is SMILES string
        # Generate a Mol object from the SMILES string if possible
        try:
            mol = AllChem.MolFromSmiles(split_text[2])
            if not mol:
                raise ValueError
            # Generate SMILES string with stereochemistry taken into account
            smi = AllChem.MolToSmiles(mol, True)
        except (IndexError, ValueError):
            raise ValueError(f"Unable to load coreactant: {coreactant_text}")
        _id = self._add_compound(split_text[0], smi, mol=mol,
                                 cpd_type='Coreactant')
        # If hydrogens are to be explicitly represented, add them to the Mol
        # object
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        # If kekulization is preferred (no aromatic bonds, just 3 C=C bonds
        # in a 6-membered aromatic ring for example)
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        # Store coreactant in a coreactants dictionary with the Mol object
        # and hashed id as values (coreactant name as key)
        self.coreactants[split_text[0]] = (mol, _id,)

    def _load_operators(self, rule_path):
        """Loads all reaction rules from file_path into rxn_rule dict.

        :param rule_path: path to file
        :type rule_path: str
        """
        skipped = 0
        with open(rule_path) as infile:
            # Get all reaction rules from tsv file and store in dict (rdr)
            rdr = csv.DictReader((row for row in infile if not
                                  row.startswith('#')), delimiter='\t')
            for rule in rdr:
                try:
                    # Get reactants and products for each reaction into list
                    # form (not ; delimited string)
                    rule['Reactants'] = rule['Reactants'].split(';')
                    rule['Products'] = rule['Products'].split(';')
                    # Ensure that all coreactants are known and accounted for
                    all_rules = rule['Reactants'] + rule['Products']
                    for coreactant_name in all_rules:
                        if ((coreactant_name not in self.coreactants
                             and coreactant_name != 'Any')):
                            raise ValueError(f"Undefined coreactant:{coreactant_name}")
                    # Create ChemicalReaction object from SMARTS string
                    rxn = AllChem.ReactionFromSmarts(rule['SMARTS'])
                    rule.update({'_id': rule['Name'],
                                 'Reactions_predicted': 0,
                                 'SMARTS': rule['SMARTS']})
                    # Ensure that we have number of expected reactants for
                    # each rule
                    if rxn.GetNumReactantTemplates() != len(rule['Reactants'])\
                            or rxn.GetNumProductTemplates() != \
                            len(rule['Products']):
                        skipped += 1
                        print("The number of coreactants does not match the "
                              "number of compounds in the SMARTS for reaction "
                              "rule: " + rule['Name'])
                    if rule['Name'] in self.operators:
                        raise ValueError("Duplicate reaction rule name")
                    # Update reaction rules dictionary
                    self.operators[rule['Name']] = (rxn, rule)
                except Exception as e:
                    raise ValueError(str(e) + f"\nFailed to parse {rule['Name']}")
        if skipped:
            print("WARNING: {skipped} rules skipped")  

    def _filter_by_tani(self, num_workers=1):
        """ 
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a tanimoto similarity score to a target compound
        greater than or equal to the crit_tani, for expansion.
        """
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.05 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total * 100)} percent complete")

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = [cpd for cpd in self.compounds.values() 
                                if cpd['Generation'] == self.generation
                                and cpd['Type'] not in ['Coreactant', 'Target Compound']]
        
        # Set up parallel computing of compounds to expand
        chunk_size = max(
                    [round(len(compounds_to_check) / (num_workers * 10)), 1])
        print(f'Filtering Generation {self.generation}')
        pool = multiprocessing.Pool(num_workers)
        for i, res in enumerate(pool.imap_unordered(
                self._compare_to_targets, 
                [cpd for cpd in compounds_to_check 
                    if cpd['Generation'] == self.generation], chunk_size)):
            
            # If the result of comparison is false, compound is not expanded
            # Default value for a compound is True, so no need to specify expansion
            if not res[1]:
                self.compounds[res[0]]['Expand'] = False            

            print_progress(i, len(compounds_to_check), 'Tanimoto filter progress:')
        return None
    
    def _compare_to_targets(self, cpd):
        """ 
        Helper function to allow parallel computation of tanimoto filtering.
        Works with _filter_by_tani

        Returns True if a the compound is similar enough to a target.
        """
        # Generate the fingerprint of a compound and compare to the fingerprints of the targets
        try:
            fp1 = get_fp(cpd['SMILES'])
            for fp2 in self.target_fps:
                if AllChem.DataStructs.FingerprintSimilarity(fp1, fp2) >= self.crit_tani:
                    return (cpd['_id'], True)
        except:
            pass

        return (cpd['_id'], False)

    def _mol_from_dict(self, input_dict, structure_field=None):
        # detect structure field as needed
        if not structure_field:
            if not self.structure_field:
                for field in input_dict:
                    if str(field).lower() in {'smiles', 'inchi', 'structure'}:
                        self.structure_field = field
                        break
            if not self.structure_field:
                raise ValueError('Structure field not found in input')
            structure_field = self.structure_field

        if structure_field not in input_dict:
            return
        # Generate Mol object from InChI code if present
        if 'InChI=' in input_dict[structure_field]:
            mol = AllChem.MolFromInchi(input_dict[structure_field])
        # Otherwise generate Mol object from SMILES string
        else:
            mol = AllChem.MolFromSmiles(input_dict[structure_field])
        if not mol:
            if self.errors:
                print(f"Unable to Parse {input_dict[structure_field]}")
            return
        # If compound is disconnected (determined by GetMolFrags
        # from rdkit) and loading of these molecules is not
        # allowed (i.e. fragmented_mols == 1), then don't add to
        # internal dictionary. This is most common when compounds
        # are salts.
        if not self.fragmented_mols and len(AllChem.GetMolFrags(mol)) > 1:
            return
        # If specified remove charges (before applying reaction
        # rules later on)
        if self.neutralise:
            mol = utils.neutralise_charges(mol)
        return mol

    def _add_compound(self, cpd_id, smi, mol=None, cpd_type='Predicted'):
        """Adds a compound to the internal compound dictionary"""
        _id = utils.compound_hash(smi, cpd_type)
        self._raw_compounds[smi] = _id
        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        if _id not in self.compounds:
            if not mol:
                mol = AllChem.MolFromSmiles(smi)
            i_key = AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))
            # expand only Predicted and Starting_compounds
            expand = True if cpd_type in ['Predicted', 'Starting Compound'] else False
            self.compounds[_id] = {'ID': cpd_id, '_id': _id, 'SMILES': smi,
                                   'Inchi': AllChem.MolToInchi(mol),
                                   'Inchikey': i_key, 
                                   'Type': cpd_type,
                                   'Generation': self.generation,
                                   'Formula': AllChem.CalcMolFormula(mol),
                                   '_atom_count': self._get_atom_count(mol),
                                #    'Charge': AllChem.GetFormalCharge(mol),
                                   'Reactant_in': [], 'Product_of': [],
                                #    'Sources': [],
                                   'Expand': expand}
            if _id[0] =='X':
                del(self.compounds[_id]['Reactant_in'])
                del(self.compounds[_id]['Product_of'])
            # KMS Uncomment if using sources
            ## Don't track sources of coreactants
            # if _id[0] == 'X':
            #     del self.compounds[_id]['Sources']
            # If we are building a mine and generating images, do so here
            if self.image_dir and self.mine:
                try:
                    with open(os.path.join(self.image_dir, _id + '.svg'),
                              'w') as outfile:
                        nmol = rdMolDraw2D.PrepareMolForDrawing(mol)
                        d2d = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
                        d2d.DrawMolecule(nmol)
                        d2d.FinishDrawing()
                        outfile.write(d2d.GetDrawingText())
                except OSError:
                    print(f"Unable to generate image for {smi}")
        return _id

    def _add_reaction(self, reactants, rule_name, stereo_prods):
        """Hashes and inserts reaction into reaction dictionary"""
        # Hash reaction text
        rhash = rxn2hash(reactants, stereo_prods)
        # Generate unique hash from InChI keys of reactants and products
        # inchi_rxn_hash, text_rxn = \
        #     self._calculate_rxn_hash_and_text(reactants, stereo_prods)
        text_rxn = self._calculate_rxn_text(reactants, stereo_prods)
        # Add reaction to reactions dictionary if not already there
        if rhash not in self.reactions:
            self.reactions[rhash] = {'_id': rhash,
                                     'Reactants': reactants,
                                     'Products': stereo_prods,
                                    #  'InChI_hash': inchi_rxn_hash,
                                     'Operators': {rule_name},                                    
                                     'SMILES_rxn': text_rxn,
                                     'Generation': self.generation}

        # Otherwise, update the operators
        else:
            self.reactions[rhash]['Operators'].add(rule_name)

        # Update compound tracking
        for prod_id in [x.c_id for x in stereo_prods if x.c_id[0] == 'C']:
            if rhash not in self.compounds[prod_id]['Product_of']:
                self.compounds[prod_id]['Product_of'].append(rhash)
        
        for reac_id in [x.c_id for x in reactants if x.c_id[0] == 'C']:
            if rhash not in self.compounds[reac_id]['Reactant_in']:
                self.compounds[reac_id]['Reactant_in'].append(rhash)      

        return text_rxn

    def _check_atom_balance(self, product_atoms, reactant_atoms, rule_name,
                            text_rxn):
        """If the SMARTS rule is not atom balanced, this check detects the
        accidental alchemy."""
        if reactant_atoms - product_atoms \
                or product_atoms - reactant_atoms:
            return False
            if self.quiet == False:
                print("Warning: Unbalanced Reaction produced by "
                    + rule_name)
                print(text_rxn)
                print(reactant_atoms, product_atoms)
        else:
            return True
            

    def _get_atom_count(self, mol):
        """Takes a set of mol objects and returns a counter with each element
        type in the set"""
        atoms = collections.Counter()
        # Find all strings of the form A# in the molecular formula where A
        # is the element (e.g. C) and # is the number of atoms of that
        # element in the molecule. Pair is of form [A, #]
        for pair in re.findall(r'([A-Z][a-z]*)(\d*)',
                               AllChem.CalcMolFormula(mol)):
            # Add # to atom count, unless there is no # (in which case
            # there is just one of that element, as ones are implicit in
            # chemical formulas)
            if pair[1]:
                atoms[pair[0]] += int(pair[1])
            else:
                atoms[pair[0]] += 1
        if self.radical_check:
            radical = any([atom.GetNumRadicalElectrons()
                           for atom in mol.GetAtoms()])
            if radical:
                atoms['*'] += 1
        return atoms    

    def _make_half_rxn(self, mol_list, rules, split_stereoisomers=False):
        """Takes a list of mol objects for a half reaction, combines like
        compounds and returns a generator for stoich tuples"""
        # Get compound ids from Mol objects, except for coreactants, in which
        #  case we look them up in the coreactant dictionary
        comps = [self._calculate_compound_information(m, split_stereoisomers)
                 if r == 'Any' else (self.coreactants[r][1],)
                 for m, r in zip(mol_list, rules)]
        # count the number of atoms on a side
        atom_count = collections.Counter()
        for x in comps:
            atom_count += self.compounds[x[0]]['_atom_count']
        # Remove duplicates from comps
        half_rxns = [collections.Counter(subrxn) for subrxn
                     in itertools.product(*comps)]
        # Yield generator with stoichiometry tuples (of form stoichiometry, id)
        for rxn in half_rxns:
            yield [StoichTuple(y, x) for x, y in rxn.items()], atom_count

    def _calculate_compound_information(self, mol_obj, split_stereoisomers):
        """Calculate the standard data for a compound & return a tuple with
        compound_ids. Memoized with _raw_compound dict"""
        # This is the raw SMILES which may have explicit hydrogen
        raw = AllChem.MolToSmiles(mol_obj, True)
        if raw not in self._raw_compounds:
            try:
                # Remove hydrogens if explicit (this step slows down the
                # process quite a bit)
                if self.explicit_h:
                    mol_obj = AllChem.RemoveHs(mol_obj)
                AllChem.SanitizeMol(mol_obj)
            except ValueError:
                if not self.quiet:
                    print("Product ERROR!: " + raw)
                raise ValueError
            # In case we want to have separate entries for stereoisomers
            if split_stereoisomers:
                processed_mols = _racemization(mol_obj)
            else:
                processed_mols = [mol_obj]
            # Get list of SMILES string(s) from Mol object(s)
            smiles = [AllChem.MolToSmiles(m, True) for m in processed_mols]
            # Add compound(s) to internal dictionary
            self._raw_compounds[raw] = tuple([self._add_compound(None, s,
                                                                 mol=m)
                                              for s, m in zip(smiles,
                                                              processed_mols)])
        return self._raw_compounds[raw] if isinstance(
            self._raw_compounds[raw], tuple) else (self._raw_compounds[raw],)

    def _calculate_rxn_text(self, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is
        connectivity only, second block is stereo only"""
        def get_blocks(tups):
            smiles = []
            for x in tups:
                comp = self.compounds[x.c_id]
                smiles.append(f"({x.stoich}) {comp['SMILES']}")
                
            return ' + '.join(smiles)

        r_s = get_blocks(reactants)
        p_s = get_blocks(products)
        smiles_rxn = r_s + ' => ' + p_s
        
        return smiles_rxn

    def _calculate_rxn_hash_and_text(self, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is
        connectivity only, second block is stereo only"""
        def __get_blocks(tups):
            first_block, second_block, smiles = [], [], []
            for x in tups:
                comp = self.compounds[x.c_id]
                smiles.append(f"({x.stoich}) {comp['SMILES']}")
                if comp['Inchikey']:
                    # InChI keys are separated by a hyphen, where the first
                    # part is derived from connectivity and the second part
                    # comes from other layers such as stereochemistry
                    split_inchikey = comp['Inchikey'].split('-')
                    if len(split_inchikey) > 1:
                        first_block.append(f'{x.stoich},{split_inchikey[0]}')
                        second_block.append(f'{x.stoich},{split_inchikey[1]}')
                else:
                    print(f"No Inchikey for {x.c_id}")
            return '+'.join(first_block), '+'.join(second_block), \
                   ' + '.join(smiles)

        reactants.sort()
        products.sort()
        r_1, r_2, r_s = __get_blocks(reactants)
        p_1, p_2, p_s = __get_blocks(products)
        first_block = r_1 + '<==>' + p_1
        second_block = r_2 + '<==>' + p_2
        smiles_rxn = r_s + ' => ' + p_s

        return hashlib.sha256(first_block.encode()).hexdigest() + '-' \
            + hashlib.md5(second_block.encode()).hexdigest(), smiles_rxn

    def assign_ids(self):
        """Assigns a numerical ID to compounds (and reactions) for ease of
        reference. Unique only to the CURRENT run."""
        # If we were running a multiprocess expansion, this removes the dicts
        # from Manager control
        self.compounds = dict(self.compounds)
        self.reactions = dict(self.reactions)
        i = 1
        for comp in sorted(self.compounds.values(),
                           key=lambda x: (x['Generation'], x['_id'])):
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            if not comp['ID']:
                comp['ID'] = 'pkc' + str(i).zfill(7)
                i += 1
                self.compounds[comp['_id']] = comp
                # If we are not loading into the mine, we generate the image
                # here.
                if self.image_dir and not self.mine:
                    mol = AllChem.MolFromSmiles(comp['SMILES'])
                    try:
                        MolToFile(
                            mol,
                            os.path.join(self.image_dir, comp['ID'] + '.png'),
                            fitImage=True, kekulize=False)
                    except OSError:
                        print(f"Unable to generate image for {comp['SMILES']}")
        i = 1
        for rxn in sorted(self.reactions.values(),
                          key=lambda x: (x['Generation'], x['_id'])):
            rxn['ID_rxn'] = ' + '.join(
                [f"({x.stoich}) {self.compounds[x.c_id]['ID']}[c0]"
                 for x in rxn['Reactants']]) + ' => ' + ' + '.join(
                     [f"({x.stoich}) {self.compounds[x.c_id]['ID']}[c0]"
                      for x in rxn['Products']])
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            rxn['ID'] = 'pkr' + str(i).zfill(7)
            i += 1
            self.reactions[rxn['_id']] = rxn

    

    def prune_network(self, white_list):
        """
        Prune the predicted reaction network to only compounds and reactions
        that terminate in a specified white list of compounds.
        :param white_list: A list of compound_ids to include (if found)
        :type white_list: list
        :return: None
        """
        n_white = len(white_list)
        comp_set, rxn_set = self.find_minimal_set(white_list)
        print(f"Pruned network to {len(comp_set)} compounds and {len(rxn_set)} reactions based on \
                {n_white} whitelisted compounds")
        self.compounds = dict([(k, v) for k, v in self.compounds.items()
                               if k in comp_set])
        self.reactions = dict([(k, v) for k, v in self.reactions.items()
                               if k in rxn_set])

    def find_minimal_set(self, white_list):
        """
        Given a whitelist this function finds the minimal set of compound and
        reactions ids that comprise the set
        :param white_list:  A list of compound_ids to include (if found)
        :type white_list: list
        :return: compound and reaction id sets
        :rtype: tuple(set, set)
        """
        white_set = set(white_list)
        comp_set = set()
        rxn_set = set()
        for c_id in white_list:
            if c_id not in self.compounds:
                continue
            for r_id in self.compounds[c_id]['Product_of']:
                rxn_set.add(r_id)
                comp_set.update([x.c_id for x
                                 in self.reactions[r_id]['Products']])
                for tup in self.reactions[r_id]['Reactants']:
                    comp_set.add(tup.c_id)
                    # do not want duplicates or cofactors in the whitelist
                    if tup.c_id[0] == 'C' and tup.c_id not in white_set:
                        white_list.append(tup.c_id)
                        white_set.add(tup.c_id)
        return comp_set, rxn_set

    def write_compound_output_file(self, path, dialect='excel-tab'):
        """Writes all compound data to the specified path.

        :param path: path to output
        :type path: str
        :param dialect: the output format for the file. Choose excel for csv
            excel-tab for tsv.
        :type dialect: str
        """
        path = utils.prevent_overwrite(path)

        columns = ('ID', 'Type', 'Generation', 'Formula', 'Inchikey',
                   'SMILES')
        with open(path, 'w') as outfile:
            writer = csv.DictWriter(outfile, columns, dialect=dialect,
                                    extrasaction='ignore', lineterminator='\n')
            writer.writeheader()
            writer.writerows(sorted(self.compounds.values(),
                                    key=lambda x: x['ID']))

    def write_reaction_output_file(self, path, delimiter='\t'):
        """Writes all reaction data to the specified path.

        :param path: path to output
        :type path: basestring
        :param delimiter: the character with which to separate data entries
        :type delimiter: basestring
        """
        path = utils.prevent_overwrite(path)
        with open(path, 'w') as outfile:
            outfile.write('ID\tName\tID equation\tSMILES equation\tRxn hash\t'
                          'Reaction rules\n')
            for rxn in sorted(self.reactions.values(), key=lambda x: x['ID']):
                outfile.write(delimiter.join([rxn['ID'], '', rxn['ID_rxn'],
                                              rxn['SMILES_rxn'], rxn['_id'],
                                              ';'.join(rxn['Operators'])])
                              + '\n')

    def _save_compound_helper(self, comp_dict):
        # Helper function to aid parallelization of saving compounds in
        # save_to_mine
        if not comp_dict['_id'].startswith('T'):
            # These functions are outside of the MINE class in order to
            # allow for parallelization. When in the MINE class it is not
            # serializable with pickle. In comparison to the class functions,
            # these return the requests instead of appending to a passed list.
            mine_req = databases.insert_mine_compound(comp_dict)
            core_up_req = databases.update_core_compound_MINES(comp_dict, self.mine)
            core_in_req = databases.insert_core_compound(comp_dict)
            return [mine_req, core_up_req, core_in_req]
        else:
            return None

    def _save_target_helper(self, comp_dict):
        # Helper function to aid parallelization of saving targets in
        # save_to_mine
        if comp_dict['_id'].startswith('T'):
            # This functions are outside of the MINE class in order to
            # allow for parallelization. When in the MINE class it is not
            # serializable with pickle. In comparison to the class functions,
            # these return the requests instead of appending to a passed list.
            target_req = databases.insert_mine_compound(comp_dict)
            return target_req
        else:
            return None

    def save_to_mine(self, num_workers, indexing=True):
        """Save compounds to a MINE database.

        :param db_id: The name of the target database
        :type db_id: basestring
        """        
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.05 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total * 100)} percent complete")
            
        # Initialize variables for saving to mongodb
        print(f'----------------------------------------')
        print(f'Saving results to {self.mine}')
        print(f'----------------------------------------\n')
        start = time.time()
        db = MINE(self.mine, self.con_string)
        # Lists of requests for bulk_write method
        core_cpd_requests = []
        core_update_mine_requests = []
        mine_cpd_requests = []
        mine_rxn_requests = []        
        target_cpd_requests = []

        
        print('--------------- Reactions ---------------')
        # Insert Reactions
        rxn_start = time.time()
        if num_workers > 1:
            # parallel insertions
            chunk_size = max(
                [round(len(self.reactions) / (num_workers * 10)), 1])
            print("Chunk Size for writing Reactions:", chunk_size)
            pool = multiprocessing.Pool(processes=num_workers)
            n = 0
            # Call insert_reaction in parallel and append results to the mine_rxn_requests
            for i, res in enumerate(pool.imap_unordered(
                    databases.insert_reaction, self.reactions.values(), chunk_size)):
                if res:
                    mine_rxn_requests.append(res)
                    if i % (round(chunk_size)) == 0:
                        print_progress(i, len(self.reactions), 'Reaction')
        else:
            # non-parallel insertions
            for rxn in self.reactions.values():
                db.insert_reaction(rxn, requests=mine_rxn_requests)
                for op in rxn['Operators']:
                    self.operators[op][1]['Reactions_predicted'] += 1

        if mine_rxn_requests:            
            db.reactions.bulk_write(mine_rxn_requests, ordered=False)
            print(f'Done with reactions--took {time.time() - rxn_start} seconds.')
            db.meta_data.insert_one({'Timestamp': datetime.datetime.now(),
                                        'Action': "Reactions Inserted"})
            del(mine_rxn_requests)
        else:
            print('No reactions inserted')
        print(f'----------------------------------------\n')

        print('--------------- Compounds --------------')
        # Insert Compoundsd
        cpd_start = time.time()
        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(self.compounds) / (num_workers * 10)), 1])
            print("Chunk Size Compound Writing:", chunk_size)

            pool = multiprocessing.Pool(processes=num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    self._save_compound_helper, self.compounds.values(), chunk_size)):
                if res:
                    mine_cpd_requests.append(res[0])
                    core_update_mine_requests.append(res[1])
                    core_cpd_requests.append(res[2])    
                    print_progress(i, len(self.compounds), 'Compounds')
        else:
            # non-parallel insertion
            # Write generated compounds to MINE and core compounds to core
            for comp_dict in self.compounds.values():
                # Write everything except for targets
                if not comp_dict['_id'].startswith('T'):
                    # These functions are in the MINE class. The request list is
                    # passed and appended in the MINE method.
                    db.insert_mine_compound(comp_dict, mine_cpd_requests)
                    db.update_core_compound_MINES(comp_dict, core_update_mine_requests)
                    db.insert_core_compound(comp_dict, core_cpd_requests)

        print(f'Done with Compounds Prep--took {time.time() - cpd_start} seconds.')
        # Insert the three types of compounds
        cpd_start = time.time()
        db.core_compounds.bulk_write(core_cpd_requests, ordered=False)
        del(core_cpd_requests)
        print(f'Done with Core Compounds Insertion--took {time.time() - cpd_start} seconds.')
        cpd_start = time.time()
        db.core_compounds.bulk_write(core_update_mine_requests, ordered=False)
        del(core_update_mine_requests)
        print(f'Done with Updating Core Compounds MINE update--took {time.time() - cpd_start} seconds.')
        db.meta_data.insert_one({'Timestamp': datetime.datetime.now(),
                                'Action': "Core Compounds Inserted"})
        cpd_start = time.time()
        db.compounds.bulk_write(mine_cpd_requests, ordered=False)     
        del(mine_cpd_requests)  
        print(f'Done with MINE Insert--took {time.time() - cpd_start} seconds.')
        db.meta_data.insert_one({'Timestamp': datetime.datetime.now(),
                                'Action': "Compounds Inserted"})
        print(f'----------------------------------------\n')

        # Insert target compounds
        target_start = time.time()
        # Write target compounds to target collection
        # Target compounds are written as mine compounds
        if self.tani_filter:
            print('--------------- Targets ----------------')
            db.meta_data.insert_one({'Tani Threshold' : self.crit_tani})    

            # Insert target compounds
            target_start = time.time()
            # get the number of targets
            # if num_workers > 1:
            #     # parallel insertion
            #     chunk_size = max(
            #         [round(len(self.target_fps) / (num_workers * 10)), 1])
            #     print("Chunk Size for Target Writing:", chunk_size)

            #     pool = multiprocessing.Pool(processes=num_workers)
            #     for i, res in enumerate(pool.imap_unordered(
            #             self._save_target_helper, self.compounds.values(), chunk_size)):
            #         if res:
            #             target_cpd_requests.append(res)  
            #             print_progress(i, len(self.compounds), 'Targets')
            # else:

            # parallel insertion is taking a weirdly long time
            # non-parallel insertion
            for comp_dict in self.compounds.values():
                if comp_dict['_id'].startswith('T'):
                    db.insert_mine_compound(comp_dict, target_cpd_requests)     
            print(f'Done with Target Prep--took {time.time() - target_start} seconds.')    
            target_start = time.time()
            db.target_compounds.bulk_write(target_cpd_requests, ordered=False)
            del(target_cpd_requests)
            print(f'Done with Target Insert--took {time.time() - target_start} seconds.')
            db.meta_data.insert_one({'Timestamp': datetime.datetime.now(),
                                'Action': "Target Compounds Inserted"})
            print(f'----------------------------------------\n')                   

        # Save operators  
        operator_start = time.time()
        if self.operators:
            print('-------------- Operators ---------------')
            db.operators.insert_many([op[1] for op in self.operators.values()])        
            db.meta_data.insert_one({'Timestamp': datetime.datetime.now(),
                                    'Action': "Operators Inserted"})  
            print(f'Done with Operators Overall--took {time.time() - operator_start} seconds.')
        print(f'----------------------------------------\n')

        if indexing:
            print('-------------- Indices ---------------')
            index_start = time.time()
            db.build_indexes()
            print(f'Done with Indices--took {time.time() - index_start} seconds.')
            print(f'----------------------------------------\n')

        print('-------------- Overall ---------------')
        print(f'Finished uploading everything in {time.time() - start} sec')
        print(f'----------------------------------------\n')


def _racemization(compound, max_centers=3, carbon_only=True):
    """Enumerates all possible stereoisomers for unassigned chiral centers.

    :param compound: A compound
    :type compound: rdMol object
    :param max_centers: The maximum number of unspecified stereocenters to
        enumerate. Sterioisomers grow 2^n_centers so this cutoff prevents lag
    :type max_centers: int
    :param carbon_only: Only enumerate unspecified carbon centers. (other
        centers are often not tautomeric artifacts)
    :type carbon_only: bool
    :return: list of stereoisomers
    :rtype: list of rdMol objects
    """
    new_comps = []
    # FindMolChiralCenters (rdkit) finds all chiral centers. We get all
    # unassigned centers (represented by '?' in the second element
    # of the function's return parameters).
    unassigned_centers = [c[0] for c in AllChem.FindMolChiralCenters(
        compound, includeUnassigned=True) if c[1] == '?']
    # Get only unassigned centers that are carbon (atomic number of 6) if
    # indicated
    if carbon_only:
        unassigned_centers = list(
            filter(lambda x: compound.GetAtomWithIdx(x).GetAtomicNum() == 6,
                unassigned_centers))
    # Return original compound if no unassigned centers exist (or if above
    # max specified (to prevent lag))
    if not unassigned_centers or len(unassigned_centers) > max_centers:
        return [compound]
    for seq in itertools.product([1, 0], repeat=len(unassigned_centers)):
        for atomid, clockwise in zip(unassigned_centers, seq):
            # Get both cw and ccw chiral centers for each center. Used
            # itertools.product to get all combinations.
            if clockwise:
                compound.GetAtomWithIdx(atomid).SetChiralTag(
                    AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            else:
                compound.GetAtomWithIdx(atomid).SetChiralTag(
                    AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        # Duplicate C++ object so that we don't get multiple pointers to
        # same object
        new_comps.append(deepcopy(compound))
    return new_comps


if __name__ == '__main__':
    # Get initial time to calculate execution time at end
    t1 = time.time()  # pylint: disable=invalid-name
    # Parse all command line arguments
    parser = ArgumentParser()  # pylint: disable=invalid-name
    parser.add_argument('-C', '--coreactant_list',
                        default="tests/data/test_coreactants.tsv",
                        help="Specify a list of coreactants as a "
                             "tab-separated file")
    parser.add_argument('-r', '--rule_list',
                        default="tests/data/test_reaction_rules.tsv",
                        help="Specify a list of reaction rules as a "
                             "tab-separated file")
    parser.add_argument('-c', '--compound_file',
                        default="tests/data/test_compounds.tsv",
                        help="Specify a list of starting compounds as a "
                             "tab-separated file")
    parser.add_argument('-p', '--pruning_whitelist', default=None,
                        help="Specify a list of target compounds to prune "
                             "reaction network down")
    parser.add_argument('-s', '--smiles', default=None,
                        help="Specify a starting compound as SMILES.")
    parser.add_argument('-o', '--output_dir', default=".",
                        help="The directory in which to place files")
    parser.add_argument('-d', '--database', default=None,
                        help="The name of the database in which to store "
                             "output. If not specified, data is still written "
                             "as tsv files")
    parser.add_argument('-R', '--racemize', action='store_true', default=False,
                        help="Enumerate the possible chiral "
                        "forms for all unassigned stereocenters in compounds &"
                        " reactions")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help="Display RDKit errors & warnings")
    parser.add_argument('--bnice', action='store_true', default=False,
                        help="Set several options to enable compatibility "
                             "with bnice operators.")
    parser.add_argument('-m', '--max_workers', default=1, type=int,
                        help="Set the nax number of processes to spawn to "
                             "perform calculations.")
    parser.add_argument('-g', '--generations', default=1, type=int,
                        help="Set the numbers of time to apply the reaction "
                             "rules to the compound set.")
    parser.add_argument('-i', '--image_dir', default=None,
                        help="Specify a directory to store images of all "
                             "created compounds")
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help="Silence warnings about imbalenced reactions")
    OPTIONS = parser.parse_args()
    pk = Pickaxe(coreactant_list=OPTIONS.coreactant_list,
                 rule_list=OPTIONS.rule_list, racemize=OPTIONS.racemize,
                 errors=OPTIONS.verbose, explicit_h=OPTIONS.bnice,
                 kekulize=OPTIONS.bnice, neutralise=OPTIONS.bnice,
                 image_dir=OPTIONS.image_dir, database=OPTIONS.database,
                 quiet=OPTIONS.quiet)
    # Create a directory for image output file if it doesn't already exist
    if OPTIONS.image_dir and not os.path.exists(OPTIONS.image_dir):
        os.mkdir(OPTIONS.image_dir)
    # If starting compound specified as SMILES string, then add it
    if OPTIONS.smiles:
        # pylint: disable=protected-access
        pk._add_compound('Start', OPTIONS.smiles, cpd_type='Starting Compound')
    else:
        pk.load_compound_set(compound_file=OPTIONS.compound_file)
    # Generate reaction network
    pk.transform_all(max_generations=OPTIONS.generations,
                     num_workers=OPTIONS.max_workers)
    if OPTIONS.pruning_whitelist:
        # pylint: disable=invalid-name,protected-access
        mols = [pk._mol_from_dict(line) for line
                in utils.file_to_dict_list(OPTIONS.pruning_whitelist)]
        pk.prune_network([utils.compound_hash(x) for x in mols if x])
    # Save to database (e.g. Mongo) if present, otherwise create output file
    if OPTIONS.database:
        print(f"Saving results to {OPTIONS.database}")
        pk.save_to_mine(OPTIONS.database)
    else:
        pk.assign_ids()
        pk.write_compound_output_file(OPTIONS.output_dir + '/compounds.tsv')
        pk.write_reaction_output_file(OPTIONS.output_dir + '/reactions.tsv')

    print(f"Execution took {time.time()-t1} seconds.")


    # TODO: balance
    # TODO: cpd in collection not same for single worker and max worker
    # TODO: fragmented outputs