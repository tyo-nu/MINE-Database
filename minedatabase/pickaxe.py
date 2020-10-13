"""Pickaxe.py: This module generates new compounds from user-specified starting
   compounds using a set of SMARTS-based reaction rules."""
import multiprocessing
import collections
import itertools
import datetime
import hashlib
import time
import csv
import os
import re

from argparse import ArgumentParser
from functools import partial
from copy import deepcopy
from sys import exit

import minedatabase.databases as databases

from minedatabase.databases import MINE
from minedatabase import utils

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Draw import MolToFile, rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit import RDLogger

class Pickaxe:
    """This class generates new compounds from user-specified starting
    compounds using a set of SMARTS-based reaction rules. It may be initialized
    with a text file containing the reaction rules and coreactants or this may
    be done on an ad hoc basis."""
    def __init__(self, rule_list=None, coreactant_list=None, explicit_h=True,
                 kekulize=True, neutralise=True, errors=True,
                 racemize=False, database=None, database_overwrite=False,
                 mongo_uri='mongodb://localhost:27017',
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
        :param database_overwrite: Force overwrite of existing DB
        :type database_overwrite: bool
        :param mongo_uri: URI of mongo deployment
        :type mongo_uri: str
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
        # For tanimoto filtering
        self.target_fps = []
        self.target_smiles = []
        self.crit_tani = None
        self.tani_filter = False
        # database info
        self.mongo_uri = mongo_uri
        
        print(f'----------------------------------------')
        print("Intializing pickaxe object") 
        if database:
            # Determine if a specified database is legal
            db = MINE(database, self.mongo_uri)
            if database in db.client.list_database_names():
                if database_overwrite:
                    # Remove db from core compounds and drop db
                    print(f"Database {database} already exists. ",
                            "Deleting database and removing from ",
                            "core compound mines.")
                    db.core_compounds.update_many({},
                        {'$pull' : {'MINES' : database}}
                    )
                    db.client.drop_database(database)
                    self.mine = database
                else:
                    print(f"Warning! Database {database} already exists."
                            "Specify database_overwrite as true to delete",
                            " old database and write new.")
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

        print("\nDone intializing pickaxe object")
        print(f'----------------------------------------\n') 

    def load_target_set(self, target_compound_file=None, crit_tani=0, 
                                structure_field=None, id_field='id'):
        """
        Loads the target list into an list of fingerprints to later compare
        to compounds to determine if those compounds should be expanded.

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
        # Update options for tanimoto filtering
        self.tani_filter = True
        self.crit_tani = crit_tani

        # Set structure field to None otherwise value determined by
        # load_structures can interfere
        self.structure_field = None

        # Load target compounds
        if target_compound_file:
            for target_dict in utils.file_to_dict_list(target_compound_file):
                mol = self._mol_from_dict(target_dict, structure_field)
                if not mol:
                    continue
                # Add compound to internal dictionary as a target
                # compound and store SMILES string to be returned
                smi = AllChem.MolToSmiles(mol, True)
                cpd_name = target_dict[id_field]
                # Only operate on organic compounds 
                if 'c' in smi.lower():
                    AllChem.SanitizeMol(mol)
                    self._add_compound(cpd_name, smi, 'Target Compound', mol)
                    self.target_smiles.append(smi)                    
                    # Generate fingerprints for tanimoto filtering
                    fp = AllChem.RDKFingerprint(mol)
                    self.target_fps.append(fp)
        else:
            raise ValueError("No input file or database specified for "
                             "target compounds")

        print(f"{len(self.target_smiles)} target compounds loaded")

        return self.target_smiles    

    def load_compound_set(self, compound_file=None, id_field='id'): 
            """If a compound file is provided, this function loads the 
            compounds into its internal dictionary. 

            :param compound_file: Path to a file containing compounds as tsv
            :type compound_file: basestring            
            :param id_field: the name of the column containing the desired
                compound ID (Default: 'id)
            :type id_field: str

            :return: compound SMILES
            :rtype: list
            """
            # TODO: support for multiple sources?
            # For example, loading from MINE and KEGG

            # Set structure field to None otherwise value determined by
            # load_targets can interfere
            self.structure_field = None

            # load compounds
            compound_smiles = []
            if compound_file:
                for cpd_dict in utils.file_to_dict_list(compound_file):
                    mol = self._mol_from_dict(cpd_dict, self.structure_field)
                    if not mol:
                        continue
                    # Add compound to internal dictionary as a starting
                    # compound and store SMILES string to be returned
                    smi = AllChem.MolToSmiles(mol, True)
                    cpd_name = cpd_dict[id_field]
                    # Do not operate on inorganic compounds
                    if 'C' in smi or 'c' in smi:
                        AllChem.SanitizeMol(mol)
                        self._add_compound(cpd_name, smi,
                            cpd_type='Starting Compound', mol=mol
                        )
                        compound_smiles.append(smi)

            # TODO: Add Support for MINE    
            else:
                raise ValueError('No input file specified for '
                                'starting compounds')

            print(f"{len(compound_smiles)} compounds loaded")

            return compound_smiles

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
            # TODO: what do do about stereochemistry? Original comment is below
            # but stereochem was taken out (isn't it removed later anyway?)
            # # Generate SMILES string with stereochemistry taken into account
            smi = AllChem.MolToSmiles(mol)
        except (IndexError, ValueError):
            raise ValueError(f"Unable to load coreactant: {coreactant_text}")
        cpd_id = self._add_compound(split_text[0], smi, 'Coreactant', mol)
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
        self.coreactants[split_text[0]] = (mol, cpd_id,)

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
                            raise ValueError(
                                "Undefined coreactant:{coreactant_name}"
                            )
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
                    raise ValueError(str(e) 
                        + f"\nFailed to parse {rule['Name']}")
        if skipped:
            print("WARNING: {skipped} rules skipped")  

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

    def _gen_compound(self, cpd_name, smi, cpd_type, mol=None):
        """Generates a compound"""
        cpd_dict = {}
        cpd_id = utils.compound_hash(smi, cpd_type)
        self._raw_compounds[smi] = cpd_id
        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        if cpd_id not in self.compounds:
            if not mol:
                mol = AllChem.MolFromSmiles(smi)
            # expand only Predicted and Starting_compounds
            if cpd_type in ['Predicted', 'Starting Compound']:
                expand = True
            else:
                expand = False

            cpd_dict = {'ID': cpd_name, '_id': cpd_id, 'SMILES': smi,
                                   'Type': cpd_type,
                                   'Generation': self.generation,
                                   'atom_count': utils._getatom_count(mol, 
                                                           self.radical_check),
                                   'Reactant_in': [], 
                                   'Product_of': [],
                                   'Expand': expand,
                                   'Formula':CalcMolFormula(mol)}
            if cpd_id.startswith('X'):
                del(cpd_dict['Reactant_in'])
                del(cpd_dict['Product_of'])
        else:
            cpd_dict = self.compounds[cpd_id]
        
        return cpd_id, cpd_dict
    
    def _insert_compound(self, cpd_dict):
        """Inserts a compound into the dictionary"""        
        cpd_id = cpd_dict['_id']
        self.compounds[cpd_id] = cpd_dict

        if self.image_dir and self.mine:
                try:
                    with open(os.path.join(self.image_dir, cpd_id + '.svg'),
                            'w') as outfile:
                        mol = AllChem.MolFromSmiles(cpd_dict['SMILES'])
                        nmol = rdMolDraw2D.PrepareMolForDrawing(mol)
                        d2d = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
                        d2d.DrawMolecule(nmol)
                        d2d.FinishDrawing()
                        outfile.write(d2d.GetDrawingText())
                except OSError:
                    print(f"Unable to generate image for {smi}")

    def _add_compound(self, cpd_name, smi, cpd_type, mol=None):
        """Adds a compound to the internal compound dictionary"""
        cpd_id = utils.compound_hash(smi, cpd_type)
        self._raw_compounds[smi] = cpd_id
        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        if cpd_id not in self.compounds:
            _, cpd_dict = self._gen_compound(cpd_name, smi, cpd_type, mol)
            self._insert_compound(cpd_dict)
        
        return cpd_id
  
        """If the SMARTS rule is not atom balanced, this check detects the
        accidental alchemy."""
        if reactant_atoms - product_atoms \
                or product_atoms - reactant_atoms:
            return False
        else:
            return True      

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
                print(f"Generation {self.generation}: {round(done/total*100)}"
                        ," percent complete")

        while self.generation < max_generations:
            print('----------------------------------------')
            print(f'Expanding Generation {self.generation}')

            if self.tani_filter == True:
                # Starting time for tani filtering
                time_tani = time.time() 
                if not self.target_fps:
                    print(f'No targets to filter for. Can\'t expand.')
                    return None

                # Flag compounds to be expanded
                if type(self.crit_tani) == list:
                    crit_tani = self.crit_tani[self.generation]
                else:
                    crit_tani = self.crit_tani

                print("Filtering out compounds with maximum tanimoto match "
                      f"< {crit_tani}")
                self._filter_by_tani(num_workers=num_workers)
                n_filtered = 0
                n_total = 0                
                for cpd_dict in self.compounds.values():
                    if cpd_dict['Generation'] == (self.generation and 
                                                  cpd_dict['_id'].startswith('C')):
                        n_total += 1
                        if cpd_dict['Expand'] == True:
                            n_filtered += 1

                print(f"{n_filtered} of {n_total} compounds remain after "
                       f"filtering generation {self.generation}--took "
                       f"{time.time() - time_tani}s.\n\nExpanding.")

            # Starting time for expansion
            time_init = time.time()
            self.generation += 1

            # Tracking compounds formed        
            n_comps = len(self.compounds)
            n_rxns = len(self.reactions)

            # Get SMILES to be expanded
            compound_smiles = [cpd['SMILES'] for cpd in self.compounds.values()
                            if cpd['Generation'] == self.generation - 1
                            and cpd['Type'] not in ['Coreactant', 'Target Compound']
                            and cpd['Expand'] == True]
            # No compounds found 
            if not compound_smiles:
                print(f"No compounds to expand in generation"
                       f"{self.generation}. Finished expanding.")
                return None

            self._transform_helper(compound_smiles, num_workers)
            
            print(f"Generation {self.generation} took {time.time()-time_init}"
                   "sec and produced:")
            print(f"\t\t{len(self.compounds) - n_comps} new compounds")
            print(f"\t\t{len(self.reactions) - n_rxns} new reactions")
            print(f'----------------------------------------\n')

    def remove_cofactor_redundancy(self):
        """
        Checks for cofactors in rxns that were generated by an any;any rules and have a C.
        First check for cpd_ids that, when replacing the C, match a cofactor id. 
        Then remove or replace reactions that contain one of these compounds.    
        """

        cofactors_as_cpds = []
        cofactor_ids = [cofactor[1] for cofactor in self.coreactants.values()]
        for cpd_id in self.compounds:
            if 'X' + cpd_id[1:] in cofactor_ids and cpd_id.startswith('C'):
                cofactors_as_cpds.append(cpd_id)

        # identify reactions that need to be deleted or updated
        for cpd_id in cofactors_as_cpds:
            rxn_ids = set(self.compounds[cpd_id]['Product_of'] + self.compounds[cpd_id]['Reactant_in'])
            for rxn_id in rxn_ids:
                # cpd cof replacement
                # TODO: fix this
                rxn = self.reactions[rxn_id]
                products = list(map(lambda x: (x[0], self.compounds[x[1]]), rxn['Products']))
                reactants = list(map(lambda x: (x[0], self.compounds[x[1]]), rxn['Reactants']))
                rxn_hash, _ = utils.rxn2hash(reactants, products)
                if rxn_hash in self.reactions.keys():
                    # More specific exists so delete broader version
                    # delete old IDs
                    for _, prod in products:
                        if prod['_id'].startswith('C'):
                            self.compounds[prod['_id']]['Product_of'].remove(rxn_id)
                    
                    for _, reac in reactants:
                        if reac['_id'].startswith('C'):
                            self.compounds[reac['_id']]['Reactant_in'].remove(rxn_id)

                    # delete reaction
                    del(self.reactions[rxn_id])

                else:
                    print('err') # TODO
        
            del(self.compounds[cpd_id])

    def _filter_by_tani(self, num_workers=1):
        """ 
        Compares the current generation to the target compound fingerprints
        marking compounds, who have a tanimoto similarity score to a target compound
        greater than or equal to the crit_tani, for expansion.
        """
        # TODO: make external of class?
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
        
        if num_workers > 1:
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
                # TODO: delete these compounds instead of just labeling as false?
                if not res[1]:
                    self.compounds[res[0]]['Expand'] = False
                print_progress(i, len(compounds_to_check), 'Tanimoto filter progress:')
        else:
            print(f'Filtering Generation {self.generation}')
            cpd_to_compare = [cpd for cpd in compounds_to_check 
                        if cpd['Generation'] == self.generation]
            for i, cpd in enumerate(cpd_to_compare):
                res = self._compare_to_targets(cpd)
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
        if type(self.crit_tani) == list:
            crit_tani = self.crit_tani[self.generation]
        else:
            crit_tani = self.crit_tani

        try:
            fp1 = utils.get_fp(cpd['SMILES'])
            for fp2 in self.target_fps:
                if AllChem.DataStructs.FingerprintSimilarity(fp1, fp2) >= crit_tani:
                    return (cpd['_id'], True)
        except:
            pass

        return (cpd['_id'], False)

    def prune_network(self, white_list):
        """
        Prune the predicted reaction network to only compounds and reactions
        that terminate in a specified white list of compounds.
        :param white_list: A list of compound_ids to include (if found)
        :type white_list: list
        :return: None
        """
        white_list = set(white_list)
        n_white = len(white_list)
        comp_set, rxn_set = self.find_minimal_set(white_list)
        self.compounds = dict([(k, v) for k, v in self.compounds.items()
                               if k in comp_set])
        self.reactions = dict([(k, v) for k, v in self.reactions.items()
                               if k in rxn_set])
        print(f"""Pruned network to {len(comp_set)} compounds and {len(rxn_set)} reactions based on
                {n_white} whitelisted compounds""")

    def prune_network_to_targets(self):
        """
        Remove compounds that were unexpanded as well as reactions that ended terminally with them.
        """
        print('Pruning to target compounds')
        prune_start = time.time()
        white_list = set()
        for target_smi in self.target_smiles:
            try:
                # generate hash of predicted target compounds
                white_list.add(utils.compound_hash(target_smi, 'Predicted'))
            except:
                pass
        print(f"Identified {len(white_list)} target compounds to filter for.")
        self.prune_network(white_list)
        cpd_to_del = set()
        # Filter out upstream non-targets
        for i in reversed(range(self.generation+1)):
            for cpd_dict in self.compounds.values():
                if cpd_dict['_id'].startswith('C'):
                    if (cpd_dict['Generation'] == i and
                            not cpd_dict['Reactant_in'] and
                            cpd_dict['_id'] not in white_list):
                        # remove reactions
                        for rxn_id in cpd_dict['Product_of']:
                            if rxn_id in self.reactions:
                                del(self.reactions[rxn_id])
                        # remove compound
                        cpd_to_del.add(cpd_dict['_id'])
        for cpd in cpd_to_del:
            del(self.compounds[cpd])

        cpd_count = 0
        for cpd_dict in self.compounds.values():
            if cpd_dict['_id'].startswith('X') or cpd_dict['_id'].startswith('C'):
                cpd_count += 1
        print(f"""Removed upstream non-targets. {cpd_count} compounds remain 
                    and {len(self.reactions)} reactions remain.""")
        print(f"Pruning took {time.time() - prune_start}s")

    def find_minimal_set(self, white_set):
        """
        Given a whitelist this function finds the minimal set of compound and
        reactions ids that comprise the set
        :param white_list:  A list of compound_ids to include (if found)
        :type white_list: list
        :return: compound and reaction id sets
        :rtype: tuple(set, set)
        """
        white_list = list(white_set)
        comp_set = set()
        rxn_set = set()
        for cpd_id in white_list:
            if cpd_id not in self.compounds:
                continue
            for rxn_id in self.compounds[cpd_id]['Product_of']:
                rxn_set.add(rxn_id)
                comp_set.update([x[1] for x
                                 in self.reactions[rxn_id]['Products']])
                for reactant in self.reactions[rxn_id]['Reactants']:
                    comp_set.add(reactant[1])
                    # do not want duplicates or cofactors in the whitelist
                    if reactant[1].startswith('C') and reactant[1] not in white_set:
                        white_list.append(reactant[1])
                        white_set.add(reactant[1])

        # Save targets
        if self.tani_filter:
           for cpd_id in self.compounds:
               if cpd_id.startswith('T'):
                   comp_set.add(cpd_id)

        return comp_set, rxn_set

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
                          key=lambda x: x['_id']):
            rxn['ID_rxn'] = ' + '.join(
                [f"({x[0]}) {self.compounds[x[1]]['ID']}[c0]"
                 for x in rxn['Reactants']]) + ' => ' + ' + '.join(
                     [f"({x[0]}) {self.compounds[x[1]]['ID']}[c0]"
                      for x in rxn['Products']])
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            rxn['ID'] = 'pkr' + str(i).zfill(7)
            i += 1
            self.reactions[rxn['_id']] = rxn

    def write_compound_output_file(self, path, dialect='excel-tab'):
        """Writes all compound data to the specified path.

        :param path: path to output
        :type path: str
        :param dialect: the output format for the file. Choose excel for csv
            excel-tab for tsv.
        :type dialect: str
        """
        path = utils.prevent_overwrite(path)

        columns = ('ID', 'Type', 'Generation', 'Formula', 'InChiKey',
                   'SMILES')
        for _id, val in self.compounds.items():
            inchi_key = AllChem.MolToInchiKey(AllChem.MolFromSmiles(val['SMILES']))
            self.compounds[_id]['InChiKey'] = inchi_key

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

    def save_to_mine(self, num_workers=1, indexing=True, insert_core=True):
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

        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            n = max(n, 1)           
            for i in range(0, len(lst), n):
                yield lst[i:i + n]
            
        print(f'----------------------------------------')
        print(f'Saving results to {self.mine}')
        print(f'----------------------------------------\n')
        start = time.time()
        db = MINE(self.mine, self.mongo_uri)    

        # Insert Reactions    
        print('--------------- Reactions ---------------')
        rxn_start = time.time()
        # Due to memory concerns, reactions are chunked
        # and processed that way. Each batch is calculated
        # in parallel.
        n_rxns = len(self.reactions)
        chunk_size = max(int(n_rxns/100), 10000)     
        print(f"Reaction chunk size writing: {chunk_size}")   
        n_loop = 1
        for rxn_id_chunk in chunks(list(self.reactions.keys()), chunk_size):
            print(f"Writing Reaction Chunk {n_loop} of {round(n_rxns/chunk_size)+1}")
            n_loop += 1
            mine_rxn_requests = self._save_reactions(rxn_id_chunk, db, num_workers)
            if mine_rxn_requests:
                db.reactions.bulk_write(mine_rxn_requests, ordered=False)
        print(f'Finished Inserting Reactions in {time.time() - rxn_start} seconds.')
        print(f'----------------------------------------\n')

        print('--------------- Compounds --------------')
        cpd_start = time.time()
        # Due to memory concerns, compounds are chunked
        # and processed that way. Each batch is calculated
        # in parallel. The memory issue comes from the fact
        # that compounds are inserted in bulk in mongo.
        n_cpds = len(self.compounds)
        chunk_size = max(int(n_cpds/100), 10000)
        print(f"Compound chunk size: {chunk_size}")  
        n_loop = 1
        for cpd_id_chunk in chunks(list(self.compounds.keys()), chunk_size):
            # Insert the three types of compounds
            print(f"Writing Compound Chunk {n_loop} of {round(n_cpds/chunk_size)}")
            n_loop += 1
            core_cpd_requests, core_update_mine_requests, mine_cpd_requests = self._save_compounds(cpd_id_chunk, db, num_workers)            
            if core_cpd_requests:
                if insert_core:
                    db.core_compounds.bulk_write(core_cpd_requests, ordered=False)            
                del(core_cpd_requests)            
            if core_update_mine_requests:
                if insert_core:
                    db.core_compounds.bulk_write(core_update_mine_requests, ordered=False)            
                del(core_update_mine_requests)
            if mine_cpd_requests:
                db.compounds.bulk_write(mine_cpd_requests, ordered=False)                     
                del(mine_cpd_requests) 

        print(f"Finished inserting Compounds to the MINE in {time.time() - cpd_start} seconds.")
        if insert_core:
            db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                            "Action": "Core Compounds Inserted"})
        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
            "Action": "Mine Compounds Inserted"})
        print(f"----------------------------------------\n")

        # Insert target compounds        
        if self.tani_filter:
            target_start = time.time()
            target_cpd_requests = []
            # Write target compounds to target collection
            # Target compounds are written as mine compounds
            print("--------------- Targets ----------------")
            # Insert target compounds
            target_start = time.time()
            # non-parallel insertion
            for comp_dict in self.compounds.values():
                if comp_dict['_id'].startswith('T'):
                    db.insert_mine_compound(comp_dict, target_cpd_requests)     
            print(f"Done with Target Prep--took {time.time() - target_start} seconds.")
            if target_cpd_requests:
                target_start = time.time()
                db.target_compounds.bulk_write(target_cpd_requests, ordered=False)
                print(f"Inserted {len(target_cpd_requests)} Target Compounds in {time.time() - target_start} seconds.")
                del(target_cpd_requests)
                db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                    "Action": "Target Compounds Inserted"})
            else:
                print('No Target Compounds Inserted')
            print(f"----------------------------------------\n")                   

        # Save operators  
        operator_start = time.time()
        if self.operators:
            print("-------------- Operators ---------------")
            # update operator rxn count
            for rxn_dict in self.reactions.values():
                for op in rxn_dict['Operators']:
                    op = op.split("_")[0] # factor in bimolecular rxns
                    self.operators[op][1]['Reactions_predicted'] += 1
            db.operators.insert_many([op[1] for op in self.operators.values()])        
            db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                    "Action": "Operators Inserted"})  
            print(f"Done with Operators Overall--took {time.time() - operator_start} seconds.")
        print(f"----------------------------------------\n")

        if indexing:
            print("-------------- Indices ---------------")
            index_start = time.time()
            db.build_indexes()
            print(f"Done with Indices--took {time.time() - index_start} seconds.")
            print(f"----------------------------------------\n")

        print("-------------- Overall ---------------")
        print(f"Finished uploading everything in {time.time() - start} sec")
        print(f"----------------------------------------\n")

    def _save_compounds(self, cpd_ids, db, num_workers=1):
        # Function to save a given list of compound ids and then delete them from memory
        core_cpd_requests = []
        core_update_mine_requests = []
        mine_cpd_requests = []

        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(cpd_ids) / (num_workers * 10)), 1])
            pool = multiprocessing.Pool(processes=num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    self._save_compound_helper, cpd_ids, chunk_size)):
                if res:
                    mine_cpd_requests.append(res[0])
                    core_update_mine_requests.append(res[1])
                    core_cpd_requests.append(res[2])    
        else:
            # non-parallel insertion
            # Write generated compounds to MINE and core compounds to core
            for cpd_id in cpd_ids:
                comp_dict = self.compounds[cpd_id]
                # Write everything except for targets
                if not comp_dict['_id'].startswith('T'):
                    # These functions are in the MINE class. The request list is
                    # passed and appended in the MINE method.
                    db.insert_mine_compound(comp_dict, mine_cpd_requests)
                    db.update_core_compound_MINES(comp_dict, core_update_mine_requests)
                    db.insert_core_compound(comp_dict, core_cpd_requests)
        return core_cpd_requests, core_update_mine_requests, mine_cpd_requests
    
    def _save_compound_helper(self, cpd_id):
        # Helper function to aid parallelization of saving compounds in
        # save_to_mine
        comp_dict = self.compounds[cpd_id]
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
    
    def _save_reactions(self, rxn_ids, db, num_workers=1):
        # Function to save a given list of compound ids and then delete them from memory
        mine_rxn_requests = []
        rxns_to_write = [self.reactions[rxn_id] for rxn_id in rxn_ids]

        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(rxn_ids) / (num_workers * 10)), 1])
            pool = multiprocessing.Pool(processes=num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    databases.insert_reaction, rxns_to_write, chunk_size)):
                if res:
                    mine_rxn_requests.append(res)

        else:
            # non-parallel insertion
            # Write generated compounds to MINE and core compounds to core
            for rxn_id in rxn_ids:
                rxn = self.reactions[rxn_id]
                db.insert_reaction(rxn, requests=mine_rxn_requests)
        
        return mine_rxn_requests

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

    def _transform_helper(self, compound_smiles, num_workers):
        """Transforms compounds externally of class"""
        def update_cpds_rxns(new_cpds, new_rxns):
            # Save results to self.compounds / self.reactions 
            # ensuring there are no collisions and updating information if there are
            for cpd_id, cpd_dict in new_cpds.items():
                if cpd_id not in self.compounds:
                    self.compounds[cpd_id] = cpd_dict
            
            for rxn_id, rxn_dict in new_rxns.items():
                if rxn_id not in self.reactions:
                    self.reactions[rxn_id] = rxn_dict
                else:
                    self.reactions[rxn_id]['Operators'] = self.reactions[rxn_id]['Operators'].union(rxn_dict['Operators'])

                # Update compound tracking
                for product_id in [cpd_id for _, cpd_id in rxn_dict['Products'] if cpd_id.startswith('C')]:
                    if rxn_id not in self.compounds[product_id]['Product_of']:
                        #TODO make set
                        self.compounds[product_id]['Product_of'].append(rxn_id)
                
                for reactant_id in [cpd_id for _, cpd_id in rxn_dict['Reactants'] if cpd_id.startswith('C')]:
                    if rxn_id not in self.compounds[reactant_id]['Reactant_in']:
                        self.compounds[reactant_id]['Reactant_in'].append(rxn_id)

        # to pass coreactants externally
        coreactant_dict = {co_key: self.compounds[co_key] for _, co_key in self.coreactants.values()}

        new_cpds, new_rxns = _transform_all_compounds_with_full(compound_smiles, self.coreactants, 
            coreactant_dict, self.operators, self.generation, self.explicit_h, num_workers)

        update_cpds_rxns(new_cpds, new_rxns)
    
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

################################################################################
########## Functions to run transformations
#   1. Full Operators are the operators as loaded directly from the list of operatores.
#       These operators use a single supplied molecule for all "Any" instances in the rule.

def _run_reaction(rule_name, rule, reactant_mols, coreactant_mols, coreactant_dict, local_cpds, local_rxns, generation, explicit_h):
    # function to transform a given rule and list of reactant mols
    def _make_half_rxn(mol_list, rules):
        cpds = {}
        cpd_counter = collections.Counter()

        for mol, rule in zip(mol_list, rules):
            if rule == 'Any':
                cpd_dict = _gen_compound(mol)
                # failed compound
                if cpd_dict == None:
                    return None, None
            else:
                cpd_id = coreactant_mols[rule][1]
                cpd_dict = coreactant_dict[cpd_id]            

            cpds[cpd_dict['_id']] = cpd_dict
            cpd_counter.update({cpd_dict['_id']:1})
            
        atom_counts = collections.Counter()
        for cpd_id, cpd_dict in cpds.items():
            for atom_id, atom_count in cpd_dict['atom_count'].items():
                atom_counts[atom_id] += atom_count*cpd_counter[cpd_id]
                
        return [(stoich, cpds[cpd_id]) for cpd_id, stoich in cpd_counter.items()], atom_counts

    def _gen_compound(mol):
        try:
            if explicit_h:
                mol = AllChem.RemoveHs(mol)
            AllChem.SanitizeMol(mol)
        except:
            return None

        mol_smiles = AllChem.MolToSmiles(mol, True)     
        if '.' in mol_smiles:
            return None
                   
        cpd_id = utils.compound_hash(mol_smiles, 'Predicted')
        if cpd_id not in local_cpds:
            cpd_dict = {'ID': None, '_id': cpd_id, 'SMILES': mol_smiles,
                            'Type': 'Predicted',
                            'Generation': generation,
                            'atom_count': utils._getatom_count(mol),
                            'Reactant_in': [], 'Product_of': [],
                            'Expand': True,
                            'Formula': CalcMolFormula(mol)}
        else:
            cpd_dict = local_cpds[cpd_id]
    
        return cpd_dict

    product_sets = rule[0].RunReactants(reactant_mols)

    reactants, reactant_atoms = _make_half_rxn(reactant_mols, rule[1]['Reactants'])      
             
    if not reactants:
        return local_cpds, local_rxns

    for product_mols in product_sets:
        try:
            products, product_atoms = _make_half_rxn(product_mols, rule[1]['Products'])
            if not products:
                continue
            
            if (reactant_atoms - product_atoms or product_atoms - reactant_atoms):
                is_atom_balanced = False
            else:
                is_atom_balanced = True

            if is_atom_balanced:
                for _, cpd_dict in products:
                    if cpd_dict['_id'].startswith('C'):
                        local_cpds.update({cpd_dict['_id']:cpd_dict})

                rhash, rxn_text = utils.rxn2hash(reactants, products)
                if rhash not in local_rxns:
                    local_rxns[rhash] = {'_id': rhash,
                                            # give stoich and id of reactants/products
                                            'Reactants': [(s, r['_id']) for s, r in reactants],
                                            'Products': [(s, p['_id']) for s, p in products],
                                            'Operators': {rule_name},                                    
                                            'SMILES_rxn': rxn_text}
                else:
                    local_rxns[rhash]['Operators'].add(rule_name)                                                                                 
    
        except (ValueError, MemoryError) as e:
            continue
    # return compounds and reactions to be added into the local
    return local_cpds, local_rxns

########## Full Operators
def _transform_ind_compound_with_full(coreactant_mols, coreactant_dict, operators, generation, explicit_h, compound_smiles):
    local_cpds = dict()
    local_rxns = dict()

    mol = AllChem.MolFromSmiles(compound_smiles)
    mol = AllChem.RemoveHs(mol)
    if not mol:        
        print(f"Unable to parse: {compound_smiles}")
        return None
    AllChem.Kekulize(mol, clearAromaticFlags=True)
    if explicit_h:
        mol = AllChem.AddHs(mol)
    # Apply reaction rules to prepared compound

    # run through the single compound operatores
    for rule_name, rule in operators.items():
        # Get RDKit Mol objects for reactants
        reactant_mols = tuple([mol if x == 'Any'
                                else coreactant_mols[x][0]
                                for x in rule[1]['Reactants']])
        # Perform chemical reaction on reactants for each rule
        # try:
        generated_cpds, generated_rxns = _run_reaction(rule_name, rule, reactant_mols, coreactant_mols, coreactant_dict, local_cpds, local_rxns, generation, explicit_h)
        # This error should be addressed in a new version of RDKit
        # TODO: Check this claim
        # except:
        #     print("Runtime ERROR!" + rule_name)
        #     print(compound_smiles)
        #     continue 
        local_cpds.update(generated_cpds)
        for rxn, vals in generated_rxns.items():
            if rxn in local_rxns:
                local_rxns[rxn]['Operators'].union(vals['Operators'])

    return local_cpds,local_rxns

def _transform_all_compounds_with_full(compound_smiles, coreactants, coreactant_dict, operators, generation, explicit_h,
                                    num_workers):
    """
    This function is made to reduce the memory load of parallelization.
    This function accepts in a list of cpds (cpd_list) and runs the transformation in parallel of these.
    """
    def print_progress(done, total):
            # Use print_on to print % completion roughly every 2.5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(.025 * total), 1)
            if not done % print_on:
                print(f"Generation {generation}: {round(done / total * 100)} percent complete")

    # First transform 
    new_cpds_master = {}
    new_rxns_master = {}

    transform_compound_partial = partial(_transform_ind_compound_with_full, coreactants, coreactant_dict, operators, generation, explicit_h)
    # par loop
    if num_workers > 1:
        # TODO chunk size?
        chunk_size = max(
                [round(len(compound_smiles) / (num_workers)), 1])
        print(f'Chunk size = {chunk_size}')
        pool = multiprocessing.Pool(processes=num_workers)
        for i, res in enumerate(pool.imap_unordered(
                            transform_compound_partial, compound_smiles, chunk_size)):
            new_cpds, new_rxns = res
            new_cpds_master.update(new_cpds)
            
            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    new_rxns_master[rxn]['Operators'].union(rxn_dict['Operators'])
                else:
                    new_rxns_master.update({rxn:rxn_dict})
            print_progress(i, len(compound_smiles))
        

    else:
        for i, smiles in enumerate(compound_smiles):
            new_cpds, new_rxns = transform_compound_partial(smiles)
            # new_cpds as cpd_id:cpd_dict
            # new_rxns as rxn_id:rxn_dict
            new_cpds_master.update(new_cpds)
            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    new_rxns_master[rxn]['Operators'].union(rxn_dict['Operators'])
                else:
                    new_rxns_master.update({rxn:rxn_dict})            
            print_progress(i, len(compound_smiles))

    return new_cpds_master, new_rxns_master

if __name__ == "__main__":
    # Get initial time to calculate execution time at end
    t1 = time.time()  # pylint: disable=invalid-name
    # Parse all command line arguments
    parser = ArgumentParser()  # pylint: disable=invalid-name
    # Core args
    parser.add_argument('-C', '--coreactant_list',
                        default="./tests/data/test_coreactants.tsv",
                        help="Specify a list of coreactants as a "
                             "tab-separated file")
    parser.add_argument('-r', '--rule_list',
                        default="./tests/data/test_reaction_rules.tsv",
                        help="Specify a list of reaction rules as a "
                             "tab-separated file")
    parser.add_argument('-c', '--compound_file',
                        default="./tests/data/test_compounds.tsv",
                        help="Specify a list of starting compounds as a "
                             "tab-separated file")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help="Display RDKit errors & warnings")
    # parser.add_argument('--bnice', action='store_true', default=False,
    #                     help="Set several options to enable compatibility "
    #                          "with bnice operators.")
    parser.add_argument('-H', '--explicit_h', action='store_true', default=True,
                        help="Specify explicit hydrogen for use in reaction rules.")
    parser.add_argument('-k', '--kekulize', action='store_true', default=True,
                        help="Specify whether to kekulize compounds.")
    parser.add_argument('-n', '--neutralise', action='store_true', default=True,
                        help="Specify whether to kekulize compounds.")

    parser.add_argument('-m', '--max_workers', default=1, type=int,
                        help="Set the nax number of processes to spawn to "
                             "perform calculations.")
    parser.add_argument('-g', '--generations', default=1, type=int,
                        help="Set the numbers of time to apply the reaction "
                             "rules to the compound set.")
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help="Silence warnings about imbalenced reactions")
    parser.add_argument('-s', '--smiles', default=None,
                        help="Specify a starting compound as SMILES.")
    # Result args
    parser.add_argument('-p', '--pruning_whitelist', default=None,
                        help="Specify a list of target compounds to prune "
                             "reaction network down")    
    parser.add_argument('-o', '--output_dir', default=".",
                        help="The directory in which to place files")
    parser.add_argument('-d', '--database', default=None,
                        help="The name of the database in which to store "
                             "output. If not specified, data is still written "
                             "as tsv files")    
    parser.add_argument('-i', '--image_dir', default=None,
                        help="Specify a directory to store images of all "
                             "created compounds")
    
    OPTIONS = parser.parse_args()
    pk = Pickaxe(coreactant_list=OPTIONS.coreactant_list,
                 rule_list=OPTIONS.rule_list,
                 errors=OPTIONS.verbose, explicit_h=OPTIONS.explicit_h,
                 kekulize=OPTIONS.kekulize, neutralise=OPTIONS.neutralise,
                 image_dir=OPTIONS.image_dir, database=OPTIONS.database,
                 quiet=OPTIONS.quiet)
    # Create a directory for image output file if it doesn't already exist
    if OPTIONS.image_dir and not os.path.exists(OPTIONS.image_dir):
        os.mkdir(OPTIONS.image_dir)
    # If starting compound specified as SMILES string, then add it
    if OPTIONS.smiles:
        # pylint: disable=protected-access
        pk._add_compound("Start", OPTIONS.smiles, cpd_type='Starting Compound')
    else:
        pk.load_compound_set(compound_file=OPTIONS.compound_file)
    # Generate reaction network
    pk.transform_all(num_workers=OPTIONS.max_workers,
                max_generations=OPTIONS.generations)
    if OPTIONS.pruning_whitelist:
        # pylint: disable=invalid-name,protected-access
        mols = [pk._mol_from_dict(line) for line
                in utils.file_to_dict_list(OPTIONS.pruning_whitelist)]
        pk.prune_network([utils.compound_hash(x) for x in mols if x])
    
    pk.assign_ids()
    pk.write_compound_output_file(OPTIONS.output_dir + '/compounds.tsv')
    pk.write_reaction_output_file(OPTIONS.output_dir + '/reactions.tsv')

    print("Execution took %s seconds." % (time.time() - t1))