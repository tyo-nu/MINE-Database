"""Pickaxe.py: This module generates new compounds from user-specified starting
   compounds using a set of SMARTS-based reaction rules."""
import multiprocessing
import collections
import datetime
import copy
import time
import csv
import os

from argparse import ArgumentParser
from functools import partial
from sys import exit

import minedatabase.databases as databases

from minedatabase.databases import MINE
from minedatabase import utils

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.Draw import MolToFile, rdMolDraw2D
from rdkit.Chem.AllChem import (
    MolFromInchi, MolToInchiKey, MolToSmiles, MolFromSmiles, ReactionFromSmarts,
    GetMolFrags, SanitizeMol, AddHs, Kekulize, RemoveHs, RDKFingerprint
    )

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
                 image_dir=None, quiet=False, retro=False, react_targets=True,
                 filter_after_final_gen=True):
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
        # Main pickaxe properties
        self.operators = {}
        self.coreactants = {}
        self._raw_compounds = {}
        self.compounds = {}
        self.reactions = {}
        self.generation = 0
        # Chemistry options
        self.explicit_h = explicit_h
        self.kekulize = kekulize
        self.racemize = racemize
        self.neutralise = neutralise
        self.fragmented_mols = False
        self.radical_check = False
        # Other options
        self.structure_field = None
        self.image_dir = image_dir
        self.errors = errors
        self.quiet = quiet
        # For filtering
        self.filters = []
        self.targets = dict()
        self.target_smiles = []
        self.target_fps = []
        self.retro = retro
        self.react_targets = react_targets
        self.filter_after_final_gen = filter_after_final_gen
        # database info
        self.mongo_uri = mongo_uri
        # partial_operators
        self.use_partial = False
        self.partial_operators = dict()

        print("----------------------------------------")
        print("Intializing pickaxe object")
        if database:
            # Determine if a specified database is legal
            db = MINE(database, self.mongo_uri)
            if database in db.client.list_database_names():
                if database_overwrite:
                    # If db exists, remove db from all of core compounds
                    # and drop db
                    print((f"Database {database} already exists. "
                           "Deleting database and removing from core compound"
                           " mines."))
                    db.core_compounds.update_many(
                            {},
                            {'$pull': {'MINES': database}}
                        )
                    db.client.drop_database(database)
                    self.mine = database
                else:
                    print((f"Warning! Database {database} already exists."
                           "Specify database_overwrite as true to delete "
                           "old database and write new."))
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
        print("----------------------------------------\n")

    def load_targets(self, target_compound_file, structure_field=None,
                     id_field='id', calc_fp=True):
        """Loads the target list into an list of fingerprints to later compare to
        compounds to determine if those compounds should be expanded.

        :param structure_field: the name of the column containing the
            structure incarnation as Inchi or SMILES (Default:'structure')
        :type structure_field: str
        :param id_field: the name of the column containing the desired
            compound ID (Default: 'id)
        :type id_field: str
        """
        for target_dict in utils.file_to_dict_list(target_compound_file):
            mol = self._mol_from_dict(target_dict, structure_field)
            if not mol:
                continue
            # Add compound to internal dictionary as a target
            # compound and store SMILES string to be returned
            smi = MolToSmiles(mol, True)
            cpd_name = target_dict[id_field]
            # Only operate on organic compounds
            if 'c' in smi.lower():
                SanitizeMol(mol)
                self._add_compound(cpd_name, smi, 'Target Compound', mol)
                self.target_smiles.append(smi)
                if calc_fp:
                    # Generate fingerprints for tanimoto filtering
                    fp = RDKFingerprint(mol)
                    self.target_fps.append(fp)

        print(f"{len(self.target_smiles)} target compounds loaded\n")

    def load_compound_set(self, compound_file=None, id_field='id'):
        """If a compound file is provided, this function loads the compounds
        into its internal dictionary.

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
                smi = MolToSmiles(mol, True)
                cpd_name = cpd_dict[id_field]
                # Do not operate on inorganic compounds
                if 'C' in smi or 'c' in smi:
                    SanitizeMol(mol)
                    self._add_compound(cpd_name, smi,
                                       cpd_type='Starting Compound', mol=mol)
                    compound_smiles.append(smi)

        else:
            raise ValueError("No input file specified for "
                             "starting compounds")

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
            mol = MolFromSmiles(split_text[2])
            if not mol:
                raise ValueError
            # TODO: what do do about stereochemistry? Original comment is below
            # but stereochem was taken out (isn't it removed later anyway?)
            # # Generate SMILES string with stereochemistry taken into account
            smi = MolToSmiles(mol)
        except (IndexError, ValueError):
            raise ValueError(f"Unable to load coreactant: {coreactant_text}")
        cpd_id = self._add_compound(split_text[0], smi, 'Coreactant', mol)
        # If hydrogens are to be explicitly represented, add them to the Mol
        # object
        if self.explicit_h:
            mol = AddHs(mol)
        # If kekulization is preferred (no aromatic bonds, just 3 C=C bonds
        # in a 6-membered aromatic ring for example)
        if self.kekulize:
            Kekulize(mol, clearAromaticFlags=True)
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
                            raise ValueError("Undefined coreactant:"
                                             f"{coreactant_name}")
                    # Create ChemicalReaction object from SMARTS string
                    rxn = ReactionFromSmarts(rule['SMARTS'])
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
                    raise ValueError(f"{str(e)}\nFailed to parse"
                                     f"{rule['Name']}")
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
            mol = MolFromInchi(input_dict[structure_field])
        # Otherwise generate Mol object from SMILES string
        else:
            mol = MolFromSmiles(input_dict[structure_field])
        if not mol:
            if self.errors:
                print(f"Unable to Parse {input_dict[structure_field]}")
            return
        # If compound is disconnected (determined by GetMolFrags
        # from rdkit) and loading of these molecules is not
        # allowed (i.e. fragmented_mols == 1), then don't add to
        # internal dictionary. This is most common when compounds
        # are salts.
        if not self.fragmented_mols and len(GetMolFrags(mol)) > 1:
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
        if cpd_id:
            self._raw_compounds[smi] = cpd_id
            # We don't want to overwrite the same compound from a prior
            # generation so we check with hashed id from above
            if cpd_id not in self.compounds:
                if not mol:
                    mol = MolFromSmiles(smi)
                # expand only Predicted and Starting_compounds
                expand = cpd_type in ['Predicted', 'Starting Compound']
                cpd_dict = {'ID': cpd_name, '_id': cpd_id, 'SMILES': smi,
                            'Type': cpd_type,
                            'Generation': self.generation,
                            'atom_count': utils._getatom_count(
                                                    mol,
                                                    self.radical_check
                                                ),
                            'Reactant_in': [], 'Product_of': [],
                            'Expand': expand,
                            'Formula': CalcMolFormula(mol),
                            'last_tani': 0
                            }
                if cpd_id.startswith('X'):
                    del(cpd_dict['Reactant_in'])
                    del(cpd_dict['Product_of'])
            else:
                cpd_dict = self.compounds[cpd_id]

            return cpd_id, cpd_dict
        else:
            return

    def _add_compound(self, cpd_name, smi, cpd_type, mol=None):
        """Adds a compound to the internal compound dictionary"""
        cpd_id = utils.compound_hash(smi, cpd_type)
        if cpd_id:
            self._raw_compounds[smi] = cpd_id
            # We don't want to overwrite the same compound from a prior
            # generation so we check with hashed id from above
            if cpd_type == 'Target Compound':
                if cpd_id not in self.targets:
                    _, cpd_dict = self._gen_compound(cpd_name,
                                                     smi, cpd_type, mol)
                    self.targets[cpd_id] = cpd_dict

            else:
                if cpd_id not in self.compounds:
                    _, cpd_dict = self._gen_compound(cpd_name,
                                                     smi, cpd_type, mol)
                    self.compounds[cpd_id] = cpd_dict

                    if self.image_dir and self.mine:
                        try:
                            with open(os.path.
                                      join(self.image_dir, cpd_id + '.svg'),
                                      'w') as outfile:

                                mol = MolFromSmiles(cpd_dict['SMILES'])
                                nmol = rdMolDraw2D.PrepareMolForDrawing(mol)
                                d2d = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
                                d2d.DrawMolecule(nmol)
                                d2d.FinishDrawing()
                                outfile.write(d2d.GetDrawingText())
                        except OSError:
                            print("Unable to generate image for "
                                  f"{cpd_dict['SMILES']}")

                return cpd_id
        else:
            return None

    def transform_all(self, num_workers=1, max_generations=1):
        """This function applies all of the reaction rules to all the compounds
        until the generation cap is reached.

        :param num_workers: The number of CPUs to for the expansion process.
        :type num_workers: int
        :param max_generations: The maximum number of times a reaction rule
            may be applied
        :type max_generations: int
        """
        def print_progress(done, total):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(0.1 * total), 1)
            if not done % print_on:
                print((f"Generation {self.generation}: "
                       f"{round(done / total * 100)} percent complete"))

        while self.generation < max_generations or (self.generation == max_generations and self.filter_after_final_gen):

            for _filter in self.filters:
                _filter.apply_filter(self, num_workers)

            if self.generation < max_generations:
                print('----------------------------------------')
                print(f'Expanding Generation {self.generation + 1}\n')

                # Starting time for expansion
                time_init = time.time()

                # Tracking compounds formed
                n_comps = len(self.compounds)
                n_rxns = len(self.reactions)

                # Get SMILES to be expanded
                compound_smiles = [cpd['SMILES'] for cpd in self.compounds.values()
                                   if cpd['Generation'] == self.generation
                                   and cpd['Type'] not in ['Coreactant', 'Target Compound']
                                   and cpd['Expand']]
                # No compounds found
                if not compound_smiles:
                    print("No compounds to expand in generation "
                        f"{self.generation + 1}. Finished expanding.")
                    return None

                self._transform_helper(compound_smiles, num_workers)

                print(f"Generation {self.generation + 1} finished in {time.time()-time_init} "
                    "s and contains:")
                print(f"\t\t{len(self.compounds) - n_comps} new compounds")
                print(f"\t\t{len(self.reactions) - n_rxns} new reactions")
                print(f"\nDone expanding Generation: {self.generation + 1}.")
                print("----------------------------------------\n")

            self.generation += 1

    def load_partial_operators(self, mapped_reactions):
        """Generate set of partial operators from a list of mapped reactions
        corresponding to the reaction rules being used.

        :param mapped_reactions: A .csv file with four columns: rule id,
        source, SMARTS, mapping info.
        :type mapped_reactions: file
        """
        # generate partial operators as done in ipynb
        if not self.operators:
            print("Load reaction rules before loading partial operators")
        else:
            with open(mapped_reactions) as f:
                for line in f.readlines():
                    # Grab info from current mapped reaction
                    rule, source, smiles, _ = line.strip('\n').split('\t')
                    # There should be 2 or more reactants derived from
                    # the mapping code The mapped code doesn't include
                    # cofactors, so 2 or more means any;any*
                    exact_reactants = smiles.split('>>')[0]\
                                            .replace(';', '.').split('.')

                    base_rule = rule.split('_')[0]
                    # base rule must be loaded for partial operator to be uused
                    if base_rule in self.operators:
                        op_reactants = self.operators[base_rule][1]['Reactants']
                        if op_reactants.count('Any') >= 2:
                            mapped_reactants = []
                            for i, r in enumerate(op_reactants):
                                if r == 'Any':
                                    mapped_reactants.append(
                                        exact_reactants.pop(0)
                                        )
                                else:
                                    mapped_reactants.append(r)

                            ind_SMARTS = self.operators[base_rule][1]['SMARTS']
                            ind_SMARTS = (ind_SMARTS.split('>>')[0].
                                          split('>>')[0].replace('(', '').
                                          replace(')', '').split('.'))
                            # now loop through and generate dictionary entries
                            for i, r in enumerate(op_reactants):
                                if r != 'Any':
                                    pass
                                else:
                                    # Build entries
                                    fixed_reactants = [
                                        fr if i != j else 'SMARTS_match'
                                        for j, fr in enumerate(mapped_reactants)
                                    ]

                                    bi_rule =  {
                                        'rule': base_rule,
                                        'rule_reaction': rule,
                                        'reactants': fixed_reactants
                                    }
                                    if ind_SMARTS[i] in self.partial_operators:
                                        self.partial_operators[ind_SMARTS[i]].append(bi_rule)
                                    else:
                                        self.partial_operators[ind_SMARTS[i]] = [bi_rule]

    def _filter_partial_operators(self):
        # generate the reactions to specifically expand
        # based on current compounds
        def partial_reactants_exist(partial_rule):
            try:
                rule_reactants = self.operators[partial_rule['rule']][1]['Reactants']
                cofactor = [False if r == 'Any' else True
                            for r in rule_reactants]

                reactant_ids = []
                for is_cofactor, smi in zip(cofactor, partial_rule['reactants']):
                    if is_cofactor:
                        reactant_ids.append(self.coreactants[smi][1])
                    elif smi == 'SMARTS_match':
                        continue
                    else:
                        reactant_ids.append(utils.compound_hash(smi))

                reactants_exist = [r in self.compounds for r in reactant_ids]
                if all(reactants_exist):
                    return True
                else:
                    return False
            except:
                return False

        filtered_partials = dict()
        for SMARTS_match, rules in self.partial_operators.items():
            for rule in rules:
                if partial_reactants_exist(rule):
                    if SMARTS_match in filtered_partials:
                        filtered_partials[SMARTS_match].append(rule)
                    else:
                        filtered_partials[SMARTS_match] = [rule]

        return filtered_partials

    def remove_cofactor_redundancy(self):
        """Checks for cofactors in rxns that were generated by an any;any rules
        and are specified as generated compounds. Removes redundant reactions.
        """
        # Identify compounds who are really cofactors
        cofactors_as_cpds = []
        cofactor_ids = [cofactor[1] for cofactor in self.coreactants.values()]
        for cpd_id in self.compounds:
            if 'X' + cpd_id[1:] in cofactor_ids and cpd_id.startswith('C'):
                cofactors_as_cpds.append(cpd_id)

        # Loop through identified compounds and update
        # reactions/compounds accordingly
        rxns_to_del = set()
        for cpd_id in cofactors_as_cpds:
            rxn_ids = set(self.compounds[cpd_id]['Product_of']
                          + self.compounds[cpd_id]['Reactant_in'])
            rxns_to_del = rxns_to_del.union(rxn_ids)
            # Check and fix reactions as needed
            for rxn_id in rxn_ids:
                rxn = self.reactions[rxn_id]
                # generate products list with replacements
                reactants = []
                products = []
                for s, reactant in rxn['Reactants']:
                    if reactant in cofactors_as_cpds:
                        reactants.append(
                            (s, self.compounds['X' + reactant[1:]])
                        )
                    else:
                        reactants.append((s, self.compounds[reactant]))

                for s, product in rxn['Products']:
                    if product in cofactors_as_cpds:
                        products.append((s, self.compounds['X' + product[1:]]))
                    else:
                        products.append((s, self.compounds[product]))

                cofactor_rxn_id, rxn_text = utils.rxn2hash(reactants, products)

                # Check to see if reaction makes changes
                # Reactions such as
                # '(1) NCCCNc1cc(C(=O)O)ccc1O + (1) O =>
                #       (1) O + (1) NCCCNc1cc(C(=O)O)ccc1O'
                # are possible. Filter out
                sorted_reactants = sorted(reactants, key=lambda x: x[1]['_id'])
                sorted_products = sorted(products, key=lambda x: x[1]['_id'])

                # Update newly calculated reaction info in self.reactions.
                # Three possibilites:
                # 1. Reaction has same products as reactants -- Remove completely
                # 2. Reaction exists already -- append to existing reaction
                # 3. Reaction doesn't exist -- make new reaction

                # After this, delete redudnant compounds and update
                # other compounds in reaction

                append_new = False
                if sorted_reactants == sorted_products:
                    # Reaction is something like CoF + CoF -> Product + Product
                    # Where products are actually just cofactors
                    pass # No need to update anything in self.reactions

                elif cofactor_rxn_id in self.reactions:
                    # Reaction already exists, update reaction operators
                    ops = self.reactions[cofactor_rxn_id]['Operators']
                    ops = ops | self.reactions[rxn_id]['Operators']
                    self.reactions[cofactor_rxn_id]['Operators'] = ops
                    append_new = True

                else:
                    # Reaction does not exist, generate new reaction
                    cofactor_rxn = {
                        '_id': cofactor_rxn_id,
                        # give stoich and id of reactants/products
                        'Reactants': [(s, r['_id']) for s, r in reactants],
                        'Products': [(s, p['_id']) for s, p in products],
                        'Operators': rxn['Operators'],
                        'SMILES_rxn': rxn_text
                        }

                    # Assign reaction to reactions dict
                    self.reactions[cofactor_rxn_id] = cofactor_rxn
                    append_new = True

                for _, cpd in rxn['Reactants']:
                    if cpd.startswith('C'):
                        if rxn_id in self.compounds[cpd]['Reactant_in']:
                            self.compounds[cpd]['Reactant_in'].remove(rxn_id)

                        if (append_new
                            and cofactor_rxn_id
                                not in self.compounds[cpd]['Reactant_in']):
                            self.compounds[cpd]['Reactant_in'].append(cofactor_rxn_id)

                for _, cpd in rxn['Products']:
                    if cpd.startswith('C'):
                        if rxn_id in self.compounds[cpd]['Product_of']:
                            self.compounds[cpd]['Product_of'].remove(rxn_id)

                        if (append_new
                            and cofactor_rxn_id
                                not in self.compounds[cpd]['Product_of']):
                            self.compounds[cpd]['Product_of'].append(cofactor_rxn_id)

        if rxns_to_del:
            for rxn_id in rxns_to_del:
                del(self.reactions[rxn_id])

            for cpd_id in cofactors_as_cpds:
                del(self.compounds[cpd_id])

    def prune_network(self, white_list):
        """
        Prune the predicted reaction network to only compounds and reactions
        that terminate in a specified white list of compounds.

        :param white_list: A list of compound_ids to include (if found)
        :type white_list: list
        :return: None
        """
        n_white = len(white_list)
        cpd_set, rxn_set = self.find_minimal_set(white_list)
        self.compounds = dict([(k, v) for k, v in self.compounds.items()
                               if k in cpd_set])
        self.reactions = dict([(k, v) for k, v in self.reactions.items()
                               if k in rxn_set])
        print(f"Pruned network to {len(cpd_set)} compounds and "
              f"{len(rxn_set)} reactions based on "
              f"{n_white} whitelisted compounds.")

    def prune_network_to_targets(self):
        """
        Prune the predicted reaction network to only compounds and reactions
        that terminate in the target compounds.
        """
        print('----------------------------------------')
        prune_start = time.time()
        white_list = set()
        for target_id in self.targets:
            white_list.add('C' + target_id[1:])

        print(f'Pruning to {len(white_list)} target compounds')
        self.prune_network(white_list)

        print(f"Pruning took {time.time() - prune_start}s")
        print('----------------------------------------\n')

    def find_minimal_set(self, white_list):
        """
        Given a whitelist this function finds the minimal set of compound and
        reactions ids that comprise the set
        :param white_list:  A list of compound_ids to include (if found)
        :type white_list: set
        :return: compound and reaction id sets
        :rtype: tuple(set, set)
        """

        queue = list(white_list)
        visited = set()
        cpd_set = set()
        rxn_set = set()

        while queue:
            cpd_id = queue.pop()
            visited.add(cpd_id)
            if cpd_id not in self.compounds:
                continue

            # Add info for reactions that produce compound
            for rxn_id in self.compounds[cpd_id]['Product_of']:
                rxn_set.add(rxn_id)
                # Add products, not to be further explored
                for cpd in self.reactions[rxn_id]['Products']:
                    cpd_set.add(cpd[1])

                # Add reactants, also add to queue
                for cpd in self.reactions[rxn_id]['Reactants']:
                    cpd_set.add(cpd[1])
                    if cpd[1].startswith('C') and cpd[1] not in visited:
                        queue.append(cpd[1])

        return cpd_set, rxn_set

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
            if not comp.get('ID'):
                comp['ID'] = 'pkc' + str(i).zfill(7)
                i += 1
                self.compounds[comp['_id']] = comp
                # If we are not loading into the mine, we generate the image
                # here.
                if self.image_dir and not self.mine:
                    mol = MolFromSmiles(comp['SMILES'])
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
            inchi_key = MolToInchiKey(MolFromSmiles(val['SMILES']))
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

    def save_to_mine(self, num_workers=1, indexing=True, insert_core=True,
                     insert_targets=True):
        """Save compounds to a MINE database.

        :param num_workers: Number of processors to use.
        :type num_workers: int
        :param indexing: Should mongo indices be made.
        :type indexing: bool
        :param insert_core: Should compounds be inserted into core.
        :type insert_core: bool
        """
        def print_progress(done, total, section):
            # Use print_on to print % completion roughly every 5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(0.1 * total), 1)
            if not (done % print_on):
                print(f"{section} {round(done / total*100)} percent complete")

        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            n = max(n, 1)
            for i in range(0, len(lst), n):
                yield lst[i:i + n]

        print('\n----------------------------------------')
        print(f'Saving results to {self.mine}')
        print('----------------------------------------\n')
        start = time.time()
        db = MINE(self.mine, self.mongo_uri)

        # Insert Reactions
        print('--------------- Reactions --------------')
        rxn_start = time.time()
        # Due to memory concerns, reactions are chunked
        # and processed that way. Each batch is calculated
        # in parallel.
        n_rxns = len(self.reactions)
        chunk_size = max(int(n_rxns/(num_workers*100)), 10000)
        print(f"Reaction chunk size writing: {chunk_size}")
        n_loop = 1
        for rxn_id_chunk in chunks(list(self.reactions.keys()), chunk_size):
            print(f"Writing Reaction Chunk {n_loop} "
                  f"of {round(n_rxns/chunk_size)+1}")
            n_loop += 1
            mine_rxn_requests = self._save_reactions(rxn_id_chunk,
                                                     db, num_workers)
            if mine_rxn_requests:
                db.reactions.bulk_write(mine_rxn_requests, ordered=False)
                del(mine_rxn_requests)
        print("Finished Inserting Reactions in "
              f" {time.time() - rxn_start} seconds.")
        print('----------------------------------------\n')

        print('--------------- Compounds --------------')
        cpd_start = time.time()
        # Due to memory concerns, compounds are chunked
        # and processed that way. Each batch is calculated
        # in parallel.
        n_cpds = len(self.compounds)
        chunk_size = max(int(n_cpds/(num_workers*100)), 10000)
        print(f"Compound chunk size: {chunk_size}")
        n_loop = 1
        for cpd_id_chunk in chunks(list(self.compounds.keys()), chunk_size):
            # Insert the three types of compounds
            print(f"Writing Compound Chunk {n_loop} "
                  f"of {round(n_cpds/chunk_size)+1}")
            n_loop += 1
            res = self._save_compounds(cpd_id_chunk, db, num_workers)
            core_cpd_requests = res[0]
            core_update_mine_requests = res[1]
            mine_cpd_requests = res[2]

            if core_cpd_requests:
                if insert_core:
                    db.core_compounds.bulk_write(core_cpd_requests,
                                                 ordered=False)
                del(core_cpd_requests)
            if core_update_mine_requests:
                if insert_core:
                    db.core_compounds.bulk_write(core_update_mine_requests,
                                                 ordered=False)
                del(core_update_mine_requests)
            if mine_cpd_requests:
                db.compounds.bulk_write(mine_cpd_requests, ordered=False)
                del(mine_cpd_requests)

        print("Finished inserting Compounds to the MINE in "
              f"{time.time() - cpd_start} seconds.")
        if insert_core:
            db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                     "Action": "Core Compounds Inserted"})
        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                 "Action": "Mine Compounds Inserted"})
        print("----------------------------------------\n")

        # Insert target compounds
        if insert_targets and self.targets:
            target_start = time.time()
            target_cpd_requests = []
            # Write target compounds to target collection
            # Target compounds are written as mine compounds
            print("--------------- Targets ----------------")
            # Insert target compounds
            target_start = time.time()
            # non-parallel insertion
            for comp_dict in self.targets.values():
                db.insert_mine_compound(comp_dict, target_cpd_requests)
            print("Done with Target Prep--took "
                  f"{time.time() - target_start} seconds.")
            if target_cpd_requests:
                target_start = time.time()
                db.target_compounds.bulk_write(target_cpd_requests,
                                               ordered=False)
                print(f"Inserted {len(target_cpd_requests)}"
                      f" Target Compounds in {time.time() - target_start}"
                      " seconds.")
                del(target_cpd_requests)
                db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                         "Action": "Target Compounds Inserted"}
                                        )
            else:
                print('No Target Compounds Inserted')
            print("----------------------------------------\n")

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
            print("Done with Operators Overall--took "
                  f"{time.time() - operator_start} seconds.")
        print("----------------------------------------\n")

        if indexing:
            print("-------------- Indices ---------------")
            index_start = time.time()
            db.build_indexes()
            print("Done with Indices--took "
                  f"{time.time() - index_start} seconds.")
            print("----------------------------------------\n")

        print("-------------- Overall ---------------")
        print(f"Finished uploading everything in {time.time() - start} sec")
        print("----------------------------------------\n")

    def _save_compounds(self, cpd_ids, db, num_workers=1):
        # Function to save a given list of compound ids and then
        # delete them from memory
        core_cpd_requests = []
        core_update_mine_requests = []
        mine_cpd_requests = []

        cpd_dicts = [self.compounds[cpd_id] for cpd_id in cpd_ids
                     if not self.compounds[cpd_id]['_id'].startswith('T')]
        _save_compound_helper_partial = partial(_save_compound_helper,
                                                self.mine)
        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(cpd_ids) / (num_workers * 4)), 1])
            pool = multiprocessing.Pool(processes=num_workers)
            for i, res in enumerate(pool.imap_unordered(
                    _save_compound_helper_partial, cpd_dicts, chunk_size)):
                if res:
                    mine_cpd_requests.append(res[0])
                    core_update_mine_requests.append(res[1])
                    core_cpd_requests.append(res[2])
        else:
            # non-parallel insertion
            # Write generated compounds to MINE and core compounds to core
            for cpd_dict in cpd_dicts:
                # These functions are in the MINE database class.
                # The request list is passed and appended in the MINE method.
                db.insert_mine_compound(cpd_dict, mine_cpd_requests)
                db.update_core_compound_MINES(cpd_dict,
                                              core_update_mine_requests)
                db.insert_core_compound(cpd_dict, core_cpd_requests)
        return core_cpd_requests, core_update_mine_requests, mine_cpd_requests

    def _save_reactions(self, rxn_ids, db, num_workers=1):
        # Function to save a given list of compound ids and then delete
        # them from memory
        mine_rxn_requests = []
        rxns_to_write = [self.reactions[rxn_id] for rxn_id in rxn_ids]

        if num_workers > 1:
            # parallel insertion
            chunk_size = max(
                [round(len(rxn_ids) / (num_workers * 4)), 1])
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

    def _transform_helper(self, compound_smiles, num_workers):
        """Transforms compounds externally of class"""
        def update_cpds_rxns(new_cpds, new_rxns):
            # Save results to self.compounds / self.reactions
            # ensuring there are no collisions and updating information
            #  if there are
            for cpd_id, cpd_dict in new_cpds.items():
                if cpd_id not in self.compounds:
                    self.compounds[cpd_id] = cpd_dict

            for rxn_id, rxn_dict in new_rxns.items():
                if rxn_id not in self.reactions:
                    self.reactions[rxn_id] = rxn_dict
                else:
                    rxn_ops = self.reactions[rxn_id]['Operators']
                    rxn_ops = rxn_ops.union(rxn_dict['Operators'])
                    if 'Partial Operators' in self.reactions[rxn_id]:
                        par_ops = self.reactions[rxn_id]['Partial Operators']
                        par_ops = par_ops.union(rxn_dict['Partial Operators'])

                # Update compound tracking
                for product_id in [cpd_id for _, cpd_id in rxn_dict['Products']
                                   if cpd_id.startswith('C')]:
                    if rxn_id not in self.compounds[product_id]['Product_of']:
                        self.compounds[product_id]['Product_of'].append(rxn_id)

                for reactant_id in [
                                  cpd_id for _, cpd_id in rxn_dict['Reactants']
                                  if cpd_id.startswith('C')
                                   ]:
                    if rxn_id not in self.compounds[reactant_id]['Reactant_in']:
                        self.compounds[reactant_id]['Reactant_in'].append(rxn_id)

        # to pass coreactants externally
        coreactant_dict = {co_key: self.compounds[co_key] for _, co_key in
                           self.coreactants.values()}

        new_cpds, new_rxns = _transform_all_compounds_with_full(
                                    compound_smiles, self.coreactants,
                                    coreactant_dict, self.operators,
                                    self.generation + 1, self.explicit_h,
                                    num_workers
                            )

        update_cpds_rxns(new_cpds, new_rxns)

        if self.partial_operators:
            print("\nGenerating partial operators...")
            partial_operators = self._filter_partial_operators()
            if partial_operators:
                print("Found partial operators, applying.")
                # transform partial
                new_cpds, new_rxns = _transform_all_compounds_with_partial(
                                            compound_smiles, self.coreactants,
                                            coreactant_dict, self.operators,
                                            self.generation, self.explicit_h,
                                            num_workers, partial_operators
                                    )

                update_cpds_rxns(new_cpds, new_rxns)
            else:
                print("No partial operators could be generated.")



# def _racemization(compound, max_centers=3, carbon_only=True):
#     """Enumerates all possible stereoisomers for unassigned chiral centers.

#     :param compound: A compound
#     :type compound: rdMol object
#     :param max_centers: The maximum number of unspecified stereocenters to
#         enumerate. Sterioisomers grow 2^n_centers so this cutoff prevents lag
#     :type max_centers: int
#     :param carbon_only: Only enumerate unspecified carbon centers. (other
#         centers are often not tautomeric artifacts)
#     :type carbon_only: bool
#     :return: list of stereoisomers
#     :rtype: list of rdMol objects
#     """
#     new_comps = []
#     # FindMolChiralCenters (rdkit) finds all chiral centers. We get all
#     # unassigned centers (represented by '?' in the second element
#     # of the function's return parameters).
#     unassigned_centers = [c[0] for c in AllChem.FindMolChiralCenters(
#         compound, includeUnassigned=True) if c[1] == '?']
#     # Get only unassigned centers that are carbon (atomic number of 6) if
#     # indicated
#     if carbon_only:
#         unassigned_centers = list(
#             filter(lambda x: compound.GetAtomWithIdx(x).GetAtomicNum() == 6,
#                 unassigned_centers))
#     # Return original compound if no unassigned centers exist (or if above
#     # max specified (to prevent lag))
#     if not unassigned_centers or len(unassigned_centers) > max_centers:
#         return [compound]
#     for seq in itertools.product([1, 0], repeat=len(unassigned_centers)):
#         for atomid, clockwise in zip(unassigned_centers, seq):
#             # Get both cw and ccw chiral centers for each center. Used
#             # itertools.product to get all combinations.
#             if clockwise:
#                 compound.GetAtomWithIdx(atomid).SetChiralTag(
#                     AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
#             else:
#                 compound.GetAtomWithIdx(atomid).SetChiralTag(
#                     AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
#         # Duplicate C++ object so that we don't get multiple pointers to
#         # same object
#         new_comps.append(deepcopy(compound))
#     return new_comps


###############################################################################
# Functions to run transformations
# There are two distinct functions to transform two flavors of operators.
#   1. Full Operators are the operators as loaded directly from the list
#       of operatores. These operators use a single supplied molecule for all
#       "Any" instances in the rule.
#   2. Partial operators are operators for reactions with more than one
#       "Any" in the reactants. These rules are derived from individually
#       mapped reactions and are called partial because only one "Any" is
#       allowed to be novel, the other "Any"s are determined by the
#       mapped reactions.

# Both operators are preprocessed slightly differently, but yield the same
# output format back to the pickaxe object.

# there are way too many variables passed... switch to **kwargs?
# Generic reaction implementation
def _run_reaction(rule_name, rule, reactant_mols, coreactant_mols,
                  coreactant_dict, local_cpds, local_rxns, generation,
                  explicit_h):
    """
    Transform list of mols and a reaction rule into a half reaction
    describing either the reactants or products.
    """
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

        cpd_returns = [(stoich, cpds[cpd_id]) for
                       cpd_id, stoich in cpd_counter.items()]

        return cpd_returns, atom_counts

    def _gen_compound(mol):
        try:
            if explicit_h:
                mol = RemoveHs(mol)
            SanitizeMol(mol)
        except:
            return None

        mol_smiles = MolToSmiles(mol, True)
        if '.' in mol_smiles:
            return None

        cpd_id = utils.compound_hash(mol_smiles, 'Predicted')
        if cpd_id:
            if cpd_id not in local_cpds:
                cpd_dict = {'ID': None, '_id': cpd_id, 'SMILES': mol_smiles,
                                'Type': 'Predicted',
                                'Generation': generation,
                                'atom_count': utils._getatom_count(mol),
                                'Reactant_in': [], 'Product_of': [],
                                'Expand': True,
                                'Formula': CalcMolFormula(mol),
                                'last_tani': 0}
            else:
                cpd_dict = local_cpds[cpd_id]

            return cpd_dict
        else:
            return None

    try:
        product_sets = rule[0].RunReactants(reactant_mols)
        reactants, reactant_atoms = _make_half_rxn(reactant_mols,
                                                   rule[1]['Reactants'])
    except:
        reactants = None

    if not reactants:
        return local_cpds, local_rxns

    for product_mols in product_sets:
        try:
            products, product_atoms = _make_half_rxn(product_mols,
                                                     rule[1]['Products'])
            if not products:
                continue

            if (
                reactant_atoms - product_atoms or
                product_atoms - reactant_atoms
               ):
                is_atom_balanced = False
            else:
                is_atom_balanced = True

            if is_atom_balanced:
                for _, cpd_dict in products:
                    if cpd_dict['_id'].startswith('C'):
                        local_cpds.update({cpd_dict['_id']: cpd_dict})

                rhash, rxn_text = utils.rxn2hash(reactants, products)
                if rhash not in local_rxns:
                    local_rxns[rhash] = {
                            '_id': rhash,
                            # give stoich and id of reactants/products
                            'Reactants': [(s, r['_id']) for s, r in reactants],
                            'Products': [(s, p['_id']) for s, p in products],
                            'Operators': {rule_name},
                            'SMILES_rxn': rxn_text
                            }
                else:
                    local_rxns[rhash]['Operators'].add(rule_name)

        except (ValueError, MemoryError) as e:
            continue
    # return compounds and reactions to be added into the local
    return local_cpds, local_rxns


# Full Operators
def _transform_ind_compound_with_full(coreactant_mols, coreactant_dict,
                                      operators, generation, explicit_h,
                                      compound_smiles):
    logger = RDLogger.logger()
    logger.setLevel(4)
    local_cpds = dict()
    local_rxns = dict()

    mol = MolFromSmiles(compound_smiles)
    mol = RemoveHs(mol)
    if not mol:
        print(f"Unable to parse: {compound_smiles}")
        return None
    Kekulize(mol, clearAromaticFlags=True)
    if explicit_h:
        mol = AddHs(mol)
    # Apply reaction rules to prepared compound

    # run through the single compound operatores
    for rule_name, rule in operators.items():
        # Get RDKit Mol objects for reactants
        reactant_mols = tuple([mol if x == 'Any'
                               else coreactant_mols[x][0]
                               for x in rule[1]['Reactants']])
        # Perform chemical reaction on reactants for each rule
        # try:
        generated_cpds, generated_rxns = _run_reaction(
            rule_name, rule, reactant_mols, coreactant_mols,
            coreactant_dict, local_cpds, local_rxns, generation, explicit_h
        )

        local_cpds.update(generated_cpds)
        for rxn, vals in generated_rxns.items():
            if rxn in local_rxns:
                local_rxns[rxn]['Operators'].union(vals['Operators'])

    return local_cpds, local_rxns


def _transform_all_compounds_with_full(compound_smiles, coreactants,
                                       coreactant_dict, operators, generation,
                                       explicit_h, num_workers):
    """
    This function is made to reduce the memory load of parallelization.
    This function accepts in a list of cpds (cpd_list) and runs the
    transformation in parallel of these.
    """
    def print_progress(done, total):
        # Use print_on to print % completion roughly every 2.5 percent
        # Include max to print no more than once per compound (e.g. if
        # less than 20 compounds)
        print_on = max(round(0.1 * total), 1)
        if not done % print_on:
            print(f"Generation {generation}: {round(done / total * 100)}"
                  " percent complete")

    # First transform
    new_cpds_master = {}
    new_rxns_master = {}

    transform_compound_partial = partial(
        _transform_ind_compound_with_full,
        copy.copy(coreactants),
        copy.copy(coreactant_dict),
        copy.copy(operators),
        copy.copy(generation),
        copy.copy(explicit_h)
    )
    # par loop
    if num_workers > 1:
        # TODO chunk size?
        # chunk_size = max(
        #         [round(len(compound_smiles) / (num_workers)), 1])
        chunk_size = 1
        # print(f'Chunk size = {chunk_size}')
        pool = multiprocessing.Pool(processes=num_workers)
        for i, res in enumerate(pool.imap_unordered(
                            transform_compound_partial,
                            compound_smiles,
                            chunk_size)):
            new_cpds, new_rxns = res
            new_cpds_master.update(new_cpds)

            # Need to check if reactions already exist to update operators list
            for rxn, rxn_dict in new_rxns.items():
                if rxn in new_rxns_master:
                    ops_set = rxn_dict["Operators"]
                    new_rxns_master[rxn]['Operators'].union(ops_set)
                else:
                    new_rxns_master.update({rxn: rxn_dict})
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
                    ops_set = rxn_dict["Operators"]
                    new_rxns_master[rxn]['Operators'].union(ops_set)
                else:
                    new_rxns_master.update({rxn: rxn_dict})
            print_progress(i, len(compound_smiles))

    return new_cpds_master, new_rxns_master


def _save_compound_helper(mine, cpd_dict):
    # Helper function to aid parallelization of saving compounds in
    # save_to_mine
    # These functions are outside of the MINE class in order to
    # allow for parallelization. When in the MINE class it is not
    # serializable with pickle. In comparison to the class functions,
    # these return the requests instead of appending to a passed list.
    mine_req = databases.insert_mine_compound(cpd_dict)
    core_up_req = databases.update_core_compound_MINES(cpd_dict, mine)
    core_in_req = databases.insert_core_compound(cpd_dict)
    return [mine_req, core_up_req, core_in_req]


# Partial Operators
def _transform_ind_compound_with_partial(coreactant_mols, coreactant_dict,
                                         operators, generation, explicit_h,
                                         partial_rules, compound_smiles):
    # 1. See if rule matches the compound passed
    #   (rule from partial_rules dict keys)
    # 2. If match apply transform_ind_compound_with_full to each
    def generate_partial_mols(partial_rule):
        def gen_mol(smi):
            mol = MolFromSmiles(smi)
            mol = RemoveHs(mol)
            Kekulize(mol, clearAromaticFlags=True)
            if explicit_h:
                mol = AddHs(mol)
            return mol

        rule_reactants = operators[partial_rule['rule']][1]['Reactants']
        cofactor = [False if r == 'Any' else True for r in rule_reactants]
        reactant_mols = []
        for is_cofactor, smi in zip(cofactor, partial_rule['reactants']):
                if is_cofactor:
                    reactant_mols.append(coreactant_mols[smi][0])
                elif smi == 'SMARTS_match':
                    reactant_mols.append(gen_mol(compound_smiles))
                else:
                    # These reactions already happen with any;any
                    if (utils.compound_hash(smi) !=
                            utils.compound_hash(compound_smiles)):
                        reactant_mols.append(gen_mol(smi))
                    else:
                        return None
        return reactant_mols

    local_cpds = dict()
    local_rxns = dict()

    mol = MolFromSmiles(compound_smiles)
    mol = RemoveHs(mol)
    if not mol:
        print(f"Unable to parse: {compound_smiles}")
        return None
    Kekulize(mol, clearAromaticFlags=True)
    if explicit_h:
        mol = AddHs(mol)
    # Apply reaction rules to prepared compound

    # run through the single compound operatores
    for ind_SMARTS, rules in partial_rules.items():
        # does mol match vs smiles match change things?
        if QuickSmartsMatch(compound_smiles, ind_SMARTS):
            for partial_rule in rules:
                # Perform chemical reaction on reactants for each rule
                # try:
                rule_name = partial_rule['rule_reaction'].split('_')[0]
                rule = operators[partial_rule['rule']]
                reactant_mols = generate_partial_mols(partial_rule)
                if reactant_mols:
                    generated_cpds, generated_rxns = _run_reaction(rule_name, rule,
                        reactant_mols, coreactant_mols, coreactant_dict, local_cpds, local_rxns, generation, explicit_h)


                    local_cpds.update(generated_cpds)
                    for rxn, vals in generated_rxns.items():
                        if rxn in local_rxns:
                            if 'Partial Operators' in local_rxns[rxn]:
                                local_rxns[rxn]['Partial Operators'].update([partial_rule['rule_reaction']])
                            else:
                                local_rxns[rxn]['Partial Operators'] = set([partial_rule['rule_reaction']])
    return local_cpds,local_rxns


def _transform_all_compounds_with_partial(compound_smiles, coreactants,
                                          coreactant_dict, operators,
                                          generation, explicit_h, num_workers,
                                          partial_rules):
    """
    This function is made to reduce the memory load of parallelization.
    This function accepts in a list of cpds (cpd_list) and runs the transformation in parallel of these.
    """
    def print_progress(done, total):
            # Use print_on to print % completion roughly every 2.5 percent
            # Include max to print no more than once per compound (e.g. if
            # less than 20 compounds)
            print_on = max(round(0.1 * total), 1)
            if not done % print_on:
                print(f"Generation {generation}: {round(done / total * 100)} percent complete")

    # First transform
    new_cpds_master = {}
    new_rxns_master = {}

    transform_compound_partial = partial(_transform_ind_compound_with_partial, coreactants,
                                            coreactant_dict, operators, generation, explicit_h, partial_rules)
    # par loop
    if num_workers > 1:
        # TODO chunk size?
        # chunk_size = max(
        #         [round(len(compound_smiles) / (num_workers)), 1])
        chunk_size = 1
        # print(f'Chunk size = {chunk_size}')
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
                    new_rxns_master[rxn]['Partial Operators'] = new_rxns_master[rxn]['Partial Operators'].union(rxn_dict['Partial Operators'])
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
