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
from argparse import ArgumentParser
from copy import deepcopy

from minedatabase import utils
from minedatabase.databases import MINE
from minedatabase.utils import rxn2hash, stoich_tuple
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToFile, rdMolDraw2D


class Pickaxe:
    def __init__(self, rule_list=None, coreactant_list=None, explicit_h=True,
                 kekulize=True, neutralise=True, errors=True,
                 racemize=False, database=None, image_dir=None):
        """This class generates new compounds from user-specified starting
        compounds using a set of SMARTS-based reaction rules. It may be
        initialized with a text file containing the reaction rules and
        coreactants or this may be done on an ad hock basis.

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
        """
        self.rxn_rules = {}
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
        self.radical_check = False
        # Make sure that if a database is to be used, that the database is empty
        if database:
            self.mine = database
            db = MINE(database)
            if db.compounds.count():
                print("Warning: expansion will overwrite existing compounds and"
                      " operators!")
        else:
            self.mine = None

        # Use RDLogger to catch errors in log file. SetLevel indicates mode (
        # 0 - debug, 1 - info, 2 - warning, 3 - critical). Default is no errors.
        from rdkit import RDLogger
        lg = RDLogger.logger()
        if not errors:
            lg.setLevel(4)

        # Load coreactants (if any) into Pickaxe object
        if coreactant_list:
            with open(coreactant_list) as infile:
                for coreactant in infile:
                    self._load_coreactant(coreactant)

        # Load rules (if any) into Pickaxe object
        if rule_list:
            self.load_rxn_rules(rule_list)

    def _load_coreactant(self, coreactant_text):
        """
        Loads a coreactant into the coreactant dictionary from a tab-delimited
            string
        :param coreactant_text: tab-delimited string with the compound name and
            SMILES
        :type coreactant_text: basestring
        """
        # If coreactant is commented out (with '#') then don't import
        if coreactant_text[0] == "#":
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
            raise ValueError("Unable to load coreactant: %s" % coreactant_text)
        _id = self._add_compound(split_text[0], smi, mol=mol, type='Coreactant')
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

    def load_rxn_rules(self, rule_path):
        """Loads all reaction rules from file_path into rxn_rule dict.
        
        :param rule_path: path to file
        :type rule_path: str
        """
        with open(rule_path) as infile:
            # Get all reaction rules from tsv file and store in dictionary (rdr)
            rdr = csv.DictReader((row for row in infile if not
                                  row.startswith('#')), delimiter='\t')
            for rule in rdr:
                try:
                    # Get reactants and products for each reaction into list
                    # form (not ; delimited string)
                    rule['Reactants'] = rule['Reactants'].split(';')
                    rule['Products'] = rule['Products'].split(';')
                    # Ensure that all coreactants are known and accounted for
                    for coreactant_name in rule['Reactants']+rule['Products']:
                        if coreactant_name not in self.coreactants and \
                                        coreactant_name != "Any":
                            raise ValueError('Undefined coreactant:%s' %
                                             coreactant_name)
                    # Create ChemicalReaction object from SMARTS string
                    rxn = AllChem.ReactionFromSmarts(rule['SMARTS'])
                    rule.update({"_id": rule["Name"], "Reactions_predicted": 0})
                    # Ensure that we have number of expected reactants for
                    # each rule
                    if rxn.GetNumReactantTemplates() != len(rule['Reactants']) \
                            or rxn.GetNumProductTemplates() != \
                            len(rule['Products']):
                        raise ValueError("Number of coreactants does not "
                                         "match supplied reaction rule")
                    if rule["Name"] in self.rxn_rules:
                        raise ValueError("Duplicate reaction rule name")
                    # Update reaction rules dictionary
                    self.rxn_rules[rule["Name"]] = (rxn, rule)
                except Exception as e:
                    raise ValueError(str(e) + "Failed to parse %s" %
                                     (rule["Name"]))

    def load_compound_set(self, compound_file=None, structure_field='structure',
                          id_field='id', fragmented_mols=False):
        """If a compound file is provided, this function loads the compounds
        into it's internal dictionary. If not, it attempts to find the
        compounds in it's associated MINE database.
        
        :param compound_file: Path to a file containing compounds as tsv
        :type compound_file: basestring
        :param structure_field: the name of the column containing the
            structure incarnation as Inchi or SMILES (Default:'structure')
        :type structure_field: str
        :param id_field: the name of the column containing the desired
            compound ID (Default: 'id)
        :type id_field: str
        :param fragmented_mols: Permit the loading of disconnected molecules
        :type fragmented_mols: bool
        :return: compound SMILES
        :rtype: list
        """
        compound_smiles = []
        if compound_file:
            with open(compound_file) as infile:
                # Get all compounds from compound file and store in dictionary
                reader = csv.DictReader(infile.readlines(), delimiter="\t")
                for line in reader:
                    # Generate Mol object from InChI code if present
                    if "InChI=" in line[structure_field]:
                        mol = AllChem.MolFromInchi(line[structure_field])
                    # Otherwise generate Mol object from SMILES string
                    else:
                        mol = AllChem.MolFromSmiles(line[structure_field])
                    if not mol:
                        if self.errors:
                            print("Unable to Parse %s" % line[structure_field])
                        continue
                    # If compound is disconnected (determined by GetMolFrags
                    # from rdkit) and loading of these molecules is not
                    # allowed (i.e. fragmented_mols == 1), then don't add to
                    # internal dictionary. This is most common when compounds
                    # are salts.
                    if not fragmented_mols and len(AllChem.GetMolFrags(mol)) \
                            > 1:
                        continue
                    # If specified remove charges (before applying reaction
                    # rules later on)
                    if self.neutralise:
                        mol = utils.neutralise_charges(mol)
                    # Add compound to internal dictionary as a starting
                    # compound and store SMILES string to be returned
                    smi = AllChem.MolToSmiles(mol, True)
                    _id = line[id_field]
                    # Do not operate on inorganic compounds
                    if "C" in smi or "c" in smi:
                        AllChem.SanitizeMol(mol)
                        self._add_compound(_id, smi, mol=mol,
                                           type='Starting Compound')
                        compound_smiles.append(smi)
        # If a MINE database is being used instead, search for compounds
        # annotated as starting compounds and return those as a list of
        # SMILES strings
        elif self.mine:
            db = MINE(self.mine)
            for compound in db.compounds.find():
                _id = compound['_id']
                smi = compound['SMILES']
                # Assume unannotated compounds are starting compounds
                if 'type' not in compound:
                    compound['Type'] = 'Starting Compound'
                self._add_compound(_id, smi, type=compound['Type'])
                compound_smiles.append(smi)
        else:
            raise ValueError('No input file or database specified for '
                             'starting compounds')
        print("%s compounds loaded" % len(compound_smiles))
        return compound_smiles

    def _add_compound(self, id, smi, mol=None, type='Predicted'):
        """Adds a compound to the internal compound dictionary"""
        _id = utils.compound_hash(smi, type == 'Coreactant')
        self._raw_compounds[smi] = _id
        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        if _id not in self.compounds:
            if not mol:
                mol = AllChem.MolFromSmiles(smi)
            i_key = AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))
            self.compounds[_id] = {'ID': id, '_id': _id, "SMILES": smi,
                                   'Inchikey': i_key, 'Type': type,
                                   'Generation': self.generation,
                                   'Formula': self._get_atom_count(mol),
                                   'Reactant_in': [], 'Product_of': [],
                                   "Sources": []}
            # Don't track sources of coreactants
            if _id[0] == 'X':
                del self.compounds[_id]['Sources']
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
                    print("Unable to generate image for %s" % smi)
        return _id

    def transform_compound(self, compound_SMILES, rules=None):
        """Perform transformations to a compound returning the products and the
        predicted reactions
        
        :param compound_SMILES: The compound on which to operate represented
            as SMILES
        :type compound_SMILES: string
        :param rules: The names of the reaction rules to apply. If none,
            all rules in the pickaxe's dict will be used.
        :type rules: list
        :return: Transformed compounds as tuple of compound and reaction dicts
        :rtype: tuple
        """
        if not rules:
            rules = self.rxn_rules.keys()
        # Create Mol object from input SMILES string and remove hydrogens
        # (rdkit)
        mol = AllChem.MolFromSmiles(compound_SMILES)
        mol = AllChem.RemoveHs(mol)
        if not mol:
            if self.errors:
                raise ValueError('Unable to parse: %s' % compound_SMILES)
            else:
                print('Unable to parse: %s' % compound_SMILES)
                return
        if self.kekulize:
            AllChem.Kekulize(mol, clearAromaticFlags=True)
        if self.explicit_h:
            mol = AllChem.AddHs(mol)
        for rule_name in rules:
            # Lookup rule in dictionary from rule name
            rule = self.rxn_rules[rule_name]
            # Get RDKit Mol objects for reactants
            reactant_mols = tuple([mol if x == 'Any' else self.coreactants[x][0]
                                   for x in rule[1]['Reactants']])
            # Perform chemical reaction on reactants for each rule
            try:
                product_sets = rule[0].RunReactants(reactant_mols)
            # This error should be addressed in a new version of RDKit
            except RuntimeError:
                print("Runtime ERROR!"+rule_name)
                print(compound_SMILES)
                continue
            # No enumeration for reactant stereoisomers. _make_half_rxn
            # returns a generator, the first element of which is the reactants.
            reactants, reactant_atoms = next(self._make_half_rxn(reactant_mols,
                                                 rule[1]['Reactants']))
            # By sorting the reactant (and later products) we ensure that
            # compound order is fixed.
            reactants.sort()
            for product_mols in product_sets:
                try:
                    for stereo_prods, product_atoms in self._make_half_rxn(
                            product_mols, rule[1]['Products'], self.racemize):
                        # Update predicted compounds list
                        stereo_prods.sort()
                        # Get reaction text (e.g. A + B <==> C + D)
                        text_rxn = self._add_reaction(reactants, rule_name,
                                                      stereo_prods)
                        # If the SMARTS rule is not atom balanced, this check
                        # detects the accidental alchemy.
                        if reactant_atoms - product_atoms \
                                or product_atoms - reactant_atoms:
                            print("Warning: Unbalanced Reaction produced by "
                                  + rule_name)
                            print(text_rxn)
                            print(reactant_atoms, product_atoms)
                except ValueError as e:
                    print(e)
                    print("Error Processing Rule: " + rule_name)
                    continue
        return self.compounds, self.reactions

    def _get_atom_count(self, mol):
        """Takes a set of mol objects and returns a counter with each element
        type in the set"""
        atoms = collections.Counter()
        # Find all strings of the form A# in the molecular formula where A
        # is the element (e.g. C) and # is the number of atoms of that
        # element in the molecule. Pair is of form [A, #]
        for pair in re.findall('([A-Z][a-z]*)(\d*)',
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
                atoms["*"] += 1
        return atoms

    def _add_reaction(self, reactants, rule_name, stereo_prods):
        """Hashes and inserts reaction into reaction dictionary"""
        # Hash reaction text
        rhash = rxn2hash(reactants, stereo_prods)
        # Generate unique hash from InChI keys of reactants and products
        inchi_rxn_hash, text_rxn = \
            self._calculate_rxn_hash_and_text(reactants, stereo_prods)
        # Add reaction to reactions dictionary if not already there
        if rhash not in self.reactions:
            self.reactions[rhash] = {"_id": rhash,
                                     "Reactants": reactants,
                                     "Products": stereo_prods,
                                     "InChI_hash": inchi_rxn_hash,
                                     "Operators": {rule_name},
                                     "Reaction_rules": {rule_name},
                                     "SMILES_rxn": text_rxn,
                                     "Generation": self.generation}
        # Otherwise, update the operators and rules
        else:
            self.reactions[rhash]['Operators'].add(rule_name)
            self.reactions[rhash]['Reaction_rules'].add(rule_name)
        for prod_id in [x.c_id for x in stereo_prods if x.c_id[0] == 'C']:
            self.compounds[prod_id]['Sources'].append(rhash)
        return text_rxn

    def _make_half_rxn(self, mols, rules, split_stereoisomers=False):
        """Takes a list of mol objects for a half reaction, combines like
        compounds and returns a generator for stoich tuples"""
        # Get compound ids from Mol objects, except for coreactants, in which
        #  case we look them up in the coreactant dictionary
        comps = [self._calculate_compound_information(m, split_stereoisomers)
                 if r == 'Any' else (self.coreactants[r][1],)
                 for m, r in zip(mols, rules)]
        # count the number of atoms on a side
        atom_count = collections.Counter()
        for x in comps:
            atom_count += self.compounds[x[0]]['Formula']
        # Remove duplicates from comps
        half_rxns = [collections.Counter(subrxn) for subrxn
                     in itertools.product(*comps)]
        # Yield generator with stoichiometry tuples (of form stoichiometry, id)
        for rxn in half_rxns:
            yield [stoich_tuple(y, x) for x, y in rxn.items()], atom_count

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
            except:
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
            self._raw_compounds[raw] = tuple([self._add_compound(None, s, mol=m)
                                              for s, m in zip(smiles,
                                                              processed_mols)])
        return self._raw_compounds[raw] if isinstance(
            self._raw_compounds[raw], tuple) else (self._raw_compounds[raw],)

    def _calculate_rxn_hash_and_text(self, reactants, products):
        """Calculates a unique reaction hash using inchikeys. First block is
        connectivity only, second block is stereo only"""
        def __get_blocks(tups):
            first_block, second_block, smiles = [], [], []
            for x in tups:
                comp = self.compounds[x.c_id]
                smiles.append('(%s) %s' % (x.stoich, comp['SMILES']))
                if comp["Inchikey"]:
                    # InChI keys are separated by a hyphen, where the first
                    # part is derived from connectivity and the second part
                    # comes from other layers such as stereochemistry
                    split_inchikey = comp["Inchikey"].split('-')
                    if len(split_inchikey) > 1:
                        first_block.append("%s,%s" % (x.stoich,
                                                      split_inchikey[0]))
                        second_block.append("%s,%s" % (x.stoich,
                                                       split_inchikey[1]))
                else:
                    print("No Inchikey for %s" % x.c_id)
            return "+".join(first_block), "+".join(second_block), \
                   " + ".join(smiles)

        reactants.sort()
        products.sort()
        r_1, r_2, r_s = __get_blocks(reactants)
        p_1, p_2, p_s = __get_blocks(products)
        first_block = r_1+'<==>'+p_1
        second_block = r_2+'<==>'+p_2
        smiles_rxn = r_s+' => '+p_s

        return hashlib.sha256(first_block.encode()).hexdigest()+"-"+hashlib.md5(
            second_block.encode()).hexdigest(), smiles_rxn

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
                comp['ID'] = 'pkc'+str(i).zfill(7)
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
                        print("Unable to generate image for %s"
                              % comp['SMILES'])
        i = 1
        for rxn in sorted(self.reactions.values(),
                          key=lambda x: (x['Generation'], x['_id'])):
            rxn['ID_rxn'] = ' + '.join(['(%s) %s[c0]' %
                                        (x.stoich, self.compounds[x.c_id]["ID"])
                                        for x in rxn["Reactants"]]) \
                            + ' => ' + \
                            ' + '.join(['(%s) %s[c0]' %
                                        (x.stoich, self.compounds[x.c_id]["ID"])
                                        for x in rxn["Products"]])
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            rxn['ID'] = 'pkr'+str(i).zfill(7)
            i += 1
            self.reactions[rxn['_id']] = rxn

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
                print("Generation %s: %s percent complete" %
                      (self.generation,
                       round(done / total * 100)))
        while self.generation < max_generations:
            self.generation += 1
            # Use to print out time per generation at end of loop
            ti = time.time()
            n_comps = len(self.compounds)
            n_rxns = len(self.reactions)
            # Get all SMILES strings for compounds
            compound_smiles = [c['SMILES'] for c in self.compounds.values()
                               if c['Generation'] == self.generation - 1
                               and c['Type'] != 'Coreactant']
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
                    print_progress((i+1), len(compound_smiles))
                self.compounds = new_comps
                self.reactions = new_rxns

            else:
                for i, smi in enumerate(compound_smiles):
                    # Perform possible reactions on compound
                    self.transform_compound(smi)
                    print_progress(i, len(compound_smiles))

            print("Generation %s produced %s new compounds and %s new "
                  "reactions in %s sec" %
                  (self.generation, len(self.compounds)-n_comps,
                   len(self.reactions) - n_rxns, time.time()-ti))

    def prune_network(self, white_list):
        comp_set, rxn_set = self.find_minimal_set(white_list)
        self.compounds = dict([(k, v) for k, v in self.compounds.items() if k in comp_set])
        self.reactions = dict([(k, v) for k, v in self.reactions.items() if k in rxn_set])

    def find_minimal_set(self, white_list):
        # make an ordered set
        white_set = set(white_list)
        comp_set = set()
        rxn_set = set()
        for c_id in white_list:
            if c_id not in self.compounds:
                continue
            for r_id in self.compounds[c_id]['Sources']:
                rxn_set.add(r_id)
                comp_set.update([x.c_id for x in self.reactions[r_id]['Products']])
                for tup in self.reactions[r_id]['Reactants']:
                    comp_set.add(tup.c_id)
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
        columns = ('ID', 'Type', 'Generation', 'Inchikey', 'SMILES')
        with open(path, 'w') as outfile:
            w = csv.DictWriter(outfile, columns, dialect=dialect,
                               extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            w.writerows(sorted(self.compounds.values(), key=lambda x: x['ID']))

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
                                              rxn["SMILES_rxn"], rxn['_id'],
                                              ';'.join(rxn['Reaction_rules'])])
                              + '\n')

    def save_to_MINE(self, db_id):
        """Save compounds to a MINE database.
        
        :param db_id: The name of the target database
        :type db_id: basestring
        """
        db = MINE(db_id)
        bulk_c = db.compounds.initialize_unordered_bulk_op()
        bulk_r = db.reactions.initialize_unordered_bulk_op()

        # This loop performs 4 functions to reactions:
        #   1. Convert stoich_tuples to dicts with hashes
        #   2. Add reaction links to compounds
        #   3. Add source information to compounds
        #   4. Iterate the reactions predicted for each relevant reaction rule
        for rxn in self.reactions.values():
            for x in rxn['Reactants']:
                self.compounds[x.c_id]['Reactant_in'].append(rxn['_id'])
            for x in rxn['Products']:
                self.compounds[x.c_id]['Product_of'].append(rxn['_id'])
                # Don't track sources of coreactants
                if x.c_id[0] == 'X':
                    continue
                self.compounds[x.c_id]['Sources'].append(
                    {"Compounds": [x.c_id for x in rxn['Reactants']],
                     "Operators": list(rxn["Operators"])})
            # Iterate the number of reactions predicted
            for op in rxn['Reaction_rules']:
                self.rxn_rules[op][1]['Reactions_predicted'] += 1
            db.insert_reaction(rxn, bulk=bulk_r)
        if self.reactions:
            bulk_r.execute()
            db.meta_data.insert({"Timestamp": datetime.datetime.now(),
                                 "Action": "Reactions Inserted"})

        for comp_dict in self.compounds.values():
            db.insert_compound(AllChem.MolFromSmiles(comp_dict['SMILES']),
                               comp_dict, bulk=bulk_c)
        bulk_c.execute()
        db.meta_data.insert({"Timestamp": datetime.datetime.now(),
                             "Action": "Compounds Inserted"})

        for x in self.rxn_rules.values():
            # There are fewer reaction rules so bulk operations are not
            # really faster.
            db.operators.save(x[1])
        db.build_indexes()


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
    # unassigned centers (represented by "?" in the second element
    # of the function's return parameters).
    unassigned_centers = [c[0] for c in AllChem.FindMolChiralCenters(
        compound, includeUnassigned=True) if c[1] == "?"]
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
        for atomId, cw in zip(unassigned_centers, seq):
            # cw - Clockwise; ccw - Counterclockwise
            # Get both cw and ccw chiral centers for each center. Used
            # itertools.product to get all combinations.
            if cw:
                compound.GetAtomWithIdx(atomId).SetChiralTag(
                    AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            else:
                compound.GetAtomWithIdx(atomId).SetChiralTag(
                    AllChem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        # Duplicate C++ object so that we don't get multiple pointers to
        # same object
        new_comps.append(deepcopy(compound))
    return new_comps


if __name__ == "__main__":
    # Get initial time to calculate execution time at end
    t1 = time.time()
    # Parse all command line arguments
    parser = ArgumentParser()
    parser.add_argument('-C', '--coreactant_list',
                        default="tests/data/test_coreactants.tsv",
                        help="Specify a list of coreactants as a tab-separated "
                             "file")
    parser.add_argument('-r', '--rule_list',
                        default="tests/data/test_reaction_rules.tsv",
                        help="Specify a list of reaction rules as a "
                             "tab-separated file")
    parser.add_argument('-c', '--compound_file',
                        default="tests/data/test_compounds.tsv",
                        help="Specify a list of starting compounds as a "
                             "tab-separated file")
    parser.add_argument('-p', '--pruning_whitelist', default=None,
                        help="Specify a list of target compounds to prune reaction network down")
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
                        help="Set several options to enable compatibility with "
                             "bnice operators.")
    parser.add_argument('-m', '--max_workers', default=1, type=int,
                        help="Set the nax number of processes to spawn to "
                             "perform calculations.")
    parser.add_argument('-g', '--generations', default=1, type=int,
                        help="Set the numbers of time to apply the reaction "
                             "rules to the compound set.")
    parser.add_argument('-i', '--image_dir', default=None,
                        help="Specify a directory to store images of all "
                             "created compounds")
    options = parser.parse_args()
    pk = Pickaxe(coreactant_list=options.coreactant_list,
                 rule_list=options.rule_list, racemize=options.racemize,
                 errors=options.verbose, explicit_h=options.bnice,
                 kekulize=options.bnice, neutralise=options.bnice,
                 image_dir=options.image_dir, database=options.database)
    # Create a directory for image output file if it doesn't already exist
    if options.image_dir and not os.path.exists(options.image_dir):
        os.mkdir(options.image_dir)
    # If starting compound specified as SMILES string, then add it
    if options.smiles:
        pk._add_compound("Start", options.smiles, type='Starting Compound')
    else:
        pk.load_compound_set(compound_file=options.compound_file)
    # Generate reaction network
    pk.transform_all(max_generations=options.generations,
                     num_workers=options.max_workers)
    if options.pruning_whitelist:
        pk.prune_network(utils.file_to_id_list(options.pruning_whitelist))
    # Save to database (e.g. Mongo) if present, otherwise create output file
    if options.database:
        print("Saving results to %s" % options.database)
        pk.save_to_MINE(options.database)
    else:
        pk.assign_ids()
        pk.write_compound_output_file(options.output_dir+'/compounds.tsv')
        pk.write_reaction_output_file(options.output_dir+'/reactions.tsv')

    print("Execution took %s seconds." % (time.time()-t1))
