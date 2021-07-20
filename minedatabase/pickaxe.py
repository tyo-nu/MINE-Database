"""Pickaxe.py: Create network expansions from reaction rules and compounds.

This module generates new compounds from user-specified starting
compounds using a set of SMARTS-based reaction rules.
"""
import csv
import datetime
import os
import pickle
import time
from argparse import ArgumentParser
from io import StringIO
from pathlib import Path, PosixPath, WindowsPath
from sys import exit
from typing import List, Set, Tuple, Union

from rdkit.Chem.AllChem import (
    AddHs,
    GetMolFrags,
    Kekulize,
    MolFromInchi,
    MolFromSmiles,
    MolToInchiKey,
    MolToSmiles,
    RDKFingerprint,
    ReactionFromSmarts,
    SanitizeMol,
)
from rdkit.Chem.Draw import MolToFile, rdMolDraw2D
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.RDLogger import logger

from minedatabase import utils
from minedatabase.databases import (
    MINE,
    write_compounds_to_mine,
    write_core_compounds,
    write_reactions_to_mine,
    write_targets_to_mine,
)
from minedatabase.reactions import transform_all_compounds_with_full


# Default to no errors
lg = logger()
lg.setLevel(4)


class Pickaxe:
    """Class to generate expansions with compounds and reaction rules.

    This class generates new compounds from user-specified starting
    compounds using a set of SMARTS-based reaction rules. It may be initialized
    with a text file containing the reaction rules and coreactants or this may
    be done on an ad hoc basis.

    Parameters
    ----------
    rule_list : str
        Filepath of rules.
    coreactant_list : str
        Filepath of coreactants.
    explicit_h : bool, optional
        Whether rules utilize explicit hydrogens, by default True.
    kekulize : bool, optional
        Whether or not to kekulize compounds before reaction,
        by default True.
    neutralise : bool, optional
        Whether or not to neutralise compounds, by default True.
    errors : bool, optional
        Whether or not to print errors to stdout, by default True.
    inchikey_blocks_for_cid : int, optional
        How many blocks of the InChI key to use for the compound id,
        by default 1.
    database : str, optional
        Name of the database where to save results, by default None.
    database_overwrite : bool, optional
        Whether or not to erase existing database in event of a collision,
        by default False.
    mongo_uri : bool, optional
        uri for the mongo client, by default 'mongodb://localhost:27017'.
    image_dir : str, optional
        Filepath where images should be saved, by default None.
    quiet : bool, optional
        Whether to silence warnings, by default False.
    react_targets : bool, optional
        Whether or not to apply reactions to generated compounds that match
        targets, by default True.
    filter_after_final_gen : bool, optional
        Whether to apply filters after final expansion, by default True.
    prune_between_gens : bool, optional
        Whether to prune network between generations if using filters

    Attributes
    ----------
    operators: dict
        Reaction operators to transform compounds with.
    coreactants: dict
        Coreactants required by the operators.
    compounds: dict
        Compounds in the pickaxe network.
    reactions: dict
        Reactions in the pickaxe network.
    generation: int
        The current generation
    explicit_h : bool
        Whether rules utilize explicit hydrogens.
    kekulize : bool
        Whether or not to kekulize compounds before reaction.
    neutralise : bool
        Whether or not to neutralise compounds.
    fragmented_mols : bool
        Whether or not to allow fragmented molecules.
    radical_check : bool
        Whether or not to check and remove radicals.
    image_dir : str, optional
        Filepath where images should be saved.
    errors : bool
        Whether or not to print errors to stdout.
    quiet : bool
        Whether or not to silence warnings.
    filters: List[object]
        A list of filters to apply during the expansion.
    targets : dict
        Molecules to be targeted during expansions.
    target_smiles: List[str]
        The SMILES of all the targets.
    react_targets : bool
        Whether or not to react targets when generated.
    filter_after_final_gen : bool
        Whether or not to filter after the last expansion.
    prune_between_gens : bool, optional
        Whether to prune network between generations if using filters.
    mongo_uri : str
        The connection string to the mongo database.
    cid_num_inchi_blocks : int
        How many blocks of the inchi-blocks to use to generate the compound id.
    """

    def __init__(
        self,
        rule_list: str = None,
        coreactant_list: str = None,
        explicit_h: bool = False,
        kekulize: bool = True,
        neutralise: bool = True,
        errors: bool = True,
        inchikey_blocks_for_cid: int = 1,
        database: str = None,
        database_overwrite: bool = False,
        mongo_uri: bool = "mongodb://localhost:27017",
        image_dir: str = None,
        quiet: bool = True,
        react_targets: bool = True,
        filter_after_final_gen: bool = True,
        prune_between_gens: bool = False,
    ) -> None:
        # Main pickaxe properties
        self.operators = {}
        self.coreactants = {}
        self.compounds = {}
        self.reactions = {}
        self.generation = 0
        # Chemistry options
        self.explicit_h = explicit_h
        self.kekulize = kekulize
        self.neutralise = neutralise
        self.fragmented_mols = False
        self.radical_check = False
        # Other options
        self.image_dir = image_dir
        self.errors = errors
        self.quiet = quiet
        # For filtering
        self.filters = []
        self.targets = dict()
        self.target_smiles = []
        self.react_targets = react_targets
        self.filter_after_final_gen = filter_after_final_gen
        self.prune_between_gens = prune_between_gens
        # database info
        self.mongo_uri = mongo_uri
        # partial_operators
        self.use_partial = False
        self.partial_operators = dict()
        # cid options
        self.cid_num_inchi_blocks = inchikey_blocks_for_cid

        print("----------------------------------------")
        print("Intializing pickaxe object")
        if database:
            # Determine if a specified database is legal
            db = MINE(database, self.mongo_uri)
            if database in db.client.list_database_names():
                if database_overwrite:
                    # If db exists, remove db from all of core compounds
                    # and drop db
                    print(
                        (
                            f"Database {database} already exists. "
                            "Deleting database and removing from core compound"
                            " mines."
                        )
                    )
                    db.core_compounds.update_many({}, {"$pull": {"MINES": database}})
                    db.client.drop_database(database)
                    self.mine = database
                else:
                    print(
                        (
                            f"Warning! Database {database} already exists."
                            "Specify database_overwrite as true to delete "
                            "old database and write new."
                        )
                    )
                    exit("Exiting due to database name collision.")
                    self.mine = None
            else:
                self.mine = database
            del db
        else:
            self.mine = None

        # Use RDLogger to catch errors in log file. SetLevel indicates mode (
        # 0 - debug, 1 - info, 2 - warning, 3 - critical). Default is no errors
        if errors:
            lg.setLevel(0)

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

    def load_targets(
        self,
        target_compound_file: Union[str, None],
        id_field: str = "id",
    ) -> None:
        """Load targets into pickaxe.

        Parameters
        ----------
        target_compound_file : str
            Filepath of target compounds.
        id_field : str, optional
            Header value of compound id in input file, by default 'id'.
        """
        if not target_compound_file:
            print("No target file given")
            return None

        for target_dict in utils.file_to_dict_list(target_compound_file):
            mol = self._mol_from_dict(target_dict)
            if not mol:
                continue
            # Add compound to internal dictionary as a target
            # compound and store SMILES string to be returned
            smi = MolToSmiles(mol, True)
            cpd_name = target_dict[id_field]
            # Only operate on organic compounds
            if "c" in smi.lower():
                SanitizeMol(mol)
                self._add_compound(cpd_name, smi, "Target Compound", mol)
                self.target_smiles.append(smi)

        print(f"{len(self.target_smiles)} target compounds loaded\n")

    def load_compound_set(self, compound_file: str = None, id_field: str = "id") -> str:
        """Load compounds for expansion into pickaxe.

        Parameters
        ----------
        compound_file : str, optional
            Filepath of compounds, by default None.
        id_field : str, optional
            Header value of compound id in input file, by default 'id'.

        Returns
        -------
        str
            List of SMILES that were succesfully loaded into pickaxe.

        Raises
        ------
        ValueError
            No file specified for loading.
        """

        # load compounds
        compound_smiles = []
        if compound_file:
            for cpd_dict in utils.file_to_dict_list(compound_file):
                mol = self._mol_from_dict(cpd_dict)
                if not mol:
                    continue
                # Add compound to internal dictionary as a starting
                # compound and store SMILES string to be returned
                smi = MolToSmiles(mol, True)
                cpd_name = cpd_dict[id_field]
                # Do not operate on inorganic compounds
                if "C" in smi or "c" in smi:
                    SanitizeMol(mol)
                    self._add_compound(
                        cpd_name, smi, cpd_type="Starting Compound", mol=mol
                    )
                    compound_smiles.append(smi)

        else:
            raise ValueError("No input file specified for starting compounds")

        print(f"{len(compound_smiles)} compounds loaded")

        return compound_smiles

    def _load_coreactant(self, coreactant_text: str) -> None:
        """load_coreactants into pickaxe.

        Parameters
        ----------
        coreactant_text : str
            A line read in from the coreactant data file.
        """

        # If coreactant is commented out (with '#') then don't import
        if coreactant_text[0] == "#":
            return
        split_text = coreactant_text.strip().split("\t")
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
        cpd_id = self._add_compound(split_text[0], smi, "Coreactant", mol)
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
        self.coreactants[split_text[0]] = (
            mol,
            cpd_id,
        )

    def _load_operators(self, rule_path: str) -> None:
        """Load reaction rules into pickaxe.

        Parameters
        ----------
        rule_path : str
            Filepath of reaction rules.
        """
        skipped = 0

        # Get the stream for rule input
        if type(rule_path) in [str, Path, PosixPath, WindowsPath]:
            infile = open(rule_path)
        elif type(rule_path) == StringIO:
            infile = rule_path

        with infile:
            # Get all reaction rules from tsv file and store in dict (rdr)
            rdr = csv.DictReader(
                (row for row in infile if not row.startswith("#")), delimiter="\t"
            )
            for rule in rdr:
                try:
                    # Get reactants and products for each reaction into list
                    # form (not ; delimited string)
                    rule["Reactants"] = rule["Reactants"].split(";")
                    rule["Products"] = rule["Products"].split(";")
                    # Ensure that all coreactants are known and accounted for
                    all_rules = rule["Reactants"] + rule["Products"]
                    for coreactant_name in all_rules:
                        if (
                            coreactant_name not in self.coreactants
                            and coreactant_name != "Any"
                        ):
                            raise ValueError(
                                "Undefined coreactant:" f"{coreactant_name}"
                            )
                    # Create ChemicalReaction object from SMARTS string
                    rxn = ReactionFromSmarts(rule["SMARTS"])
                    rule.update(
                        {
                            "_id": rule["Name"],
                            "Reactions_predicted": 0,
                            "SMARTS": rule["SMARTS"],
                        }
                    )
                    # Ensure that we have number of expected reactants for
                    # each rule
                    if rxn.GetNumReactantTemplates() != len(
                        rule["Reactants"]
                    ) or rxn.GetNumProductTemplates() != len(rule["Products"]):
                        skipped += 1
                        print(
                            "The number of coreactants does not match the "
                            "number of compounds in the SMARTS for reaction "
                            "rule: " + rule["Name"]
                        )
                    if rule["Name"] in self.operators:
                        raise ValueError("Duplicate reaction rule name")
                    # Update reaction rules dictionary
                    self.operators[rule["Name"]] = (rxn, rule)
                except Exception as e:
                    raise ValueError(f"{str(e)}\nFailed to parse" f"{rule['Name']}")
        if skipped:
            print(f"WARNING: {skipped} rules skipped")

    def _mol_from_dict(self, input_dict: dict) -> Mol:
        """Generate an RDKit mol object from a dictionary.

        Parameters
        ----------
        input_dict : dict
            Should have "smiles", "inchi", and "structure" fields.

        Returns
        -------
        mol : Mol
            Mol object created from structure string in input_dict.

        Raises
        ------
        ValueError
            Structure field not found in input.
        """
        structure_field = None
        for field in input_dict:
            if str(field).lower() in {"smiles", "inchi", "structure"}:
                structure_field = field
                break

        if not structure_field:
            raise ValueError("Structure field not found in input.")

        if structure_field not in input_dict:
            return
        # Generate Mol object from InChI code if present
        if "InChI=" in input_dict[structure_field]:
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

    def _gen_compound(
        self, cpd_name: str, smi: str, cpd_type: str, mol: Mol = None
    ) -> Tuple[str, dict]:
        """Generate a compound.

        Parameters
        ----------
        cpd_name : str
            Name of the compound to add.
        smi : str
            SMILES of the compound to add.
        cpd_type : str
            Type of compound
        mol : Mol, optional
            RDKit Molecule, by default None.

        Returns
        -------
        Tuple[str, dict]
            Compound id and compound dict.
        """
        cpd_dict = {}
        cpd_id, inchi_key = utils.get_compound_hash(
            smi, cpd_type, self.cid_num_inchi_blocks
        )
        if cpd_id:
            # We don't want to overwrite the same compound from a prior
            # generation so we check with hashed id from above
            if cpd_id not in self.compounds:
                if not mol:
                    mol = MolFromSmiles(smi)
                # expand only Predicted and Starting_compounds
                expand = cpd_type in ["Predicted", "Starting Compound"]
                cpd_dict = {
                    "ID": cpd_name,
                    "_id": cpd_id,
                    "SMILES": smi,
                    "InChI_key": inchi_key,
                    "Type": cpd_type,
                    "Generation": self.generation,
                    "atom_count": utils.get_atom_count(mol, self.radical_check),
                    "Reactant_in": [],
                    "Product_of": [],
                    "Expand": expand,
                    "Formula": CalcMolFormula(mol),
                    "last_tani": 0,
                }
                if cpd_id.startswith("X"):
                    del cpd_dict["Reactant_in"]
                    del cpd_dict["Product_of"]
            else:
                cpd_dict = self.compounds[cpd_id]

            return cpd_id, cpd_dict
        else:
            return None, None

    def _add_compound(
        self, cpd_name: str, smi: str, cpd_type: str, mol: Mol = None
    ) -> str:
        """Generate a compound.

        Parameters
        ----------
        cpd_name : str
            Name of the compound to add.
        smi : str
            SMILES of the compound to add.
        cpd_type : str
            Type of compound.
        mol : Mol, optional
            RDKit Molecule, by default None.

        Returns
        -------
        cpd_id : str
            Compound ID.
        """

        # We don't want to overwrite the same compound from a prior
        # generation so we check with hashed id from above
        cpd_id, cpd_dict = self._gen_compound(cpd_name, smi, cpd_type, mol)
        if cpd_id:
            if (cpd_type == "Target Compound") and (cpd_id not in self.targets):
                self.targets[cpd_id] = cpd_dict
            elif cpd_id not in self.compounds:
                self.compounds[cpd_id] = cpd_dict

                if self.image_dir and self.mine:
                    try:
                        with open(
                            os.path.join(self.image_dir, cpd_id + ".svg"), "w"
                        ) as outfile:

                            mol = MolFromSmiles(cpd_dict["SMILES"])
                            nmol = rdMolDraw2D.PrepareMolForDrawing(mol)
                            d2d = rdMolDraw2D.MolDraw2DSVG(1000, 1000)
                            d2d.DrawMolecule(nmol)
                            d2d.FinishDrawing()
                            outfile.write(d2d.GetDrawingText())
                    except OSError:
                        print(f"Unable to generate image for {cpd_dict['SMILES']}")

        return cpd_id

    def transform_all(self, processes: int = 1, generations: int = 1) -> None:
        """Transform compounds with reaction operators.

        Apply reaction rules to compounds and generate a specified number
        of new generations.

        Parameters
        ----------
        processes : int, optional
            Number of processes to run in parallel, by default 1.
        generations : int, optional
            Number of generations to create, by default 1.
        """

        while self.generation < generations or (
            self.generation == generations and self.filter_after_final_gen
        ):

            for _filter in self.filters:
                _filter.apply_filter(self, processes, generation=generations)

            if self.generation < generations:
                # Prune network to only things that are expanded as white list
                if self.prune_between_gens and self.filters:
                    # Find targets and compounds where expand is true
                    white_list = []
                    for cpd_id, cpd_info in self.compounds.items():
                        if cpd_id.startswith("T"):
                            continue
                        elif cpd_info["Expand"] or cpd_id.startswith("X"):
                            white_list.append(cpd_id)
                        elif f"T{cpd_id[1:]}" in self.targets:
                            white_list.append(cpd_id)

                    # Prune network to only things being expanded
                    self.prune_network(white_list, True)

                print("----------------------------------------")
                print(f"Expanding Generation {self.generation + 1}\n")

                # Starting time for expansion
                time_init = time.time()

                # Tracking compounds formed
                n_comps = len(self.compounds)
                n_rxns = len(self.reactions)

                # Get SMILES to be expanded
                compound_smiles = [
                    cpd["SMILES"]
                    for cpd in self.compounds.values()
                    if cpd["Generation"] == self.generation
                    and cpd["Type"] not in ["Coreactant", "Target Compound"]
                    and cpd["Expand"]
                ]
                # No compounds found
                if not compound_smiles:
                    print(
                        "No compounds to expand in generation "
                        f"{self.generation + 1}. Finished expanding."
                    )
                    return None

                self._transform_helper(compound_smiles, processes)
                self._remove_cofactor_redundancy()

                print(
                    f"Generation {self.generation + 1} finished in"
                    f" {time.time()-time_init} s and contains:"
                )
                print(f"\t\t{len(self.compounds) - n_comps} new compounds")
                print(f"\t\t{len(self.reactions) - n_rxns} new reactions")
                print(f"\nDone expanding Generation: {self.generation + 1}.")
                print("----------------------------------------\n")

            self.generation += 1

    # Partial operator code
    # def load_partial_operators(self, mapped_reactions):
    #     """Generate set of partial operators from a list of mapped reactions
    #     corresponding to the reaction rules being used.
    #     :param mapped_reactions: A .csv file with four columns: rule id,
    #     source, SMARTS, mapping info.
    #     :type mapped_reactions: file
    #     """
    #     # generate partial operators as done in ipynb
    #     if not self.operators:
    #         print("Load reaction rules before loading partial operators")
    #     else:
    #         with open(mapped_reactions) as f:
    #             for line in f.readlines():
    #                 # Grab info from current mapped reaction
    #                 rule, source, smiles, _ = line.strip('\n').split('\t')
    #                 # There should be 2 or more reactants derived from
    #                 # the mapping code The mapped code doesn't include
    #                 # cofactors, so 2 or more means any;any*
    #                 exact_reactants = smiles.split('>>')[0]\
    #                                         .replace(';', '.').split('.')

    #                 base_rule = rule.split('_')[0]
    #                 # base rule must be loaded for partial operator
    #                 # to be used
    #                 if base_rule in self.operators:
    #                     op_reactants = (
    #                       self.operators[base_rule][1]['Reactants']
    #                     )
    #                     if op_reactants.count('Any') >= 2:
    #                         mapped_reactants = []
    #                         for i, r in enumerate(op_reactants):
    #                             if r == 'Any':
    #                                 mapped_reactants.append(
    #                                     exact_reactants.pop(0)
    #                                     )
    #                             else:
    #                                 mapped_reactants.append(r)

    #                         ind_SMARTS = (
    #                           self.operators[base_rule][1]['SMARTS']
    #                         )
    #                         ind_SMARTS = (ind_SMARTS.split('>>')[0].
    #                                       split('>>')[0].replace('(', '').
    #                                       replace(')', '').split('.'))
    #                         # now loop through and generate dictionary entries
    #                         for i, r in enumerate(op_reactants):
    #                             if r != 'Any':
    #                                 pass
    #                             else:
    #                                 # Build entries
    #                                 fixed_reactants = [
    #                                     fr if i != j else 'SMARTS_match'
    #                                     for j, fr in
    #                                     enumerate(mapped_reactants)
    #                                 ]

    # bi_rule =  {
    #     'rule': base_rule,
    #     'rule_reaction': rule,
    #     'reactants': fixed_reactants
    # }
    # if (ind_SMARTS[i] in self.partial_operators:
    #     self.partial_operators[ind_SMARTS[i]].append(bi_rule)
    # else:
    #     self.partial_operators[ind_SMARTS[i]] = [bi_rule]

    # def _filter_partial_operators(self):
    #     # generate the reactions to specifically expand
    #     # based on current compounds
    #     def partial_reactants_exist(partial_rule):
    #         try:
    #             rule_reactants = (
    #               self.operators[partial_rule['rule']][1]['Reactants']
    #             )
    #             cofactor = [False if r == 'Any' else True
    #                         for r in rule_reactants]

    #             reactant_ids = []
    #         for is_cofactor, smi in zip(cofactor, partial_rule['reactants']):
    #                 if is_cofactor:
    #                     reactant_ids.append(self.coreactants[smi][1])
    #                 elif smi == 'SMARTS_match':
    #                     continue
    #                 else:
    #                     reactant_ids.append(utils.get_compound_hash(smi))

    #             reactants_exist = [r in self.compounds for r in reactant_ids]
    #             if all(reactants_exist):
    #                 return True
    #             else:
    #                 return False
    #         except:
    #             return False

    #     filtered_partials = dict()
    #     for SMARTS_match, rules in self.partial_operators.items():
    #         for rule in rules:
    #             if partial_reactants_exist(rule):
    #                 if SMARTS_match in filtered_partials:
    #                     filtered_partials[SMARTS_match].append(rule)
    #                 else:
    #                     filtered_partials[SMARTS_match] = [rule]

    #     return filtered_partials

    def _remove_cofactor_redundancy(self) -> None:
        """Check for reactions that are to be removed.

        Checks for cofactors in rxns that were generated by an any;any rules
        and are specified as generated compounds. Removes redundant reactions
        and ensures cofactors are all labeled as cofactors.
        """
        # Identify compounds who are really cofactors
        cofactors_as_cpds = []
        cofactor_ids = [cofactor[1] for cofactor in self.coreactants.values()]
        for cpd_id in self.compounds:
            if "X" + cpd_id[1:] in cofactor_ids and cpd_id.startswith("C"):
                cofactors_as_cpds.append(cpd_id)

        # Loop through identified compounds and update
        # reactions/compounds accordingly
        rxns_to_del = set()
        for cpd_id in cofactors_as_cpds:
            rxn_ids = set(
                self.compounds[cpd_id]["Product_of"]
                + self.compounds[cpd_id]["Reactant_in"]
            )
            rxns_to_del = rxns_to_del.union(rxn_ids)
            # Check and fix reactions as needed
            for rxn_id in rxn_ids:
                rxn = self.reactions[rxn_id]
                # generate products list with replacements
                reactants = []
                products = []
                for s, reactant in rxn["Reactants"]:
                    if reactant in cofactors_as_cpds:
                        reactants.append((s, self.compounds["X" + reactant[1:]]))
                    else:
                        reactants.append((s, self.compounds[reactant]))

                for s, product in rxn["Products"]:
                    if product in cofactors_as_cpds:
                        products.append((s, self.compounds["X" + product[1:]]))
                    else:
                        products.append((s, self.compounds[product]))

                cofactor_rxn_id, rxn_text = utils.get_reaction_hash(reactants, products)

                # Check to see if reaction makes changes
                # Reactions such as
                # '(1) NCCCNc1cc(C(=O)O)ccc1O + (1) O =>
                #       (1) O + (1) NCCCNc1cc(C(=O)O)ccc1O'
                # are possible. Filter out
                sorted_reactants = sorted(reactants, key=lambda x: x[1]["_id"])
                sorted_products = sorted(products, key=lambda x: x[1]["_id"])

                # Update newly calculated reaction info in self.reactions.
                # Four possibilites:
                # 1. Reaction has same products as reactants -- Remove all
                # 2. Cofactors only as reactants
                # 3. Reaction exists already -- append to existing reaction
                # 4. Reaction doesn't exist -- make new reaction

                # After this, delete redudnant compounds and update
                # other compounds in reaction

                append_new = False
                # Reaction is something like CoF + CoF -> Product + Product
                # Where products are actually just cofactors
                if sorted_reactants == sorted_products:
                    pass  # No need to update anything in self.reactions

                # All reactants are cofactors, don't create reaction
                elif all([r[1]["_id"].startswith("X") for r in sorted_reactants]):
                    pass

                # Reaction already exists, update reaction operators
                elif cofactor_rxn_id in self.reactions:
                    ops = self.reactions[cofactor_rxn_id]["Operators"]
                    ops = ops | self.reactions[rxn_id]["Operators"]
                    self.reactions[cofactor_rxn_id]["Operators"] = ops
                    append_new = True

                # Reaction does not exist, generate new reaction
                else:
                    cofactor_rxn = {
                        "_id": cofactor_rxn_id,
                        # give stoich and id of reactants/products
                        "Reactants": [(s, r["_id"]) for s, r in reactants],
                        "Products": [(s, p["_id"]) for s, p in products],
                        "Operators": rxn["Operators"],
                        "SMILES_rxn": rxn_text,
                    }

                    # Assign reaction to reactions dict
                    self.reactions[cofactor_rxn_id] = cofactor_rxn
                    append_new = True

                for _, cpd in rxn["Reactants"]:
                    if cpd.startswith("C"):
                        if rxn_id in self.compounds[cpd]["Reactant_in"]:
                            self.compounds[cpd]["Reactant_in"].remove(rxn_id)

                        if (
                            append_new
                            and cofactor_rxn_id
                            not in self.compounds[cpd]["Reactant_in"]
                        ):
                            self.compounds[cpd]["Reactant_in"].append(cofactor_rxn_id)

                for _, cpd in rxn["Products"]:
                    if cpd.startswith("C"):
                        if rxn_id in self.compounds[cpd]["Product_of"]:
                            self.compounds[cpd]["Product_of"].remove(rxn_id)

                        if append_new and (
                            cofactor_rxn_id not in self.compounds[cpd]["Product_of"]
                        ):
                            self.compounds[cpd]["Product_of"].append(cofactor_rxn_id)

        if rxns_to_del:
            for rxn_id in rxns_to_del:
                # Remove any links from products/reactants
                # TODO: this should have happened in above loop
                # but getting errors stemming from this
                for _, cpd in self.reactions[rxn_id]["Reactants"]:
                    if cpd.startswith("C"):
                        if rxn_id in self.compounds[cpd]["Reactant_in"]:
                            self.compounds[cpd]["Reactant_in"].remove(rxn_id)

                for _, cpd in self.reactions[rxn_id]["Products"]:
                    if cpd.startswith("C"):
                        if rxn_id in self.compounds[cpd]["Product_of"]:
                            self.compounds[cpd]["Product_of"].remove(rxn_id)

                del self.reactions[rxn_id]

            for cpd_id in cofactors_as_cpds:
                del self.compounds[cpd_id]

        orphan_cpds = []
        # Check for orphaned compounds
        for cpd_id, cpd in self.compounds.items():
            if cpd_id.startswith("C"):
                if (
                    (not cpd["Reactant_in"])
                    and (not cpd["Product_of"])
                    and (cpd["Type"] != "Starting Compound")
                ):
                    orphan_cpds.append(cpd_id)

        for cpd_id in orphan_cpds:
            del self.compounds[cpd_id]

    def prune_network(self, white_list: list, print_output: str = True) -> None:
        """Prune the reaction network to a list of targets.

        Prune the predicted reaction network to only compounds and reactions
        that terminate in a specified white list of compounds.

        Parameters
        ----------
        white_list : list
            A list of compound ids to filter the network to.
        print_output : bool
            Whether or not to print output
        """
        n_white = len(white_list)
        cpd_set, rxn_set = self.find_minimal_set(white_list)
        self.compounds = dict(
            [(k, v) for k, v in self.compounds.items() if k in cpd_set]
        )
        self.reactions = dict(
            [(k, v) for k, v in self.reactions.items() if k in rxn_set]
        )

        if print_output:
            print(
                f"Pruned network to {len(cpd_set)} compounds and "
                f"{len(rxn_set)} reactions based on "
                f"{n_white} whitelisted compounds.\n"
            )

    def prune_network_to_targets(self) -> None:
        """Prune the reaction network to the target compounds.

        Prune the predicted reaction network to only compounds and reactions
        that terminate in the target compounds.
        """
        print("----------------------------------------")
        prune_start = time.time()
        white_list = set()
        for target_id in self.targets:
            white_list.add("C" + target_id[1:])

        print(f"Pruning to {len(white_list)} target compounds")
        self.prune_network(white_list)

        n_targets = 0
        for cpd_id in self.compounds:
            if f"T{cpd_id[1:]}" in self.targets:
                n_targets += 1

        print(f"Found {n_targets} targets.")
        print(f"Pruning took {time.time() - prune_start}s")
        print("----------------------------------------\n")

    def find_minimal_set(self, white_list: Set[str]) -> Tuple[set, set]:
        """Find the minimal set of compounds and reactions given a white list.

        Given a whitelist this function finds the minimal set of compound and
        reactions ids that comprise the set.

        Parameters
        ----------
        white_list : Set[str]
            List of compound_ids to use to filter reaction network to.

        Returns
        -------
        Tuple[set, set]
            The filtered compounds and reactions.
        """

        queue = list(white_list)
        visited = set()
        cpd_set = set()
        rxn_set = set()

        while queue:
            cpd_id = queue.pop()
            visited.add(cpd_id)
            # Select compound from whitelist
            if cpd_id not in self.compounds:
                continue
            else:
                cpd_set.add(cpd_id)

            # Add info for reactions that produce compound
            for rxn_id in self.compounds[cpd_id].get("Product_of", []):
                rxn_set.add(rxn_id)
                # Add products, not to be further explored
                for cpd in self.reactions[rxn_id]["Products"]:
                    cpd_set.add(cpd[1])

                # Add reactants, also add to queue
                for cpd in self.reactions[rxn_id]["Reactants"]:
                    cpd_set.add(cpd[1])
                    # Termination conditions
                    if (
                        cpd[1].startswith("C")
                        and cpd[1] not in visited
                        and self.compounds[cpd[1]]["Type"] != "Starting Compound"
                    ):
                        queue.append(cpd[1])

        return cpd_set, rxn_set

    def assign_ids(self) -> None:
        """Assign a numerical ID to compounds (and reactions).

        Assign IDs that are unique only to the CURRENT run.
        """
        # TODO is this necessary?
        # # If we were running a multiprocess expansion, this removes the dicts
        # # from Manager control
        # self.compounds = dict(self.compounds)
        # self.reactions = dict(self.reactions)
        i = 1
        for comp in sorted(
            self.compounds.values(), key=lambda x: (x["Generation"], x["_id"])
        ):
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            if not comp.get("ID"):
                comp["ID"] = "pkc" + str(i).zfill(7)
                i += 1
                self.compounds[comp["_id"]] = comp
                # If we are not loading into the mine, we generate the image
                # here.
                if self.image_dir and not self.mine:
                    mol = MolFromSmiles(comp["SMILES"])
                    try:
                        MolToFile(
                            mol,
                            os.path.join(self.image_dir, comp["ID"] + ".png"),
                            fitImage=True,
                            kekulize=False,
                        )
                    except OSError:
                        print(f"Unable to generate image for {comp['SMILES']}")
        i = 1
        for rxn in sorted(self.reactions.values(), key=lambda x: x["_id"]):
            rxn["ID_rxn"] = (
                " + ".join(
                    [
                        f"({x[0]}) {self.compounds[x[1]]['ID']}[c0]"
                        for x in rxn["Reactants"]
                    ]
                )
                + " => "
                + " + ".join(
                    [
                        f"({x[0]}) {self.compounds[x[1]]['ID']}[c0]"
                        for x in rxn["Products"]
                    ]
                )
            )
            # Create ID of form ####### ending with i, padded with zeroes to
            # fill unused spots to the left with zfill (e.g. ID = '0003721' if
            # i = 3721).
            rxn["ID"] = "pkr" + str(i).zfill(7)
            i += 1
            self.reactions[rxn["_id"]] = rxn

    def write_compound_output_file(self, path: str, dialect: str = "excel-tab") -> None:
        """Write compounds to an output file.

        Parameters
        ----------
        path : str
            Path to write data.
        dialect : str, optional
            Dialect of the output, by default 'excel-tab'.
        """
        path = utils.prevent_overwrite(path)

        columns = ("ID", "Type", "Generation", "Formula", "InChIKey", "SMILES")
        for _id, val in self.compounds.items():
            inchi_key = MolToInchiKey(MolFromSmiles(val["SMILES"]))
            self.compounds[_id]["InChIKey"] = inchi_key

        with open(path, "w") as outfile:
            writer = csv.DictWriter(
                outfile,
                columns,
                dialect=dialect,
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(sorted(self.compounds.values(), key=lambda x: x["ID"]))

    def write_reaction_output_file(self, path: str, delimiter: str = "\t") -> None:
        r"""Write all reaction data to the specified path.

        Parameters
        ----------
        path : str
            Path to write data.
        delimiter : str, optional
            Delimiter for the output file, by default '\t'.
        """
        path = utils.prevent_overwrite(path)
        with open(path, "w") as outfile:
            outfile.write(
                "ID\tName\tID equation\tSMILES equation\tRxn hash\t" "Reaction rules\n"
            )
            for rxn in sorted(self.reactions.values(), key=lambda x: x["ID"]):
                outfile.write(
                    delimiter.join(
                        [
                            rxn["ID"],
                            "",
                            rxn["ID_rxn"],
                            rxn["SMILES_rxn"],
                            rxn["_id"],
                            ";".join(rxn["Operators"]),
                        ]
                    )
                    + "\n"
                )

    def save_to_mine(
        self, processes: int = 1, indexing: bool = True, write_core: bool = False
    ) -> None:
        """Save pickaxe run to MINE database.

        Parameters
        ----------
        processes : int, optional
            Number of processes to use, by default 1.
        indexing : bool, optional
            Whether or not to add indexes, by default True.
        write_core : bool, optional
            Whether or not to write to core database, by default False.
        """
        print("\n----------------------------------------")
        print(f"Writing results to {self.mine} Database")
        print("----------------------------------------\n")
        start = time.time()
        db = MINE(self.mine, self.mongo_uri)

        # Insert Reactions
        print("--------------- Reactions --------------")
        rxn_start = time.time()
        write_reactions_to_mine(self.reactions.values(), db)
        print(f"Wrote Reactions in {time.time() - rxn_start} seconds.")
        print("----------------------------------------\n")

        # Insert Reactions
        print("--------------- Compounds --------------")
        cpd_start = time.time()
        write_compounds_to_mine(self.compounds.values(), db, processes=processes)
        print(f"Wrote Compounds in {time.time() - cpd_start} seconds.")
        if write_core:
            cpd_start = time.time()
            print("\nWriting Core Compounds")
            write_core_compounds(
                self.compounds.values(), db, self.mine, processes=processes
            )
            print(f"Wrote Core Compounds in {time.time() - cpd_start} seconds.")
        print("----------------------------------------\n")

        if self.targets:
            print("--------------- Targets ----------------")
            target_start = time.time()
            write_targets_to_mine(self.targets.values(), db)
            print(f"Wrote Target Compounds in {time.time() - target_start}" " seconds.")
        else:
            print("No targets to write to MINE.")
        # Save operators
        # Operator values have to be calculated
        if self.operators:
            print("-------------- Operators ---------------")
            operator_start = time.time()
            # update operator rxn count
            for rxn_dict in self.reactions.values():
                for op in rxn_dict["Operators"]:
                    self.operators[op][1]["Reactions_predicted"] += 1

            db.operators.insert_many([op[1] for op in self.operators.values()])
            db.meta_data.insert_one(
                {"Timestamp": datetime.datetime.now(), "Action": "Operators Inserted"}
            )
            print(
                "Done with Operators Overall--took "
                f"{time.time() - operator_start} seconds."
            )
        print("----------------------------------------\n")

        if indexing:
            print("-------------- Indices ---------------")
            index_start = time.time()
            db.build_indexes()
            print(f"Built Indices--took {time.time() - index_start} seconds.")
            print("----------------------------------------\n")

        print("-------------- Overall ---------------")
        print(f"Finished uploading everything in {time.time() - start} sec")
        print("----------------------------------------\n")

    def _transform_helper(self, compound_smiles: List[str], processes: int) -> None:
        """Help to transform reactions external to pickaxe class

        Parameters
        ----------
        compound_smiles : List[str]
            A list of SMILES to expand.
        processes : int
            Number of processes to run.
        """

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
                    rxn_ops = self.reactions[rxn_id]["Operators"]
                    rxn_ops = rxn_ops.union(rxn_dict["Operators"])
                    # if "Partial Operators" in self.reactions[rxn_id]:
                    #     par_ops = self.reactions[rxn_id]["Partial Operators"]
                    #     par_ops = par_ops.union(rxn_dict["Partial Operators"])

                # Update compound tracking
                for product_id in [
                    cpd_id
                    for _, cpd_id in rxn_dict["Products"]
                    if cpd_id.startswith("C")
                ]:
                    if rxn_id not in self.compounds[product_id]["Product_of"]:
                        self.compounds[product_id]["Product_of"].append(rxn_id)

                for reactant_id in [
                    cpd_id
                    for _, cpd_id in rxn_dict["Reactants"]
                    if cpd_id.startswith("C")
                ]:
                    if rxn_id not in self.compounds[reactant_id]["Reactant_in"]:
                        self.compounds[reactant_id]["Reactant_in"].append(rxn_id)

        # to pass coreactants externally
        coreactant_dict = {
            co_key: self.compounds[co_key] for _, co_key in self.coreactants.values()
        }

        new_cpds, new_rxns = transform_all_compounds_with_full(
            compound_smiles,
            self.coreactants,
            coreactant_dict,
            self.operators,
            self.generation + 1,
            self.explicit_h,
            processes,
        )

        update_cpds_rxns(new_cpds, new_rxns)

        # if self.partial_operators:
        #     print("\nGenerating partial operators...")
        #     partial_operators = self._filter_partial_operators()
        #     if partial_operators:
        #         print("Found partial operators, applying.")
        #         # transform partial
        #         new_cpds, new_rxns = _transform_all_compounds_with_partial(
        #                                     compound_smiles, self.coreactants,
        #                                     coreactant_dict, self.operators,
        #                                     self.generation, self.explicit_h,
        #                                     processes, partial_operators
        #                             )

        #         update_cpds_rxns(new_cpds, new_rxns)
        #     else:
        #         print("No partial operators could be generated.")

    def pickle_pickaxe(self, fname: str) -> None:
        """Pickle key pickaxe items.

        Pickle pickaxe object to be loaded in later.

        Parameters
        ----------
        fname : str
            filename to save (must be .pk).
        """
        dict_to_pickle = {
            "compounds": self.compounds,
            "reactions": self.reactions,
            "operators": self.operators,
            "targets": self.targets,
        }

        with open(fname, "wb") as f:
            pickle.dump(dict_to_pickle, f)

    def load_pickled_pickaxe(self, fname: str) -> None:
        """Load pickaxe from pickle.

        Load pickled pickaxe object.

        Parameters
        ----------
        fname : str
            filename to read (must be .pk).
        """
        start_load = time.time()
        print(f"Loading {fname} pickled data.")
        with open(fname, "rb") as f:
            pickle_d = pickle.load(f)
            self.compounds = pickle_d["compounds"]
            self.reactions = pickle_d["reactions"]
            self.operators = pickle_d["operators"]
            self.targets = pickle_d["targets"]
            self.target_smiles = [target["SMILES"] for target in self.targets.values()]

            for key in pickle_d:
                var = getattr(self, key)
                if var:
                    print(f"Loaded {len(var)} {key}")

            print(f"Took {time.time() - start_load}")


if __name__ == "__main__":
    print(os.getcwd())
    # Get initial time to calculate execution time at end
    t1 = time.time()  # pylint: disable=invalid-name
    # Parse all command line arguments
    parser = ArgumentParser()  # pylint: disable=invalid-name
    # Core args
    parser.add_argument(
        "-C",
        "--coreactant_list",
        default="./tests/data/test_coreactants.tsv",
        help="Specify a list of coreactants as a .tsv",
    )
    parser.add_argument(
        "-r",
        "--rule_list",
        default="./tests/data/test_reaction_rules.tsv",
        help="Specify a list of reaction rules as a .tsv",
    )
    parser.add_argument(
        "-c",
        "--compound_file",
        default="./tests/data/test_compounds.tsv",
        help="Specify a list of starting compounds as .tsv or .csv",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Display RDKit errors & warnings",
    )
    # parser.add_argument('--bnice', action='store_true', default=False,
    #                     help="Set several options to enable compatibility "
    #                          "with bnice operators.")
    parser.add_argument(
        "-H",
        "--explicit_h",
        action="store_true",
        default=False,
        help="Specify explicit hydrogen for use in reaction rules.",
    )
    parser.add_argument(
        "-k",
        "--kekulize",
        action="store_true",
        default=True,
        help="Specify whether to kekulize compounds.",
    )
    parser.add_argument(
        "-n",
        "--neutralise",
        action="store_true",
        default=True,
        help="Specify whether to neturalise compounds.",
    )

    parser.add_argument(
        "-m",
        "--processes",
        default=1,
        type=int,
        help="Set the max number of processes.",
    )

    parser.add_argument(
        "-g",
        "--generations",
        default=1,
        type=int,
        help="Set the numbers of time to apply the reaction rules to the compound set.",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=True,
        help="Silence warnings about imbalanced reactions",
    )

    parser.add_argument(
        "-s", "--smiles", default=None, help="Specify a starting compound SMILES."
    )

    # Result args
    parser.add_argument(
        "-p",
        "--pruning_whitelist",
        default=None,
        help="Specify a list of target compounds to prune reaction network down to.",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default=None,
        help="The directory in which to write files.",
    )
    parser.add_argument(
        "-d",
        "--database",
        default=None,
        help="The URI of the database in which to store "
        "output. If not specified, data is written "
        "as tsv files",
    )
    parser.add_argument(
        "-u",
        "--mongo_uri",
        default="mongodb://localhost:27017",
        help="The URI of the mongo database to connect to. Defaults to"
        " mongodb://localhost:27017",
    )
    parser.add_argument(
        "-i",
        "--image_dir",
        default=None,
        help="Specify a directory to store images of all created compounds",
    )
    parser.add_argument(
        "-wc",
        "--write_core",
        default=False,
        help="Whether or not to write results into core database.",
    )

    OPTIONS = parser.parse_args()

    if not any([OPTIONS.database, OPTIONS.output_dir]):
        exit("No output selected, terminating run.")

    pk = Pickaxe(
        coreactant_list=OPTIONS.coreactant_list,
        rule_list=OPTIONS.rule_list,
        errors=OPTIONS.verbose,
        explicit_h=OPTIONS.explicit_h,
        kekulize=OPTIONS.kekulize,
        neutralise=OPTIONS.neutralise,
        image_dir=OPTIONS.image_dir,
        quiet=OPTIONS.quiet,
        database=OPTIONS.database,
        mongo_uri=OPTIONS.mongo_uri,
    )
    # Create a directory for image output file if it doesn't already exist
    if OPTIONS.image_dir and not os.path.exists(OPTIONS.image_dir):
        os.mkdir(OPTIONS.image_dir)
    # If starting compound specified as SMILES string, then add it
    if OPTIONS.smiles:
        pk._add_compound("Start", OPTIONS.smiles, cpd_type="Starting Compound")
    else:
        pk.load_compound_set(compound_file=OPTIONS.compound_file)
    # Generate reaction network
    pk.transform_all(processes=OPTIONS.processes, generations=OPTIONS.generations)
    if OPTIONS.pruning_whitelist:
        mols = [
            pk._mol_from_dict(line)
            for line in utils.file_to_dict_list(OPTIONS.pruning_whitelist)
        ]
        pk.prune_network([utils.get_compound_hash(x) for x in mols if x])

    if OPTIONS.output_dir:
        pk.assign_ids()
        pk.write_compound_output_file(OPTIONS.output_dir + "/compounds.tsv")
        pk.write_reaction_output_file(OPTIONS.output_dir + "/reactions.tsv")

    if OPTIONS.database:
        pk.save_to_mine(processes=OPTIONS.processes, write_core=OPTIONS.write_core)

    print(f"Execution took {time.time() - t1} seconds.")
