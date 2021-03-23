"""Compound_io.py: Functions to load MINE databases from and dump compounds
into common cheminformatics formats"""
import collections
import csv
import datetime
import os
import sys
from typing import List, Tuple, Union

from rdkit.Chem import AllChem

from minedatabase import utils
from minedatabase.databases import MINE


def export_sdf(mine_db: MINE, dir_path: str, max_compounds: int = None) -> None:
    """Exports compounds from the database as an MDL SDF file.

    Parameters
    ----------
    mine_db : MINE
        MINE object that contains the database.
    dir_path : str
        Directory for files.
    max_compounds : int, optional
        Maximum number of compounds per file, by default None.
    """

    # Make sure that all compounds point to all their reactants
    if not mine_db.compounds.find_one({"Product_of": {"$exists": 1}}):
        mine_db.add_rxn_pointers()

    print(
        f"Exporting {mine_db.compounds.count()} compounds from {mine_db.name}"
        " as an SDF file"
    )
    target = utils.prevent_overwrite(os.path.join(dir_path, mine_db.name) + "_1.sdf")
    # SDWriter (rdkit) writes Mol objects to SD files
    writer = AllChem.SDWriter(target)
    writer.SetKekulize(True)
    n_files = 1
    for compound in mine_db.compounds.find():
        # Convert SMILES string to Mol object, replacing 'CoA' and 'R' by '*'
        mol = AllChem.MolFromSmiles(compound["SMILES"], True, {"CoA": "*", "R": "*"})
        # if Mol object successfully generated, annotate properties
        if mol:
            mol.SetProp("_id", compound["_id"])
            mol.SetProp("Generation", str(compound["Generation"]))
            if "Reactant_in" in compound:
                mol.SetProp("Reactant_in", str(compound["Reactant_in"]))
            if "Product_of" in compound:
                mol.SetProp("Product_of", str(compound["Product_of"]))
            writer.write(mol)
            # Start writing a new sdf file if the maximum (set by user) has
            # been reached for the current file
            if max_compounds and (writer.NumMols() >= max_compounds):
                n_files += 1
                target = utils.prevent_overwrite(
                    os.path.join(dir_path, mine_db.name) + f"_(n_files).sdf"
                )
                writer = AllChem.SmilesWriter(target)
    writer.close()


def export_smiles(mine_db: MINE, dir_path: str, max_compounds: int = None) -> None:
    """Exports compounds from the database as a SMILES file.

    Parameters
    ----------
    mine_db : MINE
        MINE object that contains the database.
    dir_path : str
        Directory for files.
    max_compounds : int, optional
        Maximum number of compounds per file, by default None.
    """
    header = ["SMILES", "_id", "Generation", "Reactant_in", "Product_of"]
    # Make sure that all compounds point to all their reactants
    if not mine_db.compounds.find_one({"Product_of": {"$exists": 1}}):
        mine_db.add_rxn_pointers()

    print(
        f"Exporting {mine_db.compounds.count()} compounds from {mine_db.name()}"
        " as SMILES file"
    )
    target = open(
        utils.prevent_overwrite(os.path.join(dir_path, mine_db.name) + "_1.smiles"), "w"
    )

    # DictWriter allows for each key:value pair of a dictionary to be written
    # on its own row (by writerow)
    writer = csv.DictWriter(target, fieldnames=header, dialect="excel-tab")
    n_files = 1
    i = 0
    for compound in mine_db.compounds.find({}, dict([(x, 1) for x in header])):
        writer.writerow(compound)
        i += 1
        # If max compounds per file has been set by user and our number of
        # compounds that we have written so far is divisible by the max number,
        # then we start a new file
        if max_compounds and not i % max_compounds:
            n_files += 1
            target = open(
                utils.prevent_overwrite(
                    os.path.join(dir_path, mine_db.name) + f"_{n_files}.smiles"
                ),
                "w",
            )
            writer = csv.DictWriter(target, fieldnames=header, dialect="excel-tab")


def export_mol(mine_db: MINE, target: str, name_field: str = "_id") -> None:
    """Exports compounds from the database as a MDL molfiles

    Parameters
    ----------
    mine_db : MINE
        MINE object that contains the database.
    target : str
        Directory in which to place the files.
    name_field : str, optional
        FIeld to provide names for the mol files. Must be unique and universal.
        By default, "_id".
    """
    # Create the file if it doesn't yet exist
    if not os.path.exists(target):
        os.mkdir(target)

    # Let user know if an id does not exist for every compound in database
    if (
        mine_db.compounds.find().count()
        != mine_db.compounds.find({name_field: {"$exists": 1}}).count()
    ):
        raise ValueError(
            f"{name_field} does not exist for every compound in the database"
        )

    for compound in mine_db.compounds.find({"_id": {"$regex": "^C"}}):
        # Create Mol object from SMILES code for each compound using
        # MolFromSmiles (rdkit). Take stereochemistry into account (True),
        # and replace CoA and R with *.
        mol = AllChem.MolFromSmiles(compound["SMILES"], True, {"CoA": "*", "R": "*"})
        if "." in name_field:
            compound[name_field] = utils.get_dotted_field(compound, name_field)
        # Make things more compact and look nicer
        if isinstance(compound[name_field], list):
            compound[name_field] = ",".join(compound[name_field])
        # Use MolToMolFile (rdkit) to create a mol file from the Mol object
        # with the file path specified.
        AllChem.MolToMolFile(mol, os.path.join(target, compound[name_field] + ".mol"))


def export_tsv(
    mine_db: MINE,
    target: str,
    compound_fields: Tuple[str] = (
        "_id",
        "Names",
        "Model_SEED",
        "Formula",
        "Charge",
        "Inchi",
    ),
    reaction_fields: Tuple[str] = ("_id", "SMILES_rxn", "C_id_rxn"),
) -> None:
    """Exports MINE compound and reaction data as tab-separated values files
    amenable to use in ModelSEED.

    Parameters
    ----------
    mine_db : MINE
        The database to export.
    target : str
        Directory, in which to place the files.
    compound_fields : Tuple[str], optional
        Fields to export in the compound table, by default
        ('_id', 'Names', 'Model_SEED', 'Formula', 'Charge', 'Inchi').
    reaction_fields : Tuple[str], optional
        Fields to export in the reaction table, by default
        ('_id', 'SMILES_rxn', 'C_id_rxn').
    """
    db_links = ("KEGG", "Model_SEED", "PubChem")
    print(f"Exporting {mine_db.compounds.count()} compounds from {mine_db.name} to tsv")
    with open(
        utils.prevent_overwrite(os.path.join(target, mine_db.name) + "_compounds.tsv"),
        "w",
    ) as out:
        writer = csv.DictWriter(out, fieldnames=compound_fields, dialect="excel-tab")
        writer.writeheader()
        for compound in mine_db.compounds.find(
            {},
            dict(
                [("SMILES", 1)]
                + [
                    ("DB_links." + x, 1) if x in db_links else (x, 1)
                    for x in compound_fields
                ]
            ),
        ):
            # This is a work around for supporting older MINEs which lack Inchi
            if "Inchi" in compound_fields and "Inchi" not in compound:
                compound["Inchi"] = AllChem.MolToInchi(
                    AllChem.MolFromSmiles(compound["SMILES"])
                )
            if "SMILES" not in compound_fields:
                del compound["SMILES"]
            if "DB_links" in compound:
                for k, v in compound["DB_links"].items():
                    compound[k] = ", ".join(v)
                del compound["DB_links"]
            writer.writerow(compound)

    print(f"Exporting {mine_db.reactions.count()} reactions from {mine_db.name} to tsv")
    with open(
        utils.prevent_overwrite(os.path.join(target, mine_db.name) + "_reactions.tsv"),
        "w",
    ) as out:
        writer = csv.DictWriter(out, fieldnames=reaction_fields, dialect="excel-tab")
        writer.writeheader()
        for rxn in mine_db.reactions.find(
            {},
            dict(
                [("Reactants", 1), ("Products", 1)] + [(x, 1) for x in reaction_fields]
            ),
        ):
            if "C_id_rxn" in reaction_fields:

                def to_str(half_rxn):
                    return [f"({x['stoich']}) {x['c_id']}" for x in half_rxn]

                rxn["C_id_rxn"] = (
                    " + ".join(to_str(rxn["Reactants"]))
                    + " => "
                    + " + ".join(to_str(rxn["Products"]))
                )
            if "Reactants" not in reaction_fields:
                del rxn["Reactants"]
            if "Products" not in reaction_fields:
                del rxn["Products"]
            writer.writerow(rxn)


def export_kbase(mine_db: MINE, target: str) -> None:
    """Exports MINE compound and reaction data as tab-separated values files
    amenable to use in ModelSEED.

    Parameters
    ----------
    mine_db : MINE
        The database to export.
    target : str
        Directory in which to place the files.
    """
    compound_fields = collections.OrderedDict(
        [
            ("id", "_id"),
            ("name", ""),
            ("formula", "Formula"),
            ("charge", "Charge"),
            ("aliases", "Names"),
        ]
    )
    reaction_fields = collections.OrderedDict(
        [
            ("id", "_id"),
            ("direction", ">"),
            ("compartment", "c0"),
            ("gpr", ""),
            ("name", ""),
            ("enzyme", ""),
            ("pathway", ""),
            ("reference", ""),
            ("equation", ""),
        ]
    )
    print(
        f"Exporting {mine_db.compounds.count()} compounds from {mine_db.name()} to tsv"
    )
    with open(
        utils.prevent_overwrite(os.path.join(target, mine_db.name) + "_compounds.tsv"),
        "w",
    ) as out:
        writer = csv.DictWriter(out, fieldnames=compound_fields, dialect="excel-tab")
        writer.writeheader()
        for compound in mine_db.compounds.find(
            {},
            dict(
                [("Names", 1), ("DB_links.Model_SEED", 1)]
                + [(x, 1) for x in compound_fields.values()]
            ),
        ):
            if compound["_id"][0] == "X":
                continue
            for k, v in compound_fields.items():
                if v in compound:
                    compound[k] = compound[v]
                    del compound[v]
            if "name" in compound_fields and "Names" in compound:
                compound["name"] = compound["Names"][0]
                del compound["Names"]
            if "aliases" in compound:
                compound["aliases"] = "|".join(compound["aliases"])
                if "Model_SEED" in compound["DB_links"]:
                    compound["aliases"] += "|" + "|".join(
                        sorted(compound["DB_links"]["Model_SEED"])
                    )
            if "DB_links" in compound:
                del compound["DB_links"]
            writer.writerow(compound)

    print(f"Exporting {mine_db.reactions.count()} reactions from {mine_db.name} to tsv")
    with open(
        utils.prevent_overwrite(os.path.join(target, mine_db.name) + "_reactions.tsv"),
        "w",
    ) as out:
        writer = csv.DictWriter(out, fieldnames=reaction_fields, dialect="excel-tab")
        writer.writeheader()
        for rxn in mine_db.reactions.find(
            {},
            dict(
                [("Reactants", 1), ("Products", 1)]
                + [(x, 1) for x in reaction_fields.values()]
            ),
        ):
            for k, v in reaction_fields.items():
                if v in rxn:
                    rxn[k] = rxn[v]
                    del rxn[v]
            if "equation" in reaction_fields:

                def to_str(half_rxn):
                    return [
                        f"({x['stoich']}) {x['c_id'].replace('X', 'C')}"
                        for x in half_rxn
                    ]

                rxn["equation"] = (
                    " + ".join(to_str(rxn["Reactants"]))
                    + " => "
                    + " + ".join(to_str(rxn["Products"]))
                )
            if "Reactants" not in reaction_fields:
                del rxn["Reactants"]
            if "Products" not in reaction_fields:
                del rxn["Products"]
            writer.writerow(rxn)


def export_inchi_rxns(
    mine_db: MINE, target: str, rxn_ids: Union[List[str], None] = None
) -> None:
    """Export reactions from a MINE db to a .tsv file.

    Parameters
    ----------
    mine_db : MINE
        Name of MongoDB to export reactions from.
    target : str
        Path to folder to save .tsv export file in.
    rxn_ids : Union[List[str], None], optional
        Only export reactions with these ids, by default None.
    """
    reaction_fields = collections.OrderedDict(
        [("Reaction Rule", "Operators"), ("ID", "_id"), ("Equation", "")]
    )
    comp_memo = {}

    def get_name_and_inchi(comp_id):
        if comp_id not in comp_memo:
            comp = mine_db.compounds.find_one(
                {"_id": comp_id}, {"Names": 1, "Inchi": 1, "MINE_id": 1}
            )
            comp_memo[comp_id] = (
                comp.get("Names", [comp["MINE_id"]])[0],
                comp.get("Inchi"),
            )
        return comp_memo[comp_id]

    def to_str(half_rxn):
        lst = []
        for x in half_rxn:
            name, inchi = get_name_and_inchi(x["c_id"])
            lst.append(f"({x['stoich']}) {name}[{inchi}]")
        return lst

    with open(
        utils.prevent_overwrite(os.path.join(target, mine_db.name) + "_reactions.tsv"),
        "w",
    ) as out:
        writer = csv.DictWriter(out, fieldnames=reaction_fields, dialect="excel-tab")
        writer.writeheader()
        if rxn_ids:
            query = {"_id": {"$in": rxn_ids}}
        else:
            query = {}
        for rxn in mine_db.reactions.find(
            query,
            dict(
                [("Reactants", 1), ("Products", 1)]
                + [(x, 1) for x in reaction_fields.values()]
            ),
        ):
            for k, v in reaction_fields.items():
                if v in rxn:
                    if isinstance(rxn[v], list):
                        rxn[k] = ", ".join(rxn[v])
                    else:
                        rxn[k] = rxn[v]
                    del rxn[v]
            if "Equation" in reaction_fields:
                rxn["Equation"] = (
                    " + ".join(to_str(rxn["Reactants"]))
                    + " => "
                    + " + ".join(to_str(rxn["Products"]))
                )
            if "Reactants" not in reaction_fields:
                del rxn["Reactants"]
            if "Products" not in reaction_fields:
                del rxn["Products"]
            writer.writerow(rxn)


def import_sdf(mine_db: MINE, target: str) -> None:
    """Imports a SDF file as a MINE database.

    Parameters
    ----------
    mine_db : MINE
        The database to export.
    target : str
        Directory in which to place the files.
    """
    # SDMolSupplier (rdkit) takes entries from sdf file and returns Mol objects
    sdf_gen = AllChem.SDMolSupplier(target)
    # Go through each generated Mol object and add each to MINE database
    for mol in sdf_gen:
        mine_db.insert_compound(
            mol,
            compound_dict=mol.GetPropsAsDict(),
            pubchem_db=None,
            kegg_db=None,
            modelseed_db=None,
        )
    # Add to log file (metadata)
    mine_db.meta_data.insert(
        {
            "Timestamp": datetime.datetime.now(),
            "Action": "SDF Imported",
            "Filepath": target,
        }
    )


def import_smiles(mine_db: MINE, target: str) -> None:
    """Imports a smiles file as a MINE database.

    Parameters
    ----------
    mine_db : MINE
        The database to export.
    target : str
        Directory in which to place the files.
    """
    # SmilesMolSupplier (rdkit) generates Mol objects from smiles file (.smi)
    mols = AllChem.SmilesMolSupplier(target, delimiter="\t", nameColumn=0)
    # Go through each generated mol file and add molecule to MINE database
    # Stores compound properties in dict (GetPropsAsDict() from rdkit Mol
    # class)
    for mol in mols:
        if mol:
            mine_db.insert_compound(
                mol,
                compound_dict=mol.GetPropsAsDict(),
                pubchem_db=None,
                kegg_db=None,
                modelseed_db=None,
            )
    # Add to log file (metadata)
    mine_db.meta_data.insert(
        {
            "Timestamp": datetime.datetime.now(),
            "Action": "SDF Imported",
            "Filepath": target,
        }
    )


def import_mol_dir(
    mine_db: MINE, target: str, name_field: str = "Name", overwrite: bool = False
) -> None:
    """Imports a directory of molfiles as a MINE database.

    Parameters
    ----------
    mine_db : MINE
        The database to export.
    target : str
        Directory in which to place the files.
    name_field : str, optional
        Field for the compound name, by default "Name".
    overwrite : bool, optional
        Replace old compounds with new ones if a collision happens, by default False.
    """
    # For each .mol file in the directory of the target folder (path):
    for file in os.listdir(target):
        if ".mol" in file:
            # MolFromMolFile (rdkit) generates Mol objects from .mol files
            mol = AllChem.MolFromMolFile(target + "/" + file)
            # Mol object name becomes name of mol file without .mol extension
            name = file.rstrip(".mol")
            # Check that Mol object is successfully generated
            if mol:
                # Create hashkey for the compound
                cpdhash = utils.get_compound_hash(mol)
                # If we don't want to overwrite, and the compound (cpdhash)
                # already exists, then add an extra cpdhash for that molecule
                if not overwrite and mine_db.compounds.count({"_id": cpdhash}):
                    mine_db.compounds.update(
                        {"_id": cpdhash}, {"$addToSet": {name_field: name}}
                    )
                # If we don't care about overwriting, just insert the new
                # compound into the database
                else:
                    mine_db.insert_compound(
                        mol,
                        compound_dict={name_field: [name], "Generation": 0},
                        pubchem_db=None,
                        kegg_db=None,
                        modelseed_db=None,
                    )
    # Add to log file (metadata)
    mine_db.meta_data.insert(
        {
            "Timestamp": datetime.datetime.now(),
            "Action": "MolFiles Imported",
            "Filepath": target,
        }
    )


if __name__ == "__main__":
    # User inputs task as first argument (export-sdf, export-smi, export-mol,
    #  import-sdf, import-smi, or import-mol)
    TASK = sys.argv[1]
    # User inputs database name as second argument
    DB_NAME = sys.argv[2]
    # User inputs file path as third argument
    PATH = sys.argv[3]
    database = MINE(DB_NAME)  # pylint: disable=invalid-name
    if TASK == "export-sdf":
        # If a maximum molecules per file is specified (fourth argument
        # entered by user), then pass that to the export function.
        if len(sys.argv) == 5:
            export_sdf(database, PATH, int(sys.argv[4]))
        # Otherwise, assume an unlimited number of molecules per file
        else:
            export_sdf(database, PATH)
    elif TASK == "export-smi":
        # If a maximum molecules per file is specified (fourth argument
        # entered by user), then pass that to the export function.
        if len(sys.argv) == 5:
            export_smiles(database, PATH, int(sys.argv[4]))
        # Otherwise, assume an unlimited number of molecules per file
        else:
            export_smiles(database, PATH)
    elif TASK == "export-mol":
        # If a maximum molecules per file is specified (fourth argument
        # entered by user), then pass that to the export function.
        if len(sys.argv) == 5:
            export_mol(database, PATH, sys.argv[4])
        # Otherwise, assume an unlimited number of molecules per file
        else:
            export_mol(database, PATH)
    elif TASK == "export-tsv":
        export_tsv(database, PATH)
    elif TASK == "export-kbase":
        export_kbase(database, PATH)
    elif TASK == "import-sdf":
        import_sdf(database, PATH)
    elif TASK == "import-smi":
        import_smiles(database, PATH)
    elif TASK == "import-mol":
        import_mol_dir(database, PATH)
    else:
        print("ERROR: Unrecognised TASK")
