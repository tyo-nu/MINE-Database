"""A script to generate a metanetx database.

The purpose of the metanetx database is to provide mapping
from InChI keys to a number of database identifiers. This database
will then populate the website if there is an inchi match.

Running this script requires downloading the following from
https://www.metanetx.org/mnxdoc/mnxref.html

   1. chem_xref.tsv
   2. chem_prop.tsv

The data version used is based on the description in the header:
#Based on the following resources:
#
#RESOURCE:  MetaNetX/MNXref
#VERSION:   4.1
#DATE:      2020/09/17
#URL:       https://www.metanetx.org
"""
from pathlib import Path

import pandas as pd

from collections import defaultdict
import pymongo

pwd = Path(__file__)
pwd = pwd.parent
METANETX_PATH = (pwd / "../local_data/metanetx").resolve()

def get_cross_references(row):
    current_reference = {}

    if ":" in row["#source"]:
        current_reference["source"] = row["#source"].split(":")[0]
        current_reference["source_id"] = row["#source"].split(":")[1]
    else:
        current_reference["source"] = row["#source"]
        current_reference["source_id"] = row["#source"]
    current_reference["description"] = (
        row["description"] if not pd.isna(row["description"])
        else None
    )
    cross_ref_dict[row.ID].append(current_reference)


def get_db_entry(row):
    dict_for_db[row["#ID"]] = {
            "mnxm_id": row["#ID"],
            "Inchikey": row.InChIKey,
            "primary_reference": row.reference,
            "cross_references": cross_ref_dict[row["#ID"]]
    }


if __name__ == "__main__":
    # First step: Generate panda dfs of the xref and props
    skiprows = 347
    chem_prop_df = pd.read_csv(
        METANETX_PATH / "chem_prop.tsv",
        delimiter="\t",
        skiprows=skiprows
    )
    chem_prop_df = chem_prop_df[~chem_prop_df["InChIKey"].isna()]
    chem_prop_df = chem_prop_df[~chem_prop_df["formula"].isna()]

    chem_xref_df = pd.read_csv(
        METANETX_PATH / "chem_xref.tsv",
        delimiter="\t",
        skiprows=skiprows
    )

    # Map functions on pandas dataframes to populate dictionaries
    cross_ref_dict = defaultdict(list)
    dict_for_db = dict()

    chem_xref_df.apply(get_cross_references, axis=1)
    chem_prop_df.apply(get_db_entry, axis=1)

    print("Inserting into Mongo.")
    mongo_uri = open(pwd / "../mongo_uri.csv").readline().strip("\n")
    client = pymongo.MongoClient(mongo_uri)
    client.compound_references.data.insert_many(dict_for_db.values(), ordered=False)
    client.compound_references.data.create_index("Inchikey")
