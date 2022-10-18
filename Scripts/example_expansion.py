"""Template for a pickaxe run.

This python script provides a skeleton to build Pickaxe runs around
The general format of a script will be:
   1. Connect to mongoDB (if desired)
   2. Load reaction rules and cofactors
   3. Load starting compounds
   4. Load Filtering Options
   5. Transform Compounds
   6. Write results
"""

import datetime
import multiprocessing
import pickle
import time

import pymongo

# Uncomment to use these. Pickaxe doesn't come packaged with dependencies by default.
# from minedatabase.filters import ThermoFilter
# from minedatabase.filters import ReactionFeasibilityFilter

from minedatabase.pickaxe import Pickaxe
from minedatabase.rules import metacyc_generalized, metacyc_intermediate

start = time.time()

# Writing compound and reaction csv files locally
write_to_csv = True
output_dir = "."

# Input compounds in a csv folder with headings:
# id,smiles
input_cpds = "./example_data/starting_cpds_single.csv"

# Rule specification and generation. Rules can be manually created or
# metacyc_intermediate or metacyc_generalized can provide correctly formatted
# biological reactions derived from metacyc.
#
# See the documentation for description of options.
rule_list, coreactant_list, rule_name = metacyc_intermediate(
    fraction_coverage=0.2
)

# Core Pickaxe Run Options
generations = 1              # Total rounds of rule applications
processes = 1                # Number of processes for parallelization
verbose = False              # Display RDKit warnings and errors


if __name__ == "__main__":
    # Use "spawn" for multiprocessing
    multiprocessing.set_start_method("spawn")

    pk = Pickaxe(
        coreactant_list=coreactant_list,
        rule_list=rule_list,
        errors=verbose,
    )

    # Load compounds
    pk.load_compound_set(compound_file=input_cpds)

    # Transform compounds (the main step)
    pk.transform_all(processes, generations)
 
    if write_to_csv:
        pk.assign_ids()
        pk.write_compound_output_file(output_dir + "/compounds.tsv")
        pk.write_reaction_output_file(output_dir + "/reactions.tsv")
