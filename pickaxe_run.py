"""pickaxe_run.py
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
import time

import pymongo

from minedatabase.filters import (MCSFilter, MetabolomicsFilter,
                                  TanimotoFilter, TanimotoSamplingFilter)
from minedatabase.pickaxe import Pickaxe

# pylint: disable=invalid-name

start = time.time()

###############################################################################
##### Database and output information
# The default mongo is localhost:27017
# Connecting remotely requires the location of the database
# as well as username/password if security is being used.
# Username/password are stored in credentials.csv
# in the following format: username

# Database to write results to
write_db = True
database_overwrite = True
# database = "APAH_100Sam_50rule"
database = "test_ADP1"
# Message to insert into metadata
message = ("Example run to show how pickaxe is ran.")

# mongo DB information
use_local = True
if use_local:
    mongo_uri = 'mongodb://localhost:27017'
else:
    mongo_uri = open('mongo_uri.csv').readline().strip('\n')

# Write output .csv files locally
write_local = False
output_dir = '.'
###############################################################################

###############################################################################
#    Starting Compounds, Cofactors, and Rules
# Input compounds
input_cpds = './example_data/starting_cpds_ten.csv'

# Rules from Joseph Ni
coreactant_list = './minedatabase/data/MetaCyc_Coreactants.tsv'
# rule_list = './minedatabase/data/intermediate_rules_uniprot.tsv'
rule_list = './minedatabase/data/metacyc_272rules_90percentMapping.tsv'
#rule_list = '../Operator_Filter/all_mapped_rules.tsv'

# Input compounds
input_cpds = './local_data/ADP1_cpds_out_reduced.csv'

###############################################################################

###############################################################################
# Core Pickaxe Run Options
generations = 1
processes = 1     # Number of processes for parallelization
inchikey_blocks_for_cid = 1 # Number of blocks of the inchi key to use for the compound id
verbose = False     # Display RDKit warnings and errors
explicit_h = False
kekulize = True
neutralise = True
image_dir = None
quiet = True
indexing = False
###############################################################################

###############################################################################
##### All Filter and Sampler Options
###############################################################################
# Global Filtering Options

# Path to target cpds file (not required for metabolomics filter)
target_cpds = './local_data/ADP1_cpds_out_final.csv'

# Should targets be flagged for reaction
react_targets = False

# Prune results to remove compounds not required to produce targets
prune_to_targets = False

# Filter final generation?
filter_after_final_gen = True

##############################################################################
# Tanimoto Filtering options.
# Filters by tanimoto similarity score, using default RDKit fingerprints

# Apply this filter?
tani_filter = False

# Tanimito filter threshold. Can be single number or a list with length at least
# equal to the number of generations (+1 if filtering after expansion)
tani_threshold = [0, 0.2, 0.7]

# Make sure tani increases each generation?
increasing_tani = False

###############################################################################
# Tanimoto-based Sampling Options
# Samples by tanimoto similarity score, using default RDKit fingerprints

# Apply this sampler?
tani_sample = True

# Number of compounds per generation to sample
sample_size = 100

# weight is a function that specifies weighting of Tanimoto similarity
# weight accepts one input
# T : float in range 0-1
# and returns
# float in any range (will be rescaled later)
# weight = None will use a T^4 to weight.
def weight(T):
    return T**4

# How to represent the function in text
weight_representation = "T^4"

###############################################################################
# Maximum common substructure (MCS) filter

# Apply this filter?
mcs_filter = False

# Finds the MCS of the target and compound and identifies fraction of target the MCS composes
crit_mcs = [0.3, 0.8, 0.95]

##############################################################################
# Metabolomics Filter Options

# Apply this filter?
metabolomics_filter = True

# Path to csv with list of detected masses. For example:
# Peak ID, Retention Time, Aggregate M/Z, Polarity, Compound Name, Predicted Structure (smile), ID
# Peak1, 6.33, 74.0373, negative, propionic acid, CCC(=O)O, yes
# Peak2, 26.31, 84.06869909, positive, , , no
# ...
met_data_path = './local_data/ADP1_Metabolomics_PeakList_final.csv'

# Name of dataset
met_data_name = 'ADP1_metabolomics'

# Adducts to add to each mass in mass list to create final list of possible masses.
# See "./minedatabase/data/adducts/All adducts.txt" for options.
possible_adducts = ['[M+H]+', '[M-H]-']

# Tolerance in Da
mass_tolerance = 0.001

###############################################################################

###############################################################################
# Verbose output
print_parameters = True

def print_run_parameters():
    def print_parameter_list(plist):
        for i in plist:
            print(f"--{i}: {eval(i)}")

    print('\n-------------Run Parameters-------------')

    print ('\nRun Info')
    print_parameter_list(['coreactant_list', 'rule_list', 'input_cpds'])

    print('\nExpansion Options')
    print_parameter_list(['generations', 'processes'])

    print('\nGeneral Filter Options')
    print_parameter_list(['filter_after_final_gen', 'react_targets',
                    'prune_to_targets', 'react_targets'])
    
    if tani_sample:
        print('\nTanimoto Sampling Filter Options')
        print_parameter_list(['sample_size', 'weight_representation'])
    
    if tani_filter:
        print('\nTanimoto Threshold Filter Options')
        print_parameter_list(['tani_threshold', 'increasing_tani'])

    if mcs_filter:
        print('\nMaximum Common Substructure Filter Options')
        print_parameter_list(['crit_mcs'])

    if metabolomics_filter:
        print('\nMetabolomics Filter Options')
        print_parameter_list(['met_data_path', 'met_data_name',
                              'possible_adducts', 'mass_tolerance'])
    
    print('\nPickaxe Options')
    print_parameter_list(['verbose', 'explicit_h', 'kekulize', 'neutralise',
                'image_dir', 'quiet', 'indexing'])
    print('----------------------------------------\n')
###############################################################################

###############################################################################
##### Running pickaxe
if __name__ == '__main__':  # required for parallelization on Windows
    # Use 'spawn' for multiprocessing
    multiprocessing.set_start_method('spawn')
    # Initialize the Pickaxe class
    if write_db is False:
        database = None

    if print_parameters:
        print_run_parameters()

    pk = Pickaxe(coreactant_list=coreactant_list, rule_list=rule_list,
                 errors=verbose, explicit_h=explicit_h, kekulize=kekulize,
                 neutralise=neutralise, image_dir=image_dir,
                 inchikey_blocks_for_cid=inchikey_blocks_for_cid, database=database,
                 database_overwrite=database_overwrite, mongo_uri=mongo_uri,
                 quiet=quiet, react_targets=react_targets,
                 filter_after_final_gen=filter_after_final_gen)

    # Load compounds
    pk.load_compound_set(compound_file=input_cpds)

    # Load partial operators
    # if partial_rules:
    #     pk.load_partial_operators(mapped_rxns)

    # Load target compounds for filters
    #if (tani_filter or mcs_filter or tani_sample):
    pk.load_targets(target_cpds, structure_field="SMILES")

    # Apply filters
    if tani_filter:
        taniFilter = TanimotoFilter(filter_name="Tani", crit_tani=tani_threshold,
                                    increasing_tani=increasing_tani)
        pk.filters.append(taniFilter)

    if tani_sample:
        taniSampleFilter = TanimotoSamplingFilter(filter_name="Tani_Sample",
                                                  sample_size=sample_size,
                                                  weight=weight)
        pk.filters.append(taniSampleFilter)

    if mcs_filter:
        mcsFilter = MCSFilter(filter_name="MCS", crit_mcs=crit_mcs)
        pk.filters.append(mcsFilter)

    if metabolomics_filter:
        metFilter = MetabolomicsFilter(filter_name="ADP1_Metabolomics_Data",
                                       met_data_name=met_data_name,
                                       met_data_path=met_data_path,
                                       possible_adducts=possible_adducts,
                                       mass_tolerance=mass_tolerance)
        pk.filters.append(metFilter)

    # Transform compounds (the main step)
    pk.transform_all(processes, generations)

    # Remove cofactor redundancies
    # Eliminates cofactors that are being counted as compounds
    pk.remove_cofactor_redundancy()

    if (tani_filter or mcs_filter or tani_sample):
        if prune_to_targets:
            pk.prune_network_to_targets()

    # Write results to database
    if write_db:
        pk.save_to_mine(processes=processes, indexing=indexing)
        client = pymongo.MongoClient(mongo_uri)
        db = client[database]
        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                 "Run Time": f"{round(time.time() - start, 2)}",
                                 "Generations": f"{generations}",
                                 "Operator file": f"{rule_list}",
                                 "Coreactant file": f"{coreactant_list}",
                                 "Input compound file": f"{input_cpds}"
                                 })

        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                 "Message": message})

        if (tani_filter or mcs_filter or tani_sample):
            db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                     "React Targets": react_targets,
                                     "Tanimoto Filter": tani_filter,
                                     "Tanimoto Values": f"{tani_threshold}",
                                     "MCS Filter": mcs_filter,
                                     "MCS Values": f"{crit_mcs}",
                                     "Sample By": tani_sample,
                                     "Sample Size": sample_size,
                                     "Sample Weight": weight_representation,
                                     "Pruned": prune_to_targets
                                     })

    if write_local:
        pk.assign_ids()
        pk.write_compound_output_file(output_dir + '/compounds.tsv')
        pk.write_reaction_output_file(output_dir + '/reactions.tsv')

    print('----------------------------------------')
    print(f'Overall run took {round(time.time() - start, 2)} seconds.')
    print('----------------------------------------')
    