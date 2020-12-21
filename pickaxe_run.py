from minedatabase.pickaxe import Pickaxe
import pymongo
import datetime
import time

###############################################################################
##### pickaxe_run.py
# This python script provides a skeleton to build Pickaxe runs around
# The general format of a script will be:
#   1. Connect to mongoDB (if desired)
#   2. Load reaction rules and cofactors
#   3. Load starting compounds
#   4. Load Tanimoto filtering options
#   5. Transform Compounds
#   6. Write results
###############################################################################

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
database = "test_true"
# Message to insert into metadata
message = ("Debugging filter. Should yield no orphans. 100 Sampling. Terminal,"
           " Tani filtering.")

# mongo DB information
use_local = True
if use_local:
    mongo_uri = 'mongodb://localhost:27017'
else:
    # load file of form user,pass
    creds = open('credentials.csv').readline().split(',')
    creds = [cred.strip('\n') for cred in creds]
    # URI of remote mongo instance
    mongo_uri = f"mongodb://{creds[0]}:{creds[1]}@minedatabase.ci.northwestern.edu:27017/?authSource=admin"

# Write output .csv files locally
write_local = False
output_dir = '.'
###############################################################################

###############################################################################
#    Cofactors, rules and inputs
# Original rules derived from BNICE
# coreactant_list = './minedatabase/data/EnzymaticCoreactants.tsv'
# rule_list = './minedatabase/data/EnzymaticReactionRules.tsv'

# Rules from Joseph Ni
coreactant_list = './minedatabase/data/MetaCyc_Coreactants.tsv'
# rule_list = './minedatabase/data/intermediate_rules_uniprot.tsv'
rule_list = r'/Users/kevbot/Box Sync/Research/Projects/MINE/MINE-Database/minedatabase/data/metacyc_272rules_90percentMapping.tsv'

# Input compounds
input_cpds = './local_data/APAH.csv'

# Partial operators
# Partial operators allow use of multiple compounds in an any;any expansion
partial_rules = False
mapped_rxns = 'minedatabase/data/metacyc_mapped.tsv'
###############################################################################

###############################################################################
##### Core Pickaxe Run Options
generations = 2
num_workers = 12     # Number of processes for parallelization
verbose = False      # Display RDKit warnings and errors
explicit_h = False
kekulize = True
neutralise = True
image_dir = None
quiet = True
indexing = False
###############################################################################

###############################################################################
##### Filtering Options
target_cpds = './local_data/ADP1_cpds_out_final.csv'
# Should targets be flagged for reaction?
react_targets = True

# Specify if network expansion is done in a retrosynthetic direction
retrosynthesis = False
# Prune results to only give expanded compounds/rxns
prune_by_filter = True

##### Tanimoto Filtering options
tani_filter = False

# Tanimito filter threshold. Can be single number of a list
# of length generations.
crit_tani = [0, 0.5, 0.7]
# crit_tani = [0, 0.5] # expands first with no filter then a 0.5 filter

tani_sample = True
sample_size = 100
# Give a function that accepts a single argument and returns a single result
# Inputs are [0, 1]
# weight = None will use a f(x) = x^4 to weight.
weight = None
weight_for_db = "T^4"

##### MCS Filter options
mcs_filter = False
# MCS. Can be single number of a list
# of length generations.

# Finds the MCS of the target and compound and identifies fraction of target
# the MCS composes
crit_mcs = [0.3, 0.8, 0.95]

# TODO: MCS options
###############################################################################

###############################################################################
##### Running pickaxe
# Initialize the Pickaxe class
if write_db is False:
    database = None

pk = Pickaxe(coreactant_list=coreactant_list,
            rule_list=rule_list,
            errors=verbose, explicit_h=explicit_h,
            kekulize=kekulize, neutralise=neutralise,
            image_dir=image_dir, database=database,
            database_overwrite=database_overwrite,
            mongo_uri=mongo_uri, quiet=quiet)

# Load compounds
pk.load_compound_set(compound_file=input_cpds)

# Load partial operators
if partial_rules:
    pk.load_partial_operators(mapped_rxns)

# Initialize filter
if tani_filter or mcs_filter or tani_sample:
    pk.load_target_and_filters(
        target_compound_file=target_cpds,
        tani_filter=tani_filter, crit_tani=crit_tani,
        tani_sample=tani_sample, sample_size=sample_size, weight=weight,
        mcs_filter=mcs_filter, crit_mcs=crit_mcs,
        retrosynthesis=False, react_targets=react_targets
    )

# Transform compounds
pk.transform_all(num_workers, generations)

# Remove cofactor redundancies
# Eliminates cofactors that are being counted as compounds
pk.remove_cofactor_redundancy()

if (tani_filter or mcs_filter or tani_sample):
    if prune_by_filter:
        pk.prune_network_to_targets()

# Write results to database
if write_db:
    pk.save_to_mine(num_workers=num_workers, indexing=indexing)
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
                                 "Tanimoto Values": f"{crit_tani}",
                                 "MCS Filter": mcs_filter,
                                 "MCS Values": f"{crit_mcs}",
                                 "Sample By": tani_sample,
                                 "Sample Size": sample_size,
                                 "Sample Weight": weight_for_db,
                                 "Pruned": prune_by_filter
                                 })

if write_local:
    pk.assign_ids()
    pk.write_compound_output_file(output_dir + '/compounds.tsv')
    pk.write_reaction_output_file(output_dir + '/reactions.tsv')

print('----------------------------------------')
print(f'Overall run took {round(time.time() - start, 2)} seconds.')
print('----------------------------------------')