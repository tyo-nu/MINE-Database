from minedatabase.pickaxe import Pickaxe
import pymongo
import datetime
import time

start = time.time()
# Where are the input rxns coming from and coreactants
# Compounds that are going to be expanded
input_cpds = './example_data/starting_cpds_ten.csv'

# Cofactors and rules
# coreactant_list = './minedatabase/data/EnzymaticCoreactants.tsv'
# rule_list = './minedatabase/data/EnzymaticReactionRules.tsv'
coreactant_list = './minedatabase/data/MetaCyc_Coreactants.tsv'
rule_list = './minedatabase/data/metacyc_generalized_rules_500.tsv'

# Database to write results to
write_db = True
database_overwrite = True
database = 'KMS_devJ_single'

# creds = open('credentials.csv').readline().split(',')
# creds = [cred.strip('\n') for cred in creds]
# mongo_uri is the login information for the mongodb. The default is localhost:27017
# Connecting remotely requires the location of the database as well as username/password
# if security is being used. Username/password are stored in credentials.csv
# in the following format: username,password
# Local MINE server
mongo_uri = 'mongodb://localhost:27017'
# Connecting to the northwestern MINE server
# mongo_uri = f"mongodb://{creds[0]}:{creds[1]}@minedatabase.ci.northwestern.edu:27017/?authSource=admin"

# Pickaxe Options
generations = 1
verbose = False
explicit_h = False
kekulize = True
neutralise = True
image_dir = None
quiet = True
indexing = False
num_workers = 1

# Tanimoto Filtering options
tani_filter = False
# Prune results to only give expanded compounds/rxns
# Currently also includes all of the last generation
tani_prune = False
target_cpds = './example_data/target_list_many.csv'
# crit_tani is either a single number 
# OR a list that is the length of the number of generations
crit_tani = [0.5, 0.9]

# Running pickaxe
# Initialize the Pickaxe class instance
pk = Pickaxe(coreactant_list=coreactant_list,
            rule_list=rule_list,
            errors=verbose, explicit_h=explicit_h,
            kekulize=kekulize, neutralise=neutralise,
            image_dir=image_dir, database=database,
            database_overwrite=database_overwrite,
            mongo_uri=mongo_uri, quiet=quiet)

# Load compounds
pk.load_compound_set(compound_file=input_cpds)
pk.transform_all(num_workers, generations)
# if tani_filter:
#     pk.load_target_compounds(target_compound_file=target_cpds, crit_tani=crit_tani)

# # Transform based on reaction rules
# pk.transform_all(max_generations=generations,
#                      num_workers=max_workers)

# # Write results
if write_db:
    if tani_filter and tani_prune:
        pk.prune_network_to_targets()
    pk.save_to_mine(num_workers=num_workers, indexing=indexing)
    client = pymongo.MongoClient(mongo_uri)
    db = client[database]
    db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                    "Generations": f"{generations}",
                                    "Operator file": f"{rule_list}",
                                    "Coreactant file": f"{coreactant_list}",
                                    "Input compound file": f"{input_cpds}"}
                                    )
    db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                            "Message": ("Expansion for bioprivileged molecules."
                                        "Targeting 1k molecules identified by XZ using original 250 rules.")})


# print(f'Overall run took {round(time.time() - start, 2)} seconds.')