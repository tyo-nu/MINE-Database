from minedatabase.pickaxe import Pickaxe

# Where are the input rxns coming from and coreactants
# Compounds that are going to be expanded
input_cpds = './data/starting_cpds_ten.csv'

# Cofactors and rules
coreactant_list = './minedatabase/data/EnzymaticCoreactants.tsv'
rule_list = './minedatabase/data/EnzymaticReactionRules.tsv'

# Outputs if you are writing the results locally
write_local = False
pickaxe_rxns = 'rxns_out.tsv'
pickaxe_cpds = 'cps_out.tsv'

# Database to write results to
write_db = True
database_overwrite = True
database = 'ten_cpds'

creds = open('credentials.csv').readline().split(',')
creds = [cred.strip('\n') for cred in creds]
# con_string is the login information for the mongodb. The default is localhost:27017
# Connecting remotely requires the location of the database as well as username/password
# if security is being used. Username/password are stored in credentials.csv
# in the following format: username,password
# Local MINE server
con_string = 'mongodb://localhost:27017'
# Connecting to the northwestern MINE server
# con_string = f'mongodb://{creds[0]}:{creds[1]}@minedatabase.ci.northwestern.edu:27017/?authSource=admin'

# Pickaxe Options
generations = 1
racemize = False
verbose = False
bnice = True
kekulize = True
neutralise = True
image_dir = None
quiet = True
max_workers = 12

# Tanimoto Filtering options
tani_filter = False
target_cpds = './data/target_list_many.csv'
crit_tani = 0.9

# Running pickaxe
# Initialize the Pickaxe class instance
pk = Pickaxe(coreactant_list=coreactant_list,
            rule_list=rule_list, racemize=racemize,
            errors=verbose, explicit_h=bnice,
            kekulize=bnice, neutralise=bnice,
            image_dir=image_dir, database=database,
            database_overwrite=database_overwrite,
            con_string=con_string, quiet=quiet)

# Load compounds
pk.load_compound_set(compound_file=input_cpds)
if tani_filter:
    pk.load_target_compounds(target_compound_file=target_cpds, crit_tani=crit_tani)

# Transform based on reaction rules
pk.transform_all(max_generations=generations,
                     num_workers=max_workers)

# Write results
if write_db:
    pk.save_to_mine(num_workers=max_workers)

if write_local:
    pk.write_compound_output_file(pickaxe_cpds)
    pk.write_reaction_output_file(pickaxe_rxns)