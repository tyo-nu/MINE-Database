from minedatabase.pickaxe import Pickaxe

# Where are the input rxns coming from and coreactants
# Compounds that are going to be expanded
input_cpds = './pickaxe_tests/in_cpds_mid.csv'

# Cofactors used in the rules
coreactant_list = './minedatabase/data/EnzymaticCoreactants.tsv'
# The reaction rules themselves
rule_list = './minedatabase/data/EnzymaticReactionRules.tsv'

# Outputs if you are writing the results locally
write_local = False
pickaxe_rxns = './pickaxe_tests/test_rxns.tsv'
pickaxe_cpds = './pickaxe_tests/test_cpds.tsv'

# Database to write results to
write_db = True
database = 'mid_1gen'

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
target_cpds = './data/target_list_trunc.csv'
crit_tani = 0.0
tani_filter = False

# Run pickaxe
pk = Pickaxe(coreactant_list=coreactant_list,
            rule_list=rule_list, racemize=racemize,
            errors=verbose, explicit_h=bnice,
            kekulize=bnice, neutralise=bnice,
            image_dir=image_dir, database=database,
            quiet=quiet)

# Load compounds
pk.load_compound_set(compound_file=input_cpds)
if tani_filter:
    pk.load_target_compounds(target_compound_file=target_cpds, crit_tani=crit_tani)

# Transform based on reaction rules
pk.transform_all(max_generations=generations,
                     num_workers=max_workers)

# Write results
if write_db:
    pk.save_to_mine(database)

if write_local:
    pk.write_compound_output_file(pickaxe_cpds)
    pk.write_reaction_output_file(pickaxe_rxns)