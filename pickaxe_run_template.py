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
from rdkit import DataStructs
from rdkit.Chem import AllChem

from minedatabase.filters import (
    AtomicCompositionFilter,
    MCSFilter,
    MetabolomicsFilter,
    MWFilter,
    SimilarityFilter,
    SimilaritySamplingFilter
)
# Uncomment to use these. Pickaxe doesn't come packaged with dependencies by default.
# from minedatabase.filters import ThermoFilter
# from minedatabase.filters import ReactionFeasibilityFilter

from minedatabase.pickaxe import Pickaxe
from minedatabase.rules import metacyc_generalized, metacyc_intermediate

start = time.time()

###############################################################################
#### Database and output information
# The default mongo is localhost:27017
# Connecting remotely requires the location of the database
# as well as username/password if security is being used.
#
# The mongo_uri is read from mongo_uri.csv


# Database writing options
write_db = False
# Database name and message to print in metadata
database = "example_db"
message = "Example run to show how pickaxe is run."
# Force overwrite existing database
database_overwrite = False
# Use local DB, i.e. localhost:27017
use_local = False

# Writing compound and reaction csv files locally
write_to_csv = False
output_dir = "."
###############################################################################

###############################################################################
#### Starting Compounds, Cofactors, and Rules
# Input compounds in a csv folder with headings:
# id,smiles
input_cpds = "./example_data/starting_cpds_single.csv"

# Rule specification and generation. Rules can be manually created or
# metacyc_intermediate or metacyc_generalized can provide correctly formatted
# biological reactions derived from metacyc.
#
# See the documentation for description of options.
rule_list, coreactant_list, rule_name = metacyc_intermediate(
    n_rules=None,
    fraction_coverage=0.2,
    anaerobic=True,
    exclude_containing = ["aromatic", "halogen"]
)

###############################################################################

###############################################################################
# Core Pickaxe Run Options
generations = 1              # Total rounds of rule applications
processes = 1                # Number of processes for parallelization
verbose = False              # Display RDKit warnings and errors

# These are for MINE-Database generation and advanced options.
# Be careful changing these.
inchikey_blocks_for_cid = 1  # Number of inchi key blocks to gen cid
explicit_h = False           # use explicit hydrogens in rules
kekulize = True              # kekulize molecules
neutralise = True            # Neutralise all molecules when loading
quiet = True                 # Silence errors
indexing = False             #
###############################################################################

###############################################################################
#### Filtering and Sampling Options

#############################################
# Global Filtering Options

# Path to target cpds file (not required for all filters)
target_cpds = "./example_data/target_list_many.csv"

# Load compounds even without a filter
# Can be paired with prune_to_targets to reduce end network
load_targets_without_filter = True

# Should targets be flagged for reaction
react_targets = False

# Prune generated network to contain only compounds
# that are used to generate a target
prune_to_targets = False

# Apply filter after final generation, before pruning
filter_after_final_gen = True

##########################################
# Molecular Weight Filter options.
# Removes compounds not in range of [min_MW, max_MW]

# Apply this filter?
MW_filter = True

# Minimum MW in g/mol. None gives no lower bound.
min_MW = 100

# Maximum MW in g/mol. None gives no upper bound.
max_MW = 150

##########################################
# Atomic Composition Filter options.
# Filters compounds that do not fall within atomic composition range
# Only elements specified will be filtered

# Apply this filter?
atomic_composition_filter = False

# Atomic composition constraint specification
atomic_composition_constraints = {
    "C": [4, 7],
    "O": [5, 5]
}

##########################################
# Thermodynamics Filter options.
# Uses eQuilibrator to filter by ∆Gr

# Information for use of eQuilibrator
# URI of eQuilibrator DB, either postgres URI, or sqlite file name location
eq_uri = "compounds.sqlite"

# Maximum allowable ∆Gr in kJ/mol
dg_max = 10

# conditions
p_h = 7
p_mg = 3
ionic_strength = 0.15

# comment below line and uncomment other definition if using thermo filter
thermo_filter = None
# thermo_filter = ThermoFilter(
#             eq_uri=eq_uri,
#             dg_max=dg_max,
#             p_h=p_h,
#             p_mg=p_mg,
#             ionic_strength=ionic_strength
#         )

##########################################
# Feasibility Filter options.
# Checks Feasibility of reaction

# Apply this filter?
feasibility_filter = False

# Which generations to filter, empty list filtters all
generation_list = []
last_generation_only = True

# comment below line and uncomment other definition if using thermo filter
feasibility_filter = None
feasibility_filter = ReactionFeasibilityFilter(
    generation_list=generation_list,
    last_generation_only=last_generation_only
)


##########################################
# Similarity Filtering options.
# Filters by similarity score, uses default RDKit fingerprints and tanimoto by default

# Apply this filter?
similarity_filter = True

# Methods to calculate similarity by, default is RDkit and Tanimoto
# Supports Morgan Fingerprints and Dice similarity as well.
cutoff_fingerprint_method = "Morgan"
# arguments to pass to fingerprint_method
cutoff_fingerprint_args = {"radius": 2}
cutoff_similarity_method = "Tanimoto"

# Similarity filter threshold. Can be single number or a list with length at least
# equal to the number of generations (+1 if filtering after expansion)
similarity_threshold = [0, 0.2, 0.7]

# Only accepts compounds whose similarity is increased in comparison to their parent
increasing_similarity = False

##########################################
# Similarity Sampling Options
# Samples by similarity score
# Uses default RDKit fingerprints and tanimoto by default, but supports
# Morgan and dice

# Apply this sampler?
similarity_sample = True
# Number of compounds per generation to sample
sample_size = 100

# Default is RDKit
sample_fingerprint_method = "Morgan"
# arguments to pass to fingerprint_method
sample_fingerprint_args = {"radius": 2}
sample_similarity_method = "Tanimoto"

def weight(score):
    """weight is a function that accepts a similarity score as the sole argument
    and returns a scaled value. 
    """
    return score**4

# How to represent the function in text for database entry
weight_representation = "score^4"

##########################################
# Maximum common substructure (MCS) filter

# Apply this filter?
mcs_filter = False

# Finds the MCS of the target and compound and identifies fraction of target
# the MCS composes
crit_mcs = [0.3, 0.8, 0.95]

##########################################
# Metabolomics Filter Options

# Apply this filter?
metabolomics_filter = False

# Path to csv with list of detected masses (and optionally, retention times).
# For example: Peak ID, Retention Time, Aggregate M/Z, Polarity, Compound Name,
# Predicted Structure (smile), ID
#
# Peak1, 6.33, 74.0373, negative, propionic acid, CCC(=O)O, yes
# Peak2, 26.31, 84.06869909, positive, , , no
# ...
met_data_path = "./local_data/ADP1_Metabolomics_PeakList_final.csv"

# Name of dataset
met_data_name = "ADP1_metabolomics"

# Adducts to add to each mass in mass list to create final list of possible
# masses.
# See "./minedatabase/data/adducts/All adducts.txt" for options.
possible_adducts = ["[M+H]+", "[M-H]-"]

# Tolerance in Da
mass_tolerance = 0.001

# Retention Time Filter Options (optional but included in metabolomics filter)

# Path to pickled machine learning predictor (SMILES => RT)
rt_predictor_pickle_path = "../RT_Prediction/final_RT_model.pickle"

# Allowable deviation in predicted RT (units just have to be consistent with dataset)
rt_threshold = 4.5

# Mordred descriptors to use as input to model (must be in same order as in trained model)
# If None, will try to use all (including 3D) mordred descriptors
rt_important_features = ["nAcid", "ETA_dEpsilon_D", "NsNH2", "MDEO-11"]

###############################################################################

###############################################################################
# Verbose output
print_parameters = True


def print_run_parameters():
    """Write relevant parameters."""
    def print_parameter_list(plist):
        for i in plist:
            print(f"--{i}: {eval(i)}")

    print("\n-------------Run Parameters-------------")

    print("\nRun Info")
    print_parameter_list(["coreactant_list", "rule_name", "input_cpds"])

    print("\nExpansion Options")
    print_parameter_list(["generations", "processes"])

    print("\nGeneral Filter Options")
    print_parameter_list(
        [
            "filter_after_final_gen",
            "react_targets",
            "prune_to_targets",
        ]
    )

    if similarity_sample:
        print("\nTanimoto Sampling Filter Options")
        print_parameter_list(
            [
                "sample_size",
                "weight_representation",
                "sample_fingerprint_args",
                "sample_fingerprint_method",
                "sample_similarity_method"
            ]
        )

    if similarity_filter:
        print("\nTanimoto Threshold Filter Options")
        print_parameter_list(
            [
                "similarity_threshold",
                "increasing_similarity",
                "cutoff_fingerprint_args",
                "cutoff_fingerprint_method",
                "cutoff_similarity_method"
            ]
        )

    if mcs_filter:
        print("\nMaximum Common Substructure Filter Options")
        print_parameter_list(["crit_mcs"])

    if metabolomics_filter:
        print("\nMetabolomics Filter Options")
        print_parameter_list(["met_data_path", "met_data_name",
                              "possible_adducts", "mass_tolerance"])

    if MW_filter:
        print("\nMolecular Weight Filter Options")
        print_parameter_list(["min_MW", "max_MW"])

    if atomic_composition_filter:
        print("\nAtomic Composition Filter")
        print_parameter_list(["atomic_composition_constraints"])

    print("\nPickaxe Options")
    print_parameter_list(
        [
            "verbose",
            "explicit_h",
            "kekulize",
            "neutralise",
            "quiet",
            "indexing"
        ]
    )
    print("----------------------------------------\n")
###############################################################################


###############################################################################
#   Running pickaxe, don"t touch unless you know what you are doing
if __name__ == "__main__":
    # Use "spawn" for multiprocessing
    multiprocessing.set_start_method("spawn")

    # Define mongo_uri
    # mongo_uri definition, don't modify
    if write_db == False:
        mongo_uri = None
    elif use_local:
        mongo_uri = "mongodb://localhost:27017"
    else:
        mongo_uri = open("mongo_uri.csv").readline().strip("\n")

    # Change database to none if not writing
    if write_db is False:
        database = None

    ### Initialize the Pickaxe class
    # print parameters
    if print_parameters:
        print_run_parameters()

    pk = Pickaxe(
        coreactant_list=coreactant_list,
        rule_list=rule_list,
        errors=verbose,
        explicit_h=explicit_h,
        kekulize=kekulize,
        neutralise=neutralise,
        image_dir=None,
        inchikey_blocks_for_cid=inchikey_blocks_for_cid,
        database=database,
        database_overwrite=database_overwrite,
        mongo_uri=mongo_uri,
        quiet=quiet,
        react_targets=react_targets,
        filter_after_final_gen=filter_after_final_gen
    )

    # Load compounds
    pk.load_compound_set(compound_file=input_cpds)

    # Load target compounds for filters
    if (
        similarity_filter or mcs_filter or similarity_sample
        or load_targets_without_filter or MW_filter
        or atomic_composition_filter or thermo_filter
        or feasibility_filter
    ):
        pk.load_targets(target_cpds)

    # Apply filters in this order
    if MW_filter:
        pk.filters.append(MWFilter(min_MW, max_MW))

    if atomic_composition_filter:
        pk.filters.append(AtomicCompositionFilter(atomic_composition_constraints))

    if similarity_filter:
        taniFilter = SimilarityFilter(
            crit_similarity=similarity_threshold,
            increasing_similarity=increasing_similarity,
            fingerprint_method=cutoff_fingerprint_args,
            fingerprint_args=cutoff_fingerprint_args,
            similarity_method=cutoff_similarity_method
        )
        pk.filters.append(taniFilter)

    if similarity_sample:
        taniSampleFilter = SimilaritySamplingFilter(
            sample_size=sample_size,
            weight=weight,
            fingerprint_method=sample_fingerprint_method,
            fingerprint_args=sample_fingerprint_args,
            similarity_method=sample_similarity_method)
        pk.filters.append(taniSampleFilter)

    if mcs_filter:
        mcsFilter = MCSFilter(crit_mcs=crit_mcs)
        pk.filters.append(mcsFilter)

    if metabolomics_filter:
        if rt_predictor_pickle_path:
            with open(rt_predictor_pickle_path, "rb") as infile:
                rt_predictor = pickle.load(infile)
        else:
            rt_predictor = None

        metFilter = MetabolomicsFilter(
            filter_name="ADP1_Metabolomics_Data",
            met_data_name=met_data_name,
            met_data_path=met_data_path,
            possible_adducts=possible_adducts,
            mass_tolerance=mass_tolerance,
            rt_predictor=rt_predictor,
            rt_threshold=rt_threshold,
            rt_important_features=rt_important_features
        )
        pk.filters.append(metFilter)

    if feasibility_filter:
        pk.filers.append(feasibility_filter)

    if thermo_filter:
        pk.filters.append(thermo_filter)

    # Transform compounds (the main step)
    pk.transform_all(processes, generations)

    if pk.targets and prune_to_targets:
        pk.prune_network_to_targets()

    # Write results to database
    if write_db:
        pk.save_to_mine(processes=processes, indexing=indexing)
        client = pymongo.MongoClient(mongo_uri)
        db = client[database]
        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                 "Run Time": f"{round(time.time() - start, 2)}",
                                 "Generations": f"{generations}",
                                 "Rule Name": f"{rule_name}",
                                 "Input compound file": f"{input_cpds}"
                                 })

        db.meta_data.insert_one({"Timestamp": datetime.datetime.now(),
                                 "Message": message})

        if (similarity_filter or mcs_filter or similarity_sample):
            db.meta_data.insert_one(
                {
                    "Timestamp": datetime.datetime.now(),
                    "React Targets": react_targets,
                    "Tanimoto Filter": similarity_filter,
                    "Tanimoto Values": f"{similarity_threshold}",
                    "MCS Filter": mcs_filter,
                    "MCS Values": f"{crit_mcs}",
                    "Sample By": similarity_sample,
                    "Sample Size": sample_size,
                    "Sample Weight": weight_representation,
                    "Pruned": prune_to_targets
                }
            )

    if write_to_csv:
        pk.assign_ids()
        pk.write_compound_output_file(output_dir + "/compounds.tsv")
        pk.write_reaction_output_file(output_dir + "/reactions.tsv")

    print("----------------------------------------")
    print(f"Overall run took {round(time.time() - start, 2)} seconds.")
    print("----------------------------------------")
