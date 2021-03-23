Running Pickaxe
================

Pickaxe is the program that is used to generate the data that is stored in the MINE-Database. The database is used for
metabolomics applications, but pickaxe can be extended for use in general reaction network generation and analysis. An 
example run, pickaxe_run.py, is found in the github. This python script provides a template for producing pickaxe runs, 
exposing the key parameters for a user to modify and inputs these into the pickaxe class to run, greatly simplifying the
process. 

pickaxe_run.py highlights the key components for a pickaxe run block-by-block. This document also serves to highlight
and explain the components of running pickaxe. Generally, pickaxe_run.py operates in the following steps:

1. Specify where the output of the run will be stored 
2. Specifying the various run inputs
3. Core Pickaxe options
4. Specification of Filters

This document gives the relevant code snippets from a template and expands on existing comments. Additionally, brief 
examples of relevant inputs will be created. For more detailed descriptions please see :doc:`inputs` and :doc:`api_ref/filters`.

.. tip::
    To create custom filters, see :doc:`custom_filters`.

Example Template
----------------
This document details the specifics of a template file, pickaxe_run.py, that highlights
common Pickaxe runs.

:download:`pickaxe_run.py can be downloaded here. <_static/files/pickaxe_run.py>`

Run Output
----------
There are two ways to output data:

1. Writing to a mongo database that is specified by a `mongo uri`_, either local or in mongo_uri.csv
2. Local .tsv files

.. code-block:: python

    # Whether or not to write to a mongodb
    write_db = False
    database_overwrite = False
    # database = "APAH_100Sam_50rule"
    database = "example_pathway"
    # Message to insert into metadata
    message = ("Example run to show how pickaxe is run.")

    # mongo DB information
    use_local = False
    if write_db == False:
        mongo_uri = None
    elif use_local:
        mongo_uri = 'mongodb://localhost:27017'
    else:
        mongo_uri = open('mongo_uri1.csv').readline().strip('\n')

    # Write output .csv files locally
    write_to_csv = False
    output_dir = '.'

.. _mongo uri: https://docs.mongodb.com/manual/reference/connection-string/

Run Input
---------
There are three key inputs for a Pickaxe run to be specified:

1. **input_cpds** specifying the compounds to be reacted
2. **coreactant_list** are coreactants that are required for the reaction rules
3. **rule_list** that specifies the reaction rules to be applied

Input Compounds Example
^^^^^^^^^^^^^^^^^^^^^^^
The file specified for :code:`input_cpds` must be a .tsv or a .csv format. 
The file consists of an id and a SMILES string. An example of a .csv file is

::

    id,SMILES
    0,CC(=O)OC
    1,CCO

Coreactant and Rule lists
^^^^^^^^^^^^^^^^^^^^^^^^^
Pickaxe is provided with a default rule list generated from approximately 70,000 MetaCyc reactions.

The following code allows you to select then number of rules by either a number or by coverage:

.. code-block:: python

    from minedatabase.rules import metacyc_generalized
    # Select by number
    rule_list, coreactant_list, rule_name = metacyc_generalized(n_rules=20)

    # Select by fraction coverage
    rule_list, coreactant_list, rule_name = metacyc_generalized(fraction_coverage=0.5)



When choosing how many reactions to use, you can refer to the following table:

+-----------------+---------------------+
| Number of Rules | Percent Coverage of |
|                 | MetaCyc Reactions   |
+-----------------+---------------------+
| 20              | 50                  |
+-----------------+---------------------+
| 84              | 75                  |
+-----------------+---------------------+
| 100             | 78                  |
+-----------------+---------------------+
| 272             | 90                  |
+-----------------+---------------------+
| 500             | 95                  |
+-----------------+---------------------+
| 956             | 99                  |
+-----------------+---------------------+
| 1221            | 100                 |
+-----------------+---------------------+

.. note:: 
    Rules and coreactants can be generated manually as well, which is outlined in 
    :doc:`inputs`.

Code snippet from Pickaxe_run.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These input files are specified as follows:

.. code-block:: python

    input_cpds = './example_data/starting_cpds_single.csv'

    # Generate rules automatically from metacyc generalized. n_rules takes precedence over 
    # fraction_coverage if both specified. Passing nothing returns all rules.
    rule_list, coreactant_list, rule_name = metacyc_generalized(
        n_rules=20,
        fraction_coverage=None
    )

If you generated a file manually then specify the file directly as follows:

.. code-block:: python

    rule_list = "path/to/rules"
    coreactant_list = "path/to/coreactants"
    rule_name = "rule name"


Core Pickaxe Options
--------------------
Of these options the majority of uses will only require the changing of the following:

1. **generations** is the number of generations to expand, e.g. 2 generations will apply reaction rules twice
2. **num_works** specifies the number of processors to use

However, the remaining can be changed if needed:

3. **verbose** specifies if RDKit is suppressed or not
4. **kekulize** specifies whether or not to kekulize RDKit molecules
5. **neutralise** specifies whether or not to neutralise molecules
6. **image_dir** specifies the directory where to draw images of generated compounds
7. **quiet** specifies whether or not to suppress output
8. **indexing** specifies whether or not to index the databases 

.. code-block:: python

    generations = 1
    processes = 4     # Number of processes for parallelization
    verbose = False     # Display RDKit warnings and errors
    explicit_h = False
    kekulize = True
    neutralise = True
    image_dir = None
    quiet = True
    indexing = False

Built-In Filters
----------------
Three general filters are supplied with Pickaxe:

1. A tanimoto threshold filters
2. A tanimoto sampling filters
3. A metabolomics filters

Specified filters are applied before each generation (and at the end of the run if specified) to reduce the number of compounds
to be expanded. This allows for the removal of compounds that aren't of interest to reduce the number of non-useful compounds in the resultant network.
Additionally, custom filters can be written. To write your own filter see: 

General Filter Options
^^^^^^^^^^^^^^^^^^^^^^
These options apply to every filter and are independent of the actual filter itself.

1. **target_cpds** specifies where the target compound list is. This file is a csv
    with the header id,SMILES
2. **react_targets** specifies whether a compound generated in the expansion should be further reacted
3. **prune_to_targets** specifies whether the network should be reduced to a minimal network containing only compounds directly connected to the targets from a source
4. **filter_after_final_gen** whether to apply the filter to the final application of reaction rules

.. code-block:: python

    # Path to target cpds file (not required for metabolomics filter)
    target_cpds = './example_data/target_list_single.csv'

    # Should targets be flagged for reaction
    react_targets = True

    # Prune results to remove compounds not required to produce targets
    prune_to_targets = True

    # Filter final generation?
    filter_after_final_gen = True


Tanimoto Threshold Filter
^^^^^^^^^^^^^^^^^^^^^^^^^
The rational behind this filter is to generate a list of Tanimoto similarity scores (ranging from 0 to 1) for each generation
in comparison to the targets and use this to trim compounds to only those above a certain similarity threshold. 
The maximum similarity of a given compound compared to all the targets is used. Similarity is calculated
by using the default RDKFingerprints. 

Before each generation the maximum similarity for each compound set to be reacted is compared to a threshold. Compounds greater than or equal
to the threshold are reacted. 

1. **tani_filter** whether or not to use this filter
2. **tani_threshold** is the threshold to cut off. Can be a single value or a list. If a list then the filter will use the next value in this list for each new generation
3. **increasing_tani** specifies whether the tanimoto value of compounds must increase each generation. I.e. a child compound must be more similar to a target than at least one of its parents

.. code-block:: python

    # Apply this filter?
    tani_filter = False

    # Tanimito filter threshold. Can be single number or a list with length at least
    # equal to the number of generations (+1 if filtering after expansion)
    tani_threshold = [0, 0.2, 0.7]

    # Make sure tani increases each generation?
    increasing_tani = False

Tanimoto Sampling Filter
^^^^^^^^^^^^^^^^^^^^^^^^
For large expansions the tanimoto threshold filter does not work well. For example, expanding 10,000 compounds from KEGG with 272 rules from metacyc yields 5 million compounds. To expand this another generation
the number of compounds has to be heavily reduced for the system resources to handle it and for analysis to be reasonable. 
The threshold filter will have to be at a large value, e.g. greater than 0.9, which leads to reduced chemical diversity in the final network.

To avoid this problem, the Tanimoto Sampling Filter was implemented. The same approach as the threshold filter is taken to get a list of maximum similarity score for compounds and the list of targets.
This tanimoto score is scaled and then the distribution is sampled by inverse complementary distribution function sampling to select N compounds. This approach affords more diversity than the threshold
and can be tuned by scaling the tanimoto similarity score scaling function. By default the function is :math:`T^{4}`. 

The filter is specified as follows:

1. **tani_sample** specifies whether to use the filter
2. **sample_size** specifies the number of compounds to expand each generation. If sample_size is greater than the total number of compounds all compounds are reacted
3. **weight** specifies the weighting function for the sampling. This function accepts a float and returns a float
4. **weight_representation** specifies how to display the weighting function in the database or stdout 

.. code-block:: python

    # Apply this sampler?
    tani_sample = False

    # Number of compounds per generation to sample
    sample_size = 5

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

Metabolomics Filter
^^^^^^^^^^^^^^^^^^^
If you have a metabolomics dataset you would like to filter compounds against, you can use this filter.
It will force pickaxe to only keep compounds with masses (and, optionally, retention time (RT)) within a set
tolerance of a list of peaks. For example, if you had a dataset containing 3 peaks at 100, 200, and 300 m/z,
you could do an expansion and only keep compounds with masses within 0.001 Da of those 3 values.

This is useful for trying to annotate unknown peaks starting from a set of known compounds in a specific organism
from which metabolomics data was collected.

The filter is specified as follows. The following arguments are required:

1. **metabolomics_filter** specifies whether to use this filter

2. **met_data_path** specifies where to find your list of peaks in CSV format.

Format of CSV:

    Peak ID, Retention Time, Aggregate M/Z, Polarity, Compound Name, Predicted Structure (smile), ID
    
    Peak1, 6.33, 74.0373, negative, propionic acid, CCC(=O)O, yes
    
    Peak2, 26.31, 84.06869909, positive, , , no
    
    ...

Note that only unidentified peaks will be used by the filter.

3. **possible_adducts** specifies the possible adducts to consider when matching peaks, as different adducts cause different mass changes. For a list of options, see the first columns of  "Negative Adducts full.txt" and "Positive Adducts full.txt" in minedatabase/data/adducts.

4. **mass_tolerance** specifies (in Da) the mass tolerance to use for matching peaks. For example, if 0.001, only compounds with masses between 99.999 and 100.001 would match a peak at 100 m/z.

The following optional arguments allow you to add retention time as an extra constraint in the filter.
Note that this requires that you have built a RandomForestRegressor machine learning model to predict
retention time for arbitrary compounds, using mordred fingerprints as input.

5. **rt_predictor_pickle_path** specifies the path to the built model (pickled). Make sure this is None, if you don't want to match based on retention time.

6. **rt_threshold** specifies the retention time tolerance (in whatever units RT is in the file at met_data_path)

7. **rt_important_features** specifies which mordred descriptors to use as input into the model (must be in same order as model expects them to be). If None, will use all (including 3D) mordred descriptors.

.. code-block:: python

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