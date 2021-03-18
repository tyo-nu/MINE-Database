Running Pickaxe via Command Line
================================
Pickaxe supports running through a command line interface, but does not offer the
full functionality available through writing a python script :doc:`pickaxe_run.rst`.

Command Line Interface Features
-------------------------------

::

    $ python pickaxe.py -h
    usage: pickaxe.py [-h] [-C COREACTANT_LIST] [-r RULE_LIST] [-c COMPOUND_FILE] [-v] [-H] [-k] [-n] [-m PROCESSES] [-g GENERATIONS] [-q] [-s SMILES] [-p PRUNING_WHITELIST] [-o OUTPUT_DIR] [-d DATABASE] [-u MONGO_URI] [-i IMAGE_DIR]

    optional arguments:
    -h, --help            show this help message and exit
    -C COREACTANT_LIST, --coreactant_list COREACTANT_LIST
                            Specify a list of coreactants as a .tsv
    -r RULE_LIST, --rule_list RULE_LIST
                            Specify a list of reaction rules as a .tsv
    -c COMPOUND_FILE, --compound_file COMPOUND_FILE
                            Specify a list of starting compounds as .tsv or .csv
    -v, --verbose         Display RDKit errors & warnings
    -H, --explicit_h      Specify explicit hydrogen for use in reaction rules.
    -k, --kekulize        Specify whether to kekulize compounds.
    -n, --neutralise      Specify whether to neturalise compounds.
    -m PROCESSES, --processes PROCESSES
                            Set the max number of processes.
    -g GENERATIONS, --generations GENERATIONS
                            Set the numbers of time to apply the reaction rules to the compound set.
    -q, --quiet           Silence warnings about imbalanced reactions
    -s SMILES, --smiles SMILES
                            Specify a starting compound SMILES.
    -p PRUNING_WHITELIST, --pruning_whitelist PRUNING_WHITELIST
                            Specify a list of target compounds to prune reaction network down to.
    -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                            The directory in which to write files.
    -d DATABASE, --database DATABASE
                            The name of the database to store results.
    -u MONGO_URI, --mongo_uri MONGO_URI
                            The URI of the mongo database to connect to. Defaults to mongodb://localhost:27017
    -i IMAGE_DIR, --image_dir IMAGE_DIR
                            Specify a directory to store images of all created compounds

Examples
--------

Generate and Save Data to Local directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is the simplest example of using the command line interface. It accepts coreactant, rule, and compound files and  expands
to generations before saving the results in .tsv files in a provided directory.
::

    python pickaxe.py -r /path/to/rules.tsv -C path/to/coreactants.tsv -c /path/to/compounds.tsv -g 2 -o /path/to/output/

Generate and Save Data to a Mongo Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is possible to save to a mongo database, either locally or remotely. This option works with writing a .tsv as well, and will
write to both locations.

**Local Mongo Server**
Running the following will use, by default, mongodb://localhost:27017 as the mongo URI.
::

    python pickaxe.py -r /path/to/rules.tsv -C path/to/coreactants.tsv -c /path/to/compounds.tsv -g 2 -d database_name

**Specific Mongo Server**
Alternatively, a [specific Mongo URI can be specified](https://docs.mongodb.com/manual/reference/connection-string/), 
allowing for the use of password protected databases and remote databases.
::

    python pickaxe.py -r /path/to/rules.tsv -C path/to/coreactants.tsv -c /path/to/compounds.tsv -g 2 -d database_name -u mongodb://myDBReader:D1fficultP%40ssw0rd@mongodb0.example.com:27017/?authSource=admin

Generate with Multiple Processes and Pruning Final Network
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This example uses 4 processes to run and prunes the final network to contain only compounds
that are specified and any compounds required to generate them from the starting compounds.
::

    python pickaxe.py -r /path/to/rules.tsv -C path/to/coreactants.tsv -c /path/to/compounds.tsv -g 2 -o /path/to/output/ -m 4 -p /path/to/pruning_targets.tsv

