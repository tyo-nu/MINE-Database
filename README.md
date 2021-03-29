# MINE Databases
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Documentation](https://readthedocs.org/projects/mine-database/badge)

The MINE Database contains code for generating (through Pickaxe) and storing and retrieving compounds from a database.
Pickaxe applies reaction operators, representing reaction transformation patterns, to a list of user-specified compounds
in order to predict reactions.  

## Documentation
For general information on MINE Databases, please consult [JJeffryes et al. 2015.](http://jcheminf.springeropen.com/articles/10.1186/s13321-015-0087-1)
Documentation, hosted on read the docs at https://mine-database.readthedocs.io/en/develop/,
gives more detailed descriptions and examples uses of the software.

## Installation
MINE-Database requires the use of [rdkit](https://rdkit.org/), which currently is unavailable to install on pip. It is recommended to use the
conda environment installation.

## Conda Installation
MINE-Database can be installed via conda:
```
conda install -c condaforge mine-database
```

Alternatively MINE-Database can be installed from source:
```
conda env create --name minedatabase --file requirements.yml
```

## Pip and Conda Installation from Source
The majority of requirements can be installed through pip:
```
pip install -r requirements.txt
```
RDKit must still be installed using Anaconda
```
conda env create -c rdkit -n minedatabase rdkit
```


<!-- ## Repository Structure
This repository primarily consists of the minedatabases python module 
and its 8 submodules:
- compound_io: Contains functions to load and export chemical structures 
from MINE databases. Has command-line interface.
- databases: Contains functions which impose a schema on the underlying 
Mongo databases which store MINE data.
-filters.py: Contains filters to apply while generating reaction network
- pickaxe: Allows for the application of reaction rules to compound sets 
and the annotation of the resulting compounds an reactions. Has command-line interface.
- queries: Contains logic for text and chemical structure queries of the 
MINE database.
- reactions: Contains methods to apply reaction rules to compounds
- utils: Various utility functions such as hashing & type conversions -->

<!-- ### Compound_io command-line usage
compound_io may be called independently to import or export chemical 
strictures in MINE database format. These may be helpful for sharing 
predictions or maintaining current external databases for cross-referencing.
The call format is `python compound_io.py import-<format> <input path> 
<database>` for imports and `python compound_io.py export-<format> 
<database> <outfile path> <optionally: maximum compounds per file>` for exports.
Valid formats are:
- smi: SMILES line-code
- mol: MDL molecule files (outputs individual files in specified directory)
- sdf: Structured Data File (concatenated mol files)
- tsv: FOR EXPORT ONLY, a tab separated file compatible with ModelSEED -->

## Running Pickaxe
### Running Pickaxe through a python file (recommended)
An example file, pickaxe_run.py, provides a framework for running pickaxe through a python file.
The starting compounds, rules and cofactors, optional database information, and Pickaxe run options are specified.
After running the results are stored in a specified database or written to .tsv files. This file is
explained in more detail in the [documentation](https://mine-database.readthedocs.io/en/develop/pickaxe_run.html).

### Pickaxe command-line usage (not recommended - see above section)
Pickaxe.py can be called independently to generate predictions with or 
without database storage. To list all options call `python -m minedatabase.pickaxe -h`. Note that
due to relative imports, it needs to be run as a module (-m flag) from the MINE-Database directory.
To predict metacyc reactions for one generation on compounds in the iML1515
model one would call 
```
python pickaxe.py -C ./data/metacyc_generalized_rules.tsv -r ./data/metacyc_coreactants.tsv -g 1 -c ../example_data/iML1515_ecoli_GEM.csv
```

### Testing
`pytest` to run all tests. Ensure that pytest is installed.
To add coverage, run:
```
pytest --cov-report term --cov-report xml:tests/cov.xml --cov=minedatabase minedatabase/tests/
```
Ensure that coverage and pytest-cov are both installed.
