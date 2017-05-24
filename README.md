# MINE Databases
[![Build Status](https://travis-ci.org/JamesJeffryes/MINE-Database.svg?branch=master)](https://travis-ci.org/JamesJeffryes/MINE-Database)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/fe968c8b593d48ebbb55aaf211c66f73)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=JamesJeffryes/MINE-Database&amp;utm_campaign=Badge_Grade)

This repository contains code for the generation (though Pickaxe), 
storage and querying of MINE databases. For generalinformation on MINE 
Databases, please consult [JJeffryes et al. 2015](http://jcheminf.springeropen.com/articles/10.1186/s13321-015-0087-1)
APIs can found in [the API Repository](https://github.com/JamesJeffryes/MINE-API)

## Repository Structure
This repository primarily consists of the minedatabases python module 
and it's 5 submodules:
- compound_io: Contains functions to load and export chemical structures 
from MINE databases. Has command-line interface.
- databases: Contains functions which impose a schema on the underlying 
Mongo databases which store MINE data.
- pickaxe: Allows for the application of reaction rules to compound sets 
and the annotation of the resulting compounds an reactions. Has command-line interface.
- queries: Contains logic for text and chemical structure queries of the 
MINE database.
- utils: Various utility functions such as hashing & type conversions

### Compound_io command-line usage
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
- tsv: FOR EXPORT ONLY, a tab separated file compatible with ModelSEED

### Pickaxe command-line usage
Pickaxe.py can be called independently to generate predictions with or 
without database storage. To list all options call `python pickaxe.py -h`. 
To predict all chemical damage reactions for one generation on compounds in the iAF1260 
model one would call `python pickaxe.py -C ./data/ChemicalDamageCoreactants.tsv -r 
./data/ChemicalDamageRxnRules.tsv -g 1 -c ./data/iAF1260.tsv`
