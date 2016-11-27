#MINE Databases

This repository contains code for the generation (though Pickaxe), storage and querying of MINE databases. For general
information on MINE Databases, please consult [JJeffryes et al. 2015](http://jcheminf.springeropen.com/articles/10.1186/s13321-015-0087-1)
APIs can found in [the API Repository](https://github.com/JamesJeffryes/MINE-API)

###Pickaxe command line usage
Pickaxe.py can be called independently to generate predictions with or without database storage. To list all options
call `python pickaxe.py -h`. To predict all chemical damage reactions on compounds in the iAF1260 model one would call
`python -C ./data/ChemicalDamageCofactors.tsv -r ./data/ChemicalDamageOperators.tsv -g 1 -c ./data/iAF1260.tsv`