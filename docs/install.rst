Installation
============

MINE-Database requires the use of rdkit, which currently is unavailable to install on pip. Thus, we recommend you use conda to create a new environment and then install rdkit into that environment before proceeding:

conda create -n mine

conda activate mine

conda install -c rdkit rdkit

Then, use pip (in your conda environment) to install minedatabase:

pip install minedatabase
