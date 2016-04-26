import databases
from rdkit.Chem import AllChem
import sys

if not len(sys.argv) == 3:
    sys.exit("Usage: ImportSDF.py <database name> <sdf>")
db = databases.MINE(sys.argv[1])
if not '.sdf' in sys.argv[2]:
    sys.exit("Specified file must be SDF")
suppl = AllChem.SDMolSupplier(sys.argv[2])
for mol in suppl:
    if mol:
        try:
            db.insert_compound(mol, compound_dict=mol.GetPropsAsDict(), pubchem_db=None, kegg_db=None, modelseed_db=None)
        except ValueError:
            continue