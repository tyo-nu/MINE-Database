import os
import subprocess
import sys
import hashlib

import pandas as pd
from minedatabase import utils
from minedatabase.databases import MINE
from rdkit.Chem import AllChem

def load_cdmine_rxns(mine_db, excel_file, pic_dir=""):
    abrv = {"hn": "[*]"}
    if pic_dir and not os.path.exists(pic_dir):
        os.mkdir(pic_dir)
    compounds = pd.read_excel(excel_file, 1, skiprows=1).fillna("")
    reactions = pd.read_excel(excel_file, 0, skiprows=1).fillna("")

    for i, row in compounds.iterrows():
        if row['SMILES']:
            mol = AllChem.MolFromSmiles(row['SMILES'])
            if mol:
                c_id = mine_db.insert_compound(mol, {"Generation": 0})
                abrv[row['Abbreviation'].strip()] = c_id
                if pic_dir:
                    rc = subprocess.call("/Applications/ChemAxon/JChem/bin/molconvert -o %s/temp.png png:-a,w500 -s "
                                         "'%s'" % (pic_dir, row['SMILES'].strip()), shell=True)
                    if not rc:
                        os.rename(pic_dir + "temp.png", pic_dir + c_id + ".png")
            else:
                print("Failed to parse %s" % row['SMILES'])
        else:
            print('SMILES missing from %s' % row.name)

    reactions['Type of Reaction'].fillna('ffill', inplace=True)
    for i, row in reactions.iterrows():
        if row['Equation (Abbreviations)']:
            rxn = row[['Metabolite', 'Equation (full names)']].to_dict()
            if isinstance(row['PMID or doi'], str):
                rxn['References'] = row['PMID or doi'].strip().split('; ')
            else:
                rxn['References'] = [str(row['PMID or doi'])]

            rxn['Type'] = str(row['Type of Reaction']).strip()
            rxn['Notes'] = str(row['Comments']).strip()
            rxn['Reactants'], rxn['Products'] = utils.parse_text_rxn(row['Equation (Abbreviations)'], ' = ', ' + ', abrv)
            rxn['InChI_hash'] = utils._calculate_rxn_hash(mine_db, rxn['Reactants'], rxn['Products'])
            mine_db.insert_reaction(rxn)
        else:
            print('RXN missing from %s' % row.name)

if __name__ == '__main__':
    mine = MINE(sys.argv[1])
    load_cdmine_rxns(mine, sys.argv[2])