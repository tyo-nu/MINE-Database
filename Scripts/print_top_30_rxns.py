from minedatabase.compound_io import export_inchi_rxns
from minedatabase.databases import MINE

top_30 = [
    'cpd00229',
    'cpd00146',
    'cpd00084',
    'cpd00216',
    'cpd00040',
    'cpd00219',
    'cpd02961',
    'cpd02978',
    'cpd00095',
    'cpd00016',
    'cpd00858',
    'cpd00056',
    'cpd00005',
    'cpd00004',
    'cpd00982',
    'cpd00015',
    'cpd01777',
    'cpd01775',
    'cpd00103',
    'cpd00932',
    'cpd00931',
    'cpd00957',
    'cpd02642',
    'cpd02431',
    'cpd00236',
    'cpd03470',
    'cpd02826',
    'cpd00042',
    'cpd00010',
    'cpd00655',
    'cpd02666',
    'cpd00834',
    'cpd00135',
    'cpd00346',
    'cpd00020',
    'cpd00956',
    'cpd00330',
    'cpd00087',
    'cpd00003',
    'cpd00006',
    'cpd00102',
    'cpd00032',
    'cpd00070',
    'cpd00918',
    'cpd00203',
    'cpd02857',
    'cpd00031',
    'cpd00038',
    'cpd00126',
    'cpd00241',
    'cpd00295',
    'cpd02552',
    'cpd00338',
    'cpd00683',
    'cpd00171',
    'cpd00198',
    'cpd00238',
    'cpd01977',
    'cpd00051',
    'cpd02069',
]

db = MINE('plant_spontanious')
rxn_ids = set()
for cpd_id in top_30:
    cpd = db.compounds.find_one({"DB_links":{'Model_SEED': cpd_id}}, {'Reactant_in': 1})
    if cpd:
        print(cpd_id)
        rxn_ids.update(cpd.get('Reactant_in'))
    else:
        print("Can't find: {}".format(cpd_id))
print("Printing {} rxns".format(len(rxn_ids)))
export_inchi_rxns(db, "./", list(rxn_ids))
