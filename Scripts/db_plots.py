import matplotlib
matplotlib.use('Agg')
import seaborn
import pandas
import matplotlib.pyplot as plt
from minedatabase.databases import MINE
from collections import defaultdict, Counter
import sys


def make_violin_plots(db_list, prop_list=('Mass', 'logP', 'NP_likeness')):
    df = pandas.DataFrame()
    for db_name in db_list:
        db = MINE(db_name)
        l = []
        cursor = db.compounds.find({"Type":{'$ne': 'Coreactant'}},
                                   dict([('_id', 0), ('Type', 1)]
                                        + [(x, 1) for x in prop_list]))
        for x in cursor:
            x['DB'] = str(db_name.strip('exp2'))
            l.append(x)
        df = df.append(l)
    f, ax = plt.subplots(1, len(prop_list))
    for i, prop in enumerate(prop_list):
        seaborn.violinplot(split=True, hue='Type', x='DB', y=prop, data=df,
                           ax=ax[i])
        if i > 0:
            ax[i].legend_.remove()
    plt.tight_layout()
    plt.savefig("MINE property comparison.png")


def make_box_plots(db_list, prop_list=('Mass', 'logP', 'NP_likeness')):
        df = pandas.DataFrame()
        for db_name in db_list:
            db = MINE(db_name)
            new_name = str(db_name.replace('exp2', 'MINE').split('-')[0])
            l = []
            cursor = db.compounds.find(dict([(x, {'$exists': 1})
                                             for x in prop_list]),
                                       dict([('_id', 0)]
                                            + [(x, 1) for x in prop_list]))
            for x in cursor:
                x['DB'] = new_name
                l.append(x)
            df = df.append(l)
        f, ax = plt.subplots(1, len(prop_list))
        for i, prop in enumerate(prop_list):
            seaborn.boxplot(x='DB', y=prop, data=df, ax=ax[i],
                            showfliers=False)
        plt.tight_layout()
        plt.savefig("MINE property comparison.png")


def make_fp_heatmap(db_name, fp_type='MACCS', n_rows=25):
    db = MINE(db_name)
    data = defaultdict(Counter)
    for comp in db.compounds.find({}, {"_id": 0, "Generation": 1, fp_type: 1}):
        if fp_type in comp and int(comp['Generation']) > -1:
            data[int(comp['Generation'])].update(comp[fp_type])
    df = pandas.DataFrame(data)
    df_norm = df.div(df.max(axis=0), axis=1)
    if not n_rows:
        df_top = df_norm
    else:
        df_norm['range'] = df_norm.max(axis=1) - df_norm.min(axis=1)
        df_top = df_norm.sort_values('range', ascending=False).head(int(n_rows)).ix[:, :-1]
    hm = seaborn.heatmap(df_top)
    hm.collections[0].colorbar.set_label("Prevalence")
    plt.xlabel('Generation')
    plt.ylabel(fp_type + " bit")
    plt.yticks(rotation=0)
    plt.savefig(db_name + '_fp_heatmap.png')

if __name__ == "__main__":
    if sys.argv[1] == "violin":
        make_violin_plots(sys.argv[2:])

    elif sys.argv[1] == "boxplot":
        make_box_plots(sys.argv[2:])

    elif sys.argv[1] == 'heatmap':
        make_fp_heatmap(*sys.argv[2:])

    else:
        print("Unrecognised plot type")
