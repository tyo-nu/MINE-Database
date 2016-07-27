import seaborn
import pandas
import matplotlib.pyplot as plt
from databases import MINE
from collections import defaultdict, Counter
import sys


def make_violin_plots(db_list):
    # ToDo: Make properties editable
    df = pandas.DataFrame()
    for db_name in db_list:
        db = MINE(db_name)
        l = []
        cursor = db.compounds.find(projection={'_id': 0, 'Mass': 1, 'logP': 1, 'NP_likeness': 1})
        for x in cursor:
            x['DB'] = str(db_name.strip('exp2'))
            if "exp" in db_name:
                x['Type'] = 'MINE'
            else:
                x['Type'] = 'Source'
            l.append(x)
        df = df.append(l)
    df['logP'] = pandas.to_numeric(df['logP'])
    f, ax = plt.subplots(1, 3)
    seaborn.violinplot(split=True, hue='Type', x='DB', y='Mass', data=df, ax=ax[0])
    seaborn.violinplot(split=True, hue='Type', x='DB', y='logP', data=df, ax=ax[1])
    seaborn.violinplot(split=True, hue='Type', x='DB', y='NP_likeness', data=df, ax=ax[2])
    ax[1].legend_.remove()
    ax[2].legend_.remove()
    plt.tight_layout()
    plt.savefig("MINE property comparison.png")


def make_fp_heatmap(db_name, fp_type='MACCS', n_rows=25):
    db = MINE(db_name)
    data = defaultdict(Counter)
    for comp in db.compounds.find({}, {"_id": 0, "Generation": 1, fp_type: 1}):
        if fp_type in comp:
            data[int(comp['Generation'])].update(comp[fp_type])
    df = pandas.DataFrame(data)
    df_norm = df.div(df.max(axis=0), axis=1)
    if not n_rows:
        df_top = df_norm
    else:
        df_norm['range'] = df_norm.max(axis=1) - df_norm.min(axis=1)
        df_top = df_norm.sort_values('range', ascending=False).head(int(n_rows)).ix[:, :-1]
    seaborn.heatmap(df_top)
    plt.xlabel('Generation')
    plt.ylabel(fp_type + " bit")
    plt.yticks(rotation=0)
    plt.savefig(db_name + '_fp_heatmap.png')

if __name__ == "__main__":
    if sys.argv[1] == "violin":
        make_violin_plots(sys.argv[2:])

    if sys.argv[1] == 'heatmap':
        make_fp_heatmap(*sys.argv[2:])
