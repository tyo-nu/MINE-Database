import pandas
import seaborn
import matplotlib.pyplot as plt
import numpy
from minedatabase.databases import MINE
import sys

db = MINE(sys.argv[1])
fields = ['Compounds', 'Compound_ids', 'Reactions', 'Operators']

def pw_jaccard(series, reduce=numpy.median):
    pw = []
    for i, x in enumerate(series):
        tc = []
        for j, y in enumerate(series):
            if i != j:
                tc.append(len(x & y) / float(len(x | y)))
        pw.append(reduce(tc))
    return pw

keys = {}
results = []
for model in db.models.find():
    results.append([model['_id']]+[set([y[0] for y in model[x]])
                                   if isinstance(model[x][0], list)
                                   else set(model[x]) for x in fields])

sets = pandas.DataFrame(results, columns=['_id']+fields).set_index('_id')
sets.to_csv("model_sets.csv")
tcs = sets.apply(pw_jaccard)
tcs.to_csv("model_pairwise_jaccard.csv")
results = pandas.DataFrame.from_csv("model_pairwise_jaccard.csv",
                                    index_col='_id')
seaborn.boxplot(data=results)
plt.tight_layout()
plt.savefig("%s_boxplots.png" % db.name)