import pandas
import os
import re
from collections import deque

# translate coreactant names into cpd_ids
coreactants = {'Any': 'Any'}
for line in open('../minedatabase/data/EnzymaticCoreactants.tsv'):
    if line[0] == "#":
        continue
    sl = line.split('\t')
    coreactants[sl[1]] = sl[0]


def name_to_cid(namelist, delim=';'):
    try:
        return delim.join([coreactants[x] for x in namelist.split(delim)])
    except KeyError:
        print("Could not parse %s" % namelist)


def get_products(op_dir):
    """Get product and comment data from operator file"""
    op_data = {'Comments': {}, 'Products': {}}
    for path, dir, names in os.walk(op_dir):
        for filename in names:
            if '.dat' in filename:
                bnice_op = open(path + "/" + filename).read()
                try:
                    tail = re.split('Products:\s*\n', bnice_op)[1]
                    products, comments = re.split('Comments:\s*\n?', tail)
                    op_data['Comments'][filename[:-4]] = comments.strip()
                    op_data['Products'][filename[:-4]] = ";".join(products.strip().split())
                except (IndexError, ValueError):
                    print(filename)
    return op_data


def rotate_products(products, sep=';'):
    dq = deque(products.split(sep))
    dq.rotate(1)
    return sep.join(dq)

operators = pandas.read_csv('../minedatabase/data/EnzymaticReactionRules.tsv',
                            delimiter='\t').set_index('Name')
error_ops = set(re.findall("Warning: Unbalanced Reaction produced by "
                           "(\d.\d+.-?\d+.\w)", open('/Applications/PyCharm.app/Contents/bin/error.txt').read()))
print("%s operators will be rotated" % len(error_ops))
print(error_ops)
operators.update(operators.filter(items=error_ops, axis='index')['Products'].apply(rotate_products))
#operators["Reactants"] = operators["Reactants"].apply(name_to_cid)
#operators["Products"] = operators["Products"].apply(name_to_cid)
operators.to_csv('../minedatabase/data/EnzymaticReactionRules.tsv', sep='\t')
