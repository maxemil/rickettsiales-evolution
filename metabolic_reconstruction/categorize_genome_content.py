import pandas as pd
from collections import defaultdict
from math import modf
import ete3
import numpy as np

def parse_ALE_events(suffix):
    header = ['species_tree', 'cluster', 'node', 'duplications',
              'transfers', 'losses', 'originations', 'copies']
    df = pd.read_csv(
            '../34_ALE_updated_clusters/marineRick-aphaNOG-{}/events.txt'.format(suffix),
            sep='\t', names=header)
    del df['species_tree']
    # df = df[df['cluster'].apply(lambda x:'.clean.ufboot.ale' in x)] # remove singletons
    df['cluster'] = df['cluster'].apply(lambda x: x.replace('.clean.ufboot.ale', ''))
    return df

def get_cat_groups():
    CAT_GROUPS = {"Metabolism":['C','G','E','F','H','I','P','Q'],
                  "Cellular Processes":['D','Y','V','T','M','N','Z','W','U','O'],
                  "Information":['J','A','K','L','B'],
                  "Unknown":['S']}
    catmap = {}
    for k, value in CAT_GROUPS.items():
        for v in value:
            catmap[v] = k
    return catmap

def get_broad_group(cat, catmap):
    groups = set()
    for c in cat:
        groups.add(catmap[c])
    if len(groups) == 1:
        return groups.pop()
    else:
        return "Unknown"

def add_copies_cat(cat_node, node, bcat, copies):
    cat_node.loc[node, bcat] += copies

annot_header = ["taxid", "apronog", "protcount", "speciesCount", "category", "annotation"]
annot = pd.read_csv('aproNOG.annotations.tsv', 
        sep='\t', names=annot_header)
annot['apronog'] = annot['apronog'].apply(lambda x: x.replace('ENOG41', ''))
annot.index = annot['apronog']
cat_dict = annot['category'].to_dict()
cat_dict = defaultdict(lambda: 'S', cat_dict)

catmap = get_cat_groups()

df = parse_ALE_events('wDeia')
df['copies'] = df['copies'].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))
df['cat'] = df['cluster'].apply(lambda x: cat_dict[x])
df['bcat'] = df['cat'].apply(lambda x: get_broad_group(x, catmap))

cat_node = pd.DataFrame(index=set(df['node']), columns=set(df['bcat']))
cat_node.fillna(0, inplace=True)
df.apply(lambda row: add_copies_cat(cat_node, row['node'], row['bcat'], row['copies']), axis=1)


def layout(node):
    if not node.is_leaf():
        node.img_style['size'] = 0
    cols = ['Metabolism', 'Cellular Processes', 'Information', 'Unknown']
    percents = cat_node.loc[node.name, cols]/cat_node.loc[node.name, cols].sum() * 100
    pie = N = ete3.PieChartFace(percents, 
                    np.sqrt(df.loc[df['node'] == node.name, 'copies'].sum()/np.pi)/2, 
                    np.sqrt(df.loc[df['node'] == node.name, 'copies'].sum()/np.pi)/2, 
                    colors=['coral', 'yellowgreen', 'deepskyeblue', 'grey'])
    ete3.faces.add_face_to_node(pie, node, column=0, position='float')
    
t = ete3.Tree("../34_ALE_updated_clusters/species_trees/marineRick_b5000g20100_wDeia_clean.tree", format=1)

ts = ete3.TreeStyle()
ts.layout_fn = layout
ts.show_branch_support = True
t.ladderize(direction=1)
t.render('test_tree.pdf', tree_style=ts)
