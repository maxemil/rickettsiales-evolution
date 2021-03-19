import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from collections import defaultdict
from itertools import product
import math
import glob
import ete3

def parse_ALE_events(suffix):
    df = pd.read_csv('marineRick-aphaNOG-{}/events.txt'.format(suffix), sep='\t', header=None)
    return df

def subset_events(df, clusname, event):
    if '{}.clean.ufboot.ale'.format(clusname) in set(df[1]):
        clus = df[df[1] == '{}.clean.ufboot.ale'.format(clusname)]
    else:
        clus = df[df[1] == '{}'.format(clusname)]
        for n in set(df[2]):
            if not n in set(clus[2]):
                clus = clus.append({2:n}, ignore_index=True)
        clus = clus.fillna(0)
    clus.index = clus[2]
    clus = clus.drop([0,1,2], axis=1)
    clus.columns = ['duplications', 'transfers', 'losses', 'originations', 'copies']
    clus[clusname] = clus[event].apply(lambda x: math.modf(x)[1] + int(math.modf(x)[0] >= 0.3))
    return clus

def parse_apronogs(suffix, apronoglist, event, additional=[]):
    df = parse_ALE_events(suffix)
    types = ['duplications', 'transfers', 'losses', 'originations']
    node2event = defaultdict(lambda: {t:0 for t in types})
    all_apronogs = pd.DataFrame()
    for line in open(apronoglist):
        apronog = line.strip()
        clus = subset_events(df, apronog, event)
        all_apronogs[apronog] = clus[apronog]
    for apronog in additional:
        clus = subset_events(df, apronog, event)
        all_apronogs[apronog] = clus[apronog]
    return all_apronogs
    
def make_heatmap(all_apronogs, outfile): 
    N = all_apronogs.shape[0]
    M = all_apronogs.shape[1]
    ylabels = all_apronogs.index
    xlabels = all_apronogs.columns
    ap = all_apronogs.to_numpy()
    ap = (ap > 0).astype(int)
    x, y = np.meshgrid(np.arange(M), np.arange(N))
    fig, ax = plt.subplots(figsize=(M/7+5, N/6))
    c = np.zeros(shape=(N,M))
    circles = [plt.Circle((j,i), radius=r/2) for r, j, i in zip(ap.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=c.flatten(), facecolors=["black"] * c.size)
    ax.add_collection(col)
    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticklabels(xlabels, rotation=90, font='monospace')
    ax.set_yticklabels(ylabels, font='monospace')
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.grid(which='minor')
    # plt.show()
    plt.tight_layout()
    plt.savefig(outfile)
    plt.clf()

def inorderTraversal(node):
    res = []
    if not node.is_leaf():
        res = inorderTraversal(node.children[0])
        res.append(node.name)
        res = res + inorderTraversal(node.children[1])
    else:
        res = [node.name]
    return res

def order_df(all_apronogs, suffix):
    treefiles = {'norm':'species_trees/marineRick_b5000g15300_clean.tree', 
                 'wDeia':'species_trees/marineRick_b5000g20100_wDeia_clean.tree'}
    code2name = {line.split()[1]:line.split()[0] for line in open('species_map.txt')}
    code2inners = {'norm':{'90':'LRiCA', '76':'LMiCA', '79':'LAnCA', 
                           '91':'LRssCA', '93':'LRslCA', '94':'root', 
                           '89':'LoutCA', '86':'L125CA', '64':'LubaCA'},
                   'wDeia':{'100':'LRiCA', '86':'LMiCA', '90':'LAnCA', 
                            '101':'LRssCA', '103':'LRslCA', '85':'LDeCA',
                            '104':'root', '96':'LoutCA', '98':'L125CA', 
                            '70':'LubaCA'}}
    code2name.update(code2inners[suffix])
    tree = ete3.PhyloTree(treefiles[suffix], format=1)
    tree.ladderize(direction=1)
    sorted_nodes = inorderTraversal(tree)
    all_apronogs = all_apronogs.loc[sorted_nodes]
    all_apronogs.index = pd.Series(all_apronogs.index).apply(lambda x: code2name[x] if x in code2name else x)
    return all_apronogs

def get_color_map(suffix):
    treefiles = {'norm':'species_trees/marineRick_b5000g15300_clean.tree', 
                 'wDeia':'species_trees/marineRick_b5000g20100_wDeia_clean.tree'}
    code2name = {line.split()[1]:line.split()[0] for line in open('species_map.txt')}
    code2inners = {'norm':{'90':'LRiCA', '76':'LMiCA', '79':'LAnCA', 
                           '91':'LRssCA', '93':'LRslCA', '94':'root', 
                           '89':'LoutCA', '86':'L125CA', '64':'LubaCA'},
                   'wDeia':{'100':'LRiCA', '86':'LMiCA', '90':'LAnCA', 
                            '101':'LRssCA', '103':'LRslCA', '85':'LDeCA',
                            '104':'root', '96':'LoutCA', '98':'L125CA', 
                            '70':'LubaCA'}}
    color_map = defaultdict(lambda: 'black')
    code2name.update(code2inners[suffix])
    tree = ete3.PhyloTree(treefiles[suffix], format=1)
    for n in tree.traverse():
        if n.name in code2name:
            n.name = code2name[n.name]
    ancestor2color = {'LRiCA':'goldenrod', 'LMiCA':'olive', 'LAnCA':'chocolate',
                      'LDeCA':'salmon', 'LoutCA':'grey', 'L125CA':'cyan', 
                      'LubaCA':'teal'}
    for anc, col in ancestor2color.items():
        try:
            for n in (tree & anc).traverse():
                color_map[n.name] = col
        except:
            print(anc)
    return color_map

def seaborn_heatmap(all_apronogs, all_originations, all_losses, color_map, nog2gene, outfile):
    all_apronogs.columns = pd.Series(all_apronogs.columns).apply(lambda x: nog2gene[x] if x in nog2gene else x)
    all_originations.columns = pd.Series(all_originations.columns).apply(lambda x: nog2gene[x] if x in nog2gene else x)
    all_losses.columns = pd.Series(all_losses.columns).apply(lambda x: nog2gene[x] if x in nog2gene else x)
    fig, ax = plt.subplots(figsize=(all_apronogs.shape[1]/7+5, 
                                    all_apronogs.shape[0]/6))
    sns.heatmap((all_apronogs > 0), xticklabels=True, yticklabels=True, 
                        cmap=sns.light_palette("seagreen", as_cmap=True),
                        cbar=False, ax=ax)
    sns.heatmap((all_originations > 0), xticklabels=True, yticklabels=True, 
                        cmap=["blue"], mask=(all_originations < 1), cbar=False,
                        alpha = 0.5, ax=ax)
    sns.heatmap((all_losses > 0), xticklabels=True, yticklabels=True, 
                        cmap=["red"], mask=(all_losses < 1), cbar=False,
                        alpha = 0.5, ax=ax)
    # sns.heatmap((all_transfers > 0), xticklabels=True, yticklabels=True, 
    #                     cmap=["orange"], mask=(all_transfers < 1), cbar=False,
    #                     alpha = 0.5, ax=ax)
    ax.set_xticklabels(labels=ax.get_xticklabels(), font='monospace', rotation=90)
    ax.set_yticklabels(labels=ax.get_yticklabels(), font='monospace')
    ax.set_ylabel('')
    for tickl in ax.get_yticklabels():
        tickl.set_color(color_map[tickl.get_text()])
    plt.tight_layout()
    plt.savefig("heatmaps_GoI/" + outfile)
    plt.close()
                            
def main():
    topos = ['norm', 'wDeia']
    for pathway, suffix in product(glob.glob("visualize_GoIs/T4SS*.txt"), topos):
        prefix = pathway.replace('visualize_GoIs/', '').replace('.txt', '')
        all_apronogs = parse_apronogs(suffix, pathway, 'copies')
        all_apronogs = order_df(all_apronogs, suffix)
        all_originations = parse_apronogs(suffix, pathway, 'originations')
        all_originations = order_df(all_originations, suffix)
        all_losses = parse_apronogs(suffix, pathway, 'losses')
        all_losses = order_df(all_losses, suffix)
        color_map = get_color_map(suffix)
        nog2gene = {line.split()[0]:line.strip().split()[1] for line in open(pathway.replace('txt', 'map'))}
        # all_transfers = parse_apronogs(suffix, pathway, 'transfers')
        # all_transfers = order_df(all_transfers, suffix)
        seaborn_heatmap(all_apronogs, all_originations, all_losses, color_map, nog2gene, "{}_{}.pdf".format(prefix, suffix))
        # make_heatmap(all_apronogs, "{}_{}.pdf".format(prefix, suffix))
            
#######################

from Bio import SeqIO
import glob
import os

def get_unlab(p):
    for f in glob.glob('../29_clusters_to_trees/faa_with_Deianiraea/small_cluster/unlab_60_9*') + \
            glob.glob('../29_clusters_to_trees/faa_with_Deianiraea/unlab_60_9*'):
        for rec in SeqIO.parse(f, 'fasta'):
            if p in rec.id:
                return os.path.basename(f).replace('.faa', '')

df = pd.read_csv('../33_metabolic_reconstruction/all_annotations_marineRick_secret4.csv', sep='\t')
pf = df[(df['Pfam'] == 'PF04956') & (df['alphaNOG'].isna())]
additionals = set(pf['Sequence ID'].apply(lambda x: get_unlab(x)).dropna())

pathway = 'visualize_GoIs/T4SS.txt'
prefix = pathway.replace('visualize_GoIs/', '').replace('.txt', '')
suffix = 'wDeia'

all_apronogs = parse_apronogs(suffix, pathway, 'copies', additionals)
all_apronogs = order_df(all_apronogs, suffix)

nog2gene = {line.split()[0]:line.strip().split()[1] for line in open(pathway.replace('txt', 'map'))}
for x in additionals:
    nog2gene[x] = 'virB2'
all_apronogs.columns = pd.Series(all_apronogs.columns).apply(lambda x: nog2gene[x] if x in nog2gene else x)


summed_T4SS = all_apronogs.T.groupby([i for i in all_apronogs.T.index.values]).sum().T
fig, ax = plt.subplots(figsize=(all_apronogs.shape[1]/7+5, 
                                all_apronogs.shape[0]/6))
sns.heatmap(summed_T4SS, annot=True, 
            cmap=sns.color_palette("crest", as_cmap=True), 
            xticklabels=True, yticklabels=True, ax=ax)
plt.tight_layout()
plt.savefig("T4SS_subunits.pdf")
