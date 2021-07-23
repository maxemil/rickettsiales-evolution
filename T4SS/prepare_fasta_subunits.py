import pandas as pd
from collections import defaultdict
import glob
from Bio import SeqIO
import os
    
def get_unlab(p):
    for f in glob.glob('../../29_clusters_to_trees/faa_with_Deianiraea/small_cluster/unlab_60_9*') + \
            glob.glob('../../29_clusters_to_trees/faa_with_Deianiraea/unlab_60_9*'):
        for rec in SeqIO.parse(f, 'fasta'):
            if p in rec.id:
                return os.path.basename(f).replace('.faa', '')

df = pd.read_csv('../../33_metabolic_reconstruction/all_annotations_marineRick_secret4.csv', sep='\t')
pf = df[(df['Pfam'] == 'PF04956') & (df['alphaNOG'].isna())]
additionals = set(pf['Sequence ID'].apply(lambda x: get_unlab(x)).dropna())

pathway = '../../34_ALE_updated_clusters/visualize_GoIs/T4SS.txt'

nog2gene = {line.split()[0]:line.strip().split()[1] for line in open(pathway.replace('txt', 'map'))}
for x in additionals:
    nog2gene[x] = 'virB2'
   
gene2nogs = defaultdict(list)
for k,v in nog2gene.items():
    gene2nogs[v].append(k)

gene2nogs['virB4'].append('01UN7') 

gene2nogs = {'TLC':['01SSE','026SU','0268T','025QK','026F2','01W9K']}

for gene, nogs in gene2nogs.items():
    with open('{}.faa'.format(gene), 'w') as out:
        for nog in nogs:
            for rec in SeqIO.parse(glob.glob('../../29_clusters_to_trees/faa_with_Deianiraea/**/{}.faa'.format(nog), recursive=True)[0], 'fasta'):
                rec.id = "{}_{}".format(rec.id, nog)
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')

os.remove("virB11?.faa")
