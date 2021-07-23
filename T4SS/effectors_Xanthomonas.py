# from Bio import SeqIO
# 
# xac_loc = {'XAC0466':487,
# 'XAC3266':734,
# 'XAC2609':314,
# 'XAC0096':505,
# 'XAC1918':476,
# 'XAC0574':317,
# 'XAC3634':189,
# 'XAC0151':120,
# 'XAC1165':0,
# 'XAC0323':18,
# 'XAC2885':271,
# 'XAC4264':164}
# 
# with open("XVIPCD.faa", 'w') as out:
#     for rec in SeqIO.parse('XAC_seqs.gp', 'genbank'):
#         for f in rec.features:
#             if 'locus_tag' in f.qualifiers:
#                 rec.id = f.qualifiers['locus_tag'][0]
#                 rec.seq = rec.seq[xac_loc[rec.id]:]
#                 SeqIO.write(rec, out, 'fasta')    
# 


import glob
from collections import defaultdict
import pandas as pd

sp2domain = defaultdict(lambda: {'Peptidase_M23':0, 'PG_binding_1':0, 'DUF2974':0, 'Glyco_hydro_19':0})
for f in glob.glob('*.out'):
    sp = f.replace('.out', '')
    for line in open(f):
        if not line.startswith('#'):
            line = line.split()
            sp2domain[sp][line[2]] += 1

df = pd.DataFrame.from_dict(sp2domain).T
df = df.loc[[sp for k,g in groups.items() for sp in g]]
df.to_csv('effector_domains_summary.csv', sep='\t', header=True, index=True)
