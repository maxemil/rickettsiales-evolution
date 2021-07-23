from Bio import Entrez
from Bio import SeqIO
import pandas as pd

Entrez.email = 'max-emil.schon@icm.uu.se'

prefix = "NR_hits_RIFL"

df = pd.read_csv('{}.csv'.format(prefix), sep=',', header=None)
df = df[df[2] > 90]
df = df[df[3] > 500]
hits = set(df[1])

marine_keywords = ['marine', 'Marine', 'seawater', 'Seawater', 'pacific', 'Pacific',
            'atlantic', 'ocean', 'Ocean', 'China Sea', 'Gulf of Mexico',
            'sea water', 'GS22', 'GS19', 'GS08', 'Yellow Sea']
freshwater_keywords = ['freshwater', 'Lake', 'lake', 'GS20', 'dam reservoir',
            'spring water', 'cenote', 'Desert Springs', 'karst water', 
            'glacier ice']

hit_recs = []
for rec in SeqIO.parse('{}.fasta'.format(prefix), 'fasta'):
    if rec.id in hits:
        hit_recs.append(rec)
        
with open('{}_metadata.fasta'.format(prefix), 'w') as out:
    for i in range(0, len(hit_recs), 200):
        acc_fetch = ','.join([r.id.split('.')[0] for r in hit_recs[i:i+200]])
        handle = Entrez.efetch(db="nucleotide", id=acc_fetch, rettype="gb", retmode='text')
        gbs = handle.read().split('//')
        for gb, rec in zip(gbs,hit_recs[i:i+200]):
            if any([kw in gb for kw in marine_keywords]):
                rec.id = "{}_{}".format('marine', rec.description.replace(' ', '_'))
            elif any([kw in gb for kw in freshwater_keywords]):
                rec.id = "{}_{}".format('freshwater', rec.description.replace(' ', '_'))
            else:
                rec.id = "{}_{}".format('unknown', rec.description.replace(' ', '_'))
                print("unknown environment {}".format(rec.id))
            rec.description = ""
            SeqIO.write(rec, out, 'fasta')
# acc = ','.join(accessions)
# handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode='text')
# for rec in SeqIO.parse(handle, 'genbank'):
#         records.append(rec)

recs = {}
for rec in SeqIO.parse('67_RIFL_SILVA_hits_metadata.fasta', 'fasta'):
    recs[rec.id] = rec
for rec in SeqIO.parse('125_SILVA_hits_metadata.fasta', 'fasta'):
    recs[rec.id] = rec
with open('125_67_RIFL_SILVA_hits_metadata.fasta', 'w') as out:
    for rec in recs.values():
        SeqIO.write(rec, out, 'fasta')
