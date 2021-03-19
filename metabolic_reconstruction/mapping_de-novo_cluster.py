from Bio import SeqIO
import glob
import os

with open('mapping_de-novo_cluster.csv', 'w') as out:
    for f in glob.glob('faa_with_Deianiraea/small_cluster/[!0]*.faa') + glob.glob('faa_with_Deianiraea/[!0]*.faa'):
        clst = os.path.basename(f).replace('.faa', '')
        for rec in SeqIO.parse(f, 'fasta'):
            print(clst, rec.id, sep='\t', file=out)
