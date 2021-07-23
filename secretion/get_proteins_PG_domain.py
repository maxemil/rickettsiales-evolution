from Bio import SeqIO
from collections import defaultdict

prots = defaultdict(list)
for f in ['bin_125.out','TOBG_RS-372.out','UBA6149.out','UBA6178.out','UBA6189.out']:
    for line in open(f):
        if not line.startswith('#'):
            line = line.split()
            prots[line[0]].append(line[2])

recs = []
for f in ['bin_125.faa','TOBG_RS-372.faa','UBA6149.faa','UBA6178.faa','UBA6189.faa']:
    for rec in SeqIO.parse("../../33_metabolic_reconstruction/proteomes/{}".format(f), 'fasta'):
        if rec.id in prots:
            rec.description = "_".join(prots[rec.id])
            recs.append(rec)

with open('Mitibacteraceae_effector_domains.fasta', 'w') as out:
    for rec in recs:
        SeqIO.write(rec, out, 'fasta')
