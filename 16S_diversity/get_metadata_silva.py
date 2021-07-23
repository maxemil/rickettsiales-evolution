from Bio import Entrez
from Bio import SeqIO


Entrez.email = 'max-emil.schon@icm.uu.se'


keywords = ['marine', 'Marine', 'seawater', 'pacific', 'atlantic', 'ocean', 'aquatic']
# accessions = []
with open('67_RIFL_SILVA_hits_metadata.fasta', 'w') as out:
    for rec in SeqIO.parse('67_RIFL_SILVA_hits.fasta', 'fasta'):
        # accessions.append(rec.id.split('.')[0])
        acc = rec.id.split('.')[0]
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode='text')
        gb = handle.read()
        if any([kw in gb for kw in keywords]):
            rec.id = "{}_{}".format('marine', rec.description.replace(' ', '_'))
        else:
            rec.id = "{}_{}".format('unknown', rec.description.replace(' ', '_'))
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
