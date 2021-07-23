from Bio import SeqIO
import ete3
import re
import sys

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

def get_taxonomy(seqname):
    if re.match("^[0-9]+", seqname):
        lin = ncbi.get_lineage(int(seqname.split('.')[0]))
        tax = []
        for t in lin: 
            name = ncbi.get_taxid_translator([t])[t]
            name = re.sub(r"[\ |/|\.|'|&]", '_', name)
            if ncbi.get_rank([t])[t] in ['species', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']:
                tax.append(name)
        tax = "_".join(tax + [seqname])
        return tax

infile = sys.argv[1]
outfile = sys.argv[2]
with open(outfile, 'w') as out:
    for rec in SeqIO.parse(infile, 'fasta'):
        if not 'Rickettsiales' in get_taxonomy(rec.id):
                SeqIO.write(rec, out, 'fasta')
