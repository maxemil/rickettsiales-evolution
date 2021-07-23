import ete3
import re
import glob
import os
from Bio import SeqIO
from collections import defaultdict

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

# "virB1":"COG0741",
# "virB2":"COG3838",
# "virB6":"COG3704",
gene2cog = {"virB3":"COG3702",
"virB4":"COG3451",
"virB8":"COG3736",
"virB9":"COG3504",
"virB10":"COG2948",
"virB11":"COG0630",
"virD4":"COG3505"}

def get_taxonomy(seqid):
    lin = ncbi.get_lineage(int(seqid.split('.')[0]))
    tax = []
    for t in lin: 
        name = ncbi.get_taxid_translator([t])[t]
        name = re.sub(r"[\ |/|\.|'|&|(|)]|:", '_', name)
        if ncbi.get_rank([t])[t] in ['species', 'family', 'class', 'phylum', 'kingdom', 'superkingdom']:
            tax.append(name)
    tax = "_".join(tax + [seqid.split('.')[0]])
    return tax

taxa = set()
gene2tax = defaultdict(set)

for gene, cog in gene2cog.items():
    for rec in SeqIO.parse("../eggnog_trees/{}_noRick.faa".format(cog), 'fasta'):
# for hits in glob.glob("../eggnog_trees/vir*_EggNOG.hits"):
#     gene = os.path.basename(hits).replace('_EggNOG.hits', '')
#     for rec in SeqIO.parse(hits, 'fasta'):
        tax = get_taxonomy(rec.id)
        taxa.add(tax)
        gene2tax[gene].add(tax)


for t in sorted(list(taxa)):
    count = 0
    for k,v in gene2tax.items():
        if t in v:
            count += 1
    if count >= 6:
        print(t)


taxa_selection = [366649, 190486, 391008, 287, 272624, 491916, 311403, 283165, 
                    354242, 435832, 292, 69395, 1095747, 297246, 1122165, 
                    1357275, 306263, 13689, 283166, 573, 29486]

# matrix = defaultdict(lambda: {t:None for t in taxa_selection})

for gene, cog in gene2cog.items():
    with open("{}.fasta".format(gene), 'w') as out:
        for rec in SeqIO.parse("../eggnog_trees/{}_noRick.faa".format(cog), 'fasta'):
            taxid = int(rec.id.split('.')[0])
            if taxid in taxa_selection:
                new_id = get_taxonomy(rec.id) + "@" + rec.id.split('.')[1]
                rec.id = new_id
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')
            # if not matrix[gene][taxid]:
            #     matrix[gene][taxid] = rec
            # elif len(matrix[gene][taxid].seq) < len(rec.seq):
            #     matrix[gene][taxid] = rec
            # print(gene, get_taxonomy(rec.id))
            
for gene, seqs in matrix.items():
    with open("{}.fasta".format(gene), 'w') as out:
        for taxid, rec in seqs.items():
            if rec:
                new_id = get_taxonomy(rec.id) + "@" + rec.id.split('.')[1]
                rec.id = new_id
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')

specs = set()
for faa in glob.glob('*Ricks.faa'):
    for rec in SeqIO.parse(faa, 'fasta'):
        specs.add(rec.id.split('@')[0])

matrix = defaultdict(lambda: {t:None for t in specs})

for faa in glob.glob('*Ricks.faa'):
    gene = faa.replace("_Ricks.faa", '')
    for rec in SeqIO.parse(faa, 'fasta'):
        spec = rec.id.split('@')[0]
        if not matrix[gene][spec]:
            matrix[gene][spec] = rec
        elif len(matrix[gene][spec].seq) < len(rec.seq):
            matrix[gene][spec] = rec    

for gene, seqs in matrix.items():
    with open("{}_Ricks_clean.fasta".format(gene), 'w') as out:
        for taxid, rec in seqs.items():
            if rec:
                SeqIO.write(rec, out, 'fasta')

# for f in vir*.fasta; do mafft-einsi --thread 5 $f > "${f%.fasta}".mafft; done
# for f in vir*.mafft; do trimal -in $f -out "${f%%.*}".95.aln -gt 0.05; done
# for f in vir*.aln; do iqtree2 -s $f -pre "${f%%.*}" -bb 1000 -m LG+G+F; done
for gene in gene2cog.keys():
    if os.path.isfile("first_round_new/{}.remove".format(gene)):
        remove_seqs = [line.strip() for line in open("first_round_new/{}.remove".format(gene))]
    else:
        remove_seqs = []
    with open('cleaned_single_genes_new/{}.fasta'.format(gene), 'w') as out:
        for rec in SeqIO.parse('first_round_new/{}.fasta'.format(gene), 'fasta'):
            if not rec.id in remove_seqs:
                SeqIO.write(rec, out, 'fasta')
                
specs = set()
for f in glob.glob('first_round/*.mafft'):
    for rec in SeqIO.parse(f, 'fasta'):
        specs.add(rec.id.split('@')[0])

for f in glob.glob('cleaned_single_genes/*.fasta'):
    with open(f.replace('.fasta', '.clean_fasta'), 'w') as out:
        for rec in SeqIO.parse(f, 'fasta'):
            for sp in specs:
                if rec.id.startswith("{}_".format(sp)):
                    rec.id = rec.id.replace("{}_".format(sp), "{}@".format(sp))
                    rec.description = ""
                elif rec.id.startswith("{}..".format(sp)):
                    rec.id = rec.id.replace("{}..".format(sp), "{}@".format(sp))
                    rec.description = ""
            SeqIO.write(rec, out, 'fasta')

for f in glob.glob('cleaned_single_genes/*.fasta'):
    spec_counts = defaultdict(int)
    for rec in SeqIO.parse(f, 'fasta'):
        spec_counts[rec.id.split('@')[0]] += 1
    for sp, count in spec_counts.items():
        if count > 1:
            print(f, sp, count)

for gene in gene2cog.keys():
    if os.path.isfile("cleaned_single_genes/{}.remove".format(gene)):
        remove_seqs = [line.strip() for line in open("cleaned_single_genes/{}.remove".format(gene))]
    else:
        remove_seqs = []
    with open('cleaned_single_genes/{}.clean_dedup_fasta'.format(gene), 'w') as out:
        for rec in SeqIO.parse('cleaned_single_genes/{}.clean_fasta'.format(gene), 'fasta'):
            if not rec.id in remove_seqs:
                SeqIO.write(rec, out, 'fasta')
