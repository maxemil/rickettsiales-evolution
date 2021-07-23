cat 125_67_RIFL_SILVA_hits_metadata.fasta 16S_genomes.fasta > SILVA_GENOMES.fasta

cat /local/two/Exchange/projects/marineRick/35_16S_trees/NR_hits_125_metadata.99.fasta /local/two/Exchange/projects/marineRick/35_16S_trees/NR_hits_67-1_metadata.fasta /local/two/Exchange/projects/marineRick/35_16S_trees/NR_hits_RIFL_metadata.fasta > NR_hits_metadata.fasta

from Bio import SeqIO
seqdict = {}
for rec in SeqIO.parse('NR_hits_metadata.fasta', 'fasta'):
  seqdict[rec.id] = rec
for rec in SeqIO.parse('16S_genomes.fasta', 'fasta'):
  seqdict[rec.id] = rec
with open('NR_GENOMES.fasta', 'w') as out:
  for rec in seqdict.values():
    SeqIO.write(rec, out, 'fasta')

# cat NR_hits_metadata.fasta 16S_genomes.fasta > NR_GENOMES.fasta

mafft-einsi --thread 20 --adjustdirection NR_GENOMES_noAlpha.fasta > NR_GENOMES_noAlpha.aln
sed -i 's/_R_//g' NR_GENOMES_noAlpha.aln
trimal -in NR_GENOMES_noAlpha.aln -out NR_GENOMES_noAlpha.99.aln -gt 0.01

iqtree2 -s NR_GENOMES_noAlpha.99.aln -pre NR_GENOMES_noAlpha -nt 20 -m GTR+F+R8 -b 100 -redo


sed -i 's/unknown_MEMO01000010/rifle_MEMO01000010/' SILVA_GENOMES.treefile
sed -i 's/unknown_MEMO01000010/rifle_MEMO01000010/' SILVA_GENOMES.treefile
sed -i 's/unknown_FPLS01007638/soil_FPLS01007638/' SILVA_GENOMES.treefile
sed -i 's/unknown_FPLK01002650/freshwater_FPLK01002650/' SILVA_GENOMES.treefile
sed -i 's/unknown_FPLL01006513/soil_FPLL01006513/' SILVA_GENOMES.treefile

python /local/two/Software/misc-scripts/color_tree.py -i NR_GENOMES_noAlpha.treefile -o NR_GENOMES_noAlpha.nex -c taxa_colors.txt


sed -i 's/_partial//g' NR_GENOMES_noAlpha.nex
sed -i 's/_clone//g' NR_GENOMES_noAlpha.nex
sed -i 's/_sequence//g' NR_GENOMES_noAlpha.nex
sed -i 's/_ribosomal//g' NR_GENOMES_noAlpha.nex
sed -i 's/_RNA//g' NR_GENOMES_noAlpha.nex
sed -i 's/_gene//g' NR_GENOMES_noAlpha.nex
sed -i 's/_Uncultured//g' NR_GENOMES_noAlpha.nex
sed -i 's/_bacterium_/_/g' NR_GENOMES_noAlpha.nex
sed -i 's/_alpha_proteobacterium_/_/g' NR_GENOMES_noAlpha.nex
sed -i 's/_Alphaproteobacterium//g' NR_GENOMES_noAlpha.nex
sed -i 's/_marine//g' NR_GENOMES_noAlpha.nex
sed -i 's/_marine//g' NR_GENOMES_noAlpha.nex
sed -i 's/_Rickettsiales//g' NR_GENOMES_noAlpha.nex
sed -i 's/_16S//g' NR_GENOMES_noAlpha.nex
