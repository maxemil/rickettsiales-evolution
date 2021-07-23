from Bio import SeqIO
import glob


def get_genome_size(file):
    size = 0
    for rec in SeqIO.parse(file, 'fasta'):
        size += len(rec.seq)
    return size

def get_coding_size(file):
    size = 0
    for rec in SeqIO.parse(file, 'fasta'):
        size += len(rec.seq) * 3
    return size


for line in open('proteomes_genomes.csv'):
    line = line.strip().split()
    genome_size = get_genome_size('genomes/{}'.format(line[1]))   
    coding_size = get_coding_size('proteomes/{}'.format(line[0]))   
    print(line[1].replace('.fna', ''), coding_size/genome_size)
