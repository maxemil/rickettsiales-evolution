import sys
import os

from Bio import SeqIO
from BCBio import GFF
from Bio.Alphabet import generic_dna

def main(gff_file, fasta_file):
    out_file = "%s.gb" % os.path.basename(gff_file).split('.')[0]
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, base_dict=fasta_input)
    with open(out_file, 'w') as out:
        for rec in gff_iter:
            rec.seq.alphabet = generic_dna 
            SeqIO.write(rec, out, "genbank")

if __name__ == "__main__":
    main(*sys.argv[1:])
