import requests
import subprocess
from Bio import SeqIO 

response = requests.get('https://www.uniprot.org/uniprot/?query=virB7+AND+rickettsiales&format=fasta')
with open('uniprot_virB7.fasta', 'w') as outfile:
    outfile.write(response.text)

# response = requests.get('https://www.uniprot.org/uniprot/?query=virB7+AND+reviewed:yes&format=fasta')
# with open('uniprot_virB7.fasta', 'a') as outfile:
#     outfile.write(response.text)

specs = ['bin_125', 'bin_67-1', 'UBA6187', 'UBA6178', 'TARA_ANE_MAG_00011',
        'Rickettsia_bellii_RML369-C', 'Deianiraea_vastatrix',
        'Candidatus_Midichloria_mitochondrii_IricVA',
        'Anaplasma_phagocytophilum_str_HZ']
for spec in specs:
    print(spec)
    cmd = "tblastn -subject ../36_genome_stats/genomes/{}.fna \
    -query bins_ref_virB7.fasta -evalue 0.5 \
    -out {}_virB7.out -outfmt '6 qseqid sseqid evalue sstart send sseq'".format(spec, spec)

    subprocess.run(cmd, shell=True)

    maxline = []
    maxlen = 0
    for line in open('{}_virB7.out'.format(spec)):
        line = line.strip().split()
        if len(line[-1]) > maxlen:
            maxline = line
            maxlen = len(line[-1])
    with open('{}_virB7.fasta'.format(spec), 'w') as outhandle:
        print('>{}_{}_{}_{}'.format(spec, maxline[1], maxline[3], maxline[4]), file=outhandle)
        print(maxline[5].replace('-', ''), file=outhandle)
        print(maxline[2])

with open('bins_ref_virB7.fasta', 'w') as outfile :
    for rec in SeqIO.parse('uniprot_virB7.fasta', 'fasta'):
        SeqIO.write(rec, outfile, 'fasta')
    for spec in specs:
        for rec in SeqIO.parse('{}_virB7.fasta'.format(spec), 'fasta'):
            SeqIO.write(rec, outfile, 'fasta')
        
cmd = "mafft-einsi --thread 40 bins_ref_virB7.fasta"
with open('bins_ref_virB7.mafft', 'w') as outfile:
    subprocess.run(cmd, shell=True, stdout=outfile)
