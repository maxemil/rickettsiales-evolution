import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

# for f in glob.glob('../cluster_only/*.faa'):
#     for rec in SeqIO.parse(f, 'fasta'):
#         T4SS_seqs.append(rec.id)


def filter_gbk(infile, p, r):
    recs = []
    for rec in SeqIO.parse(infile, 'genbank'):
        for f in list(rec.features):
            if 'locus_tag' in f.qualifiers:
                f.qualifiers['locus_tag'][0] = f.qualifiers['locus_tag'][0].replace(p, r)
                if any([f.qualifiers['locus_tag'][0].split('..')[1] in s for s in T4SS_seqs]):
                    print(f.qualifiers['locus_tag'][0])
                else:
                    rec.features.remove(f)
            else:
                rec.features.remove(f)
        if rec.features:
            recs.append(rec)
    return recs

def prune_gbk(recs, outfile):
    with open(outfile, 'w') as out:
        for rec in recs:
            locs = []
            for f in rec.features:
                locs.append(int(f.location.end))
                locs.append(int(f.location.start))
            offset = min(locs)
            for f in rec.features:
                f.location = FeatureLocation(f.location.start-offset, 
                                             f.location.end-offset, 
                                             f.location.strand)
            rec.seq = rec.seq[offset:max(locs)]
            SeqIO.write(rec, out, 'genbank')
            
            
i2o = {
"TOBG_RS-372.gbk":'TOBG_RS-372.gb',
"UBA6149.gbk":'UBA6149.gb',
"UBA6178.gbk":'UBA6178.gb',
"UBA6189.gbk":'UBA6189.gb',
"UBA6187.gbk":'UBA6187.gb',
"UBA2645.gbk":'UBA2645.gb',
"TARA_ANE_MAG_00011.gbk":'TARA_ANE_MAG_00011.gb',
"bin125.gbk":'../../33_metabolic_reconstruction/prokka/bin_125/PROKKA_06182020.gbk',
"bin67-1.gbk":'../../33_metabolic_reconstruction/prokka/bin_67_1/PROKKA_06182020.gbk',
"bin67-3.gbk":'../../33_metabolic_reconstruction/prokka/bin_67_3/PROKKA_06182020.gbk'
}

gb2p = {
'UBA6149.gbk':'UBA6149',
'UBA6178.gbk':'UBA6178',
'UBA6189.gbk':'UBA6189',
'UBA6187.gbk':'UBA6187',
'UBA2645.gbk':'UBA2645',
'TARA_ANE_MAG_00011.gbk':'TARA_ANE_MAG_00011',
'TOBG_RS-372.gbk':'TOBG_RS-372',
"bin125.gbk":'',
"bin67-1.gbk":'',
"bin67-3.gbk":'',
}

gb2r = defaultdict(str,{
'UBA6149.gbk':'UBA6149..UBA6149',
'UBA6178.gbk':'UBA6178..UBA6178',
'UBA6189.gbk':'UBA6189..UBA6189',
'UBA6187.gbk':'UBA6187..UBA6187',
'UBA2645.gbk':'UBA2645..UBA2645',
'TARA_ANE_MAG_00011.gbk':'TARANE..TARANE',
'TOBG_RS-372.gbk':'TOBGRS..TOBGRS',
"bin125.gbk":'',
"bin67-1.gbk":'',
"bin67-3.gbk":''
}

t4ss_seqs = [line.strip() for line in open('T4SS_seqs.txt')]

for o, i in i2o.items():
    recs = filter_gbk(i, gb2p[o], gb2r[o])
    prune_gbk(recs, o)
