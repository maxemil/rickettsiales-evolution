from Bio import SeqIO
import glob
import os
from Bio.SeqFeature import FeatureLocation
import copy
from collections import defaultdict

def write_T4SS_per_spec(outfile):
    with open(outfile, 'w') as out:
        for f in glob.glob('faa_subunits/*.faa'):
            for rec in SeqIO.parse(f, 'fasta'):
                if rec.id.startswith(sp):
                    SeqIO.write(rec, out, 'fasta')

def diamond_blast(t4ss_quer, proteome, outfile):
    db = proteome.replace('faa', 'dmnd')
    cmd = "diamond makedb --in {} --db {} -p 10 --quiet".format(proteome, db)
    os.system(cmd)
    cmd = "diamond blastp -d {} -q {} --quiet -p 10 -o {} -f 6 \
            --ultra-sensitive -k 1".format(db, t4ss_quer, outfile)
    os.system(cmd)

def get_locus_tags(outfile, gbff):
    prots2seqid = defaultdict(list)
    for line in open(outfile):
        prots2seqid[line.split()[0]].append(line.split()[1])
    locus_tags = []
    for rec in SeqIO.parse(gbff, 'genbank'):
        for f in list(rec.features):
            if 'protein_id' in f.qualifiers:
                if f.qualifiers['protein_id'][0] in prots:
                    locus_tags.append(f.qualifiers['locus_tag'][0])
    assert len(locus_tags) == len(prots)
    return locus_tags

def filter_gbk_mags(infile, p, r, t4ss_seqs):
    recs = []
    for rec in SeqIO.parse(infile, 'genbank'):
        for f in list(rec.features):
            if 'locus_tag' in f.qualifiers:
                f.qualifiers['locus_tag'][0] = f.qualifiers['locus_tag'][0].replace(p, r)
                if any([f.qualifiers['locus_tag'][0].split('..')[1] in s for s in t4ss_seqs]):
                    print(f.qualifiers['locus_tag'][0])
                else:
                    rec.features.remove(f)
            else:
                rec.features.remove(f)
        if rec.features:
            recs.append(rec)
    return recs


def filter_gbk_refs(gbff, t4ss_locus_tags):
    recs = []
    for rec in SeqIO.parse(gbff, 'genbank'):
        for f in list(rec.features):
            if 'locus_tag' in f.qualifiers:
                if any([f.qualifiers['locus_tag'][0] in s for s in t4ss_locus_tags]):
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
                if f.type == 'CDS':
                    locs.append(int(f.location.end))
                    locs.append(int(f.location.start))
            if locs:
                regions = cut_locs_regions(locs)
            for rs, re in regions:
                removelist = []
                nrec = copy.deepcopy(rec)
                nrec.id = "{}_{}:{}".format(nrec.id, rs, re)
                nrec.seq = nrec.seq[rs:re]
                for f in nrec.features:
                    print(f.type, f.location)
                    if f.location.start >= rs and f.location.end <= re:
                        f.location = FeatureLocation(f.location.start-rs, 
                                             f.location.end-rs, 
                                             f.location.strand)
                    else:
                        removelist.append(f)
                for f in removelist:
                    nrec.features.remove(f)
                print("write ", nrec.id)
                SeqIO.write(nrec, out, 'genbank')

def cut_locs_regions(locs):
    regions = []
    locs = sorted(locs)
    start = locs[0]
    end = locs[0]
    for i, l in enumerate(locs):
        if not l > end + 10000:
            end = l
        else:
            regions.append((start, end))
            start = l
            end = l
    regions.append((start, end))
    return regions

def main():
    #for file target in $(cut -f2,3 /local/two/Exchange/projects/marineRick/38_T4SS_trees/references_gff_faa.csv);do  wget $target -O "$file".gz; done
    #for file target in $(cut -f4,5 /local/two/Exchange/projects/marineRick/38_T4SS_trees/references_gff_faa.csv);do  wget $target -O "$file".gz; done
    # sp = 'OriTsu'
    # proteome = "Orientia_tsutsugamushi_Ikeda.faa"
    # gbff = "Orientia_tsutsugamushi_Ikeda.gbff"
    all_locus_tags = []
    for line in open('/local/two/Exchange/projects/marineRick/38_T4SS_trees/references_gff_faa.csv'):
            line = line.split()
            sp = line[0]
            proteome = line[1]
            gbff = line[3]
            print(proteome)
            t4ss_quer = '{}_T4SS.faa'.format(sp)
            outfile = proteome.replace('faa', 'out')
            write_T4SS_per_spec(t4ss_quer)
            diamond_blast(t4ss_quer, proteome, outfile)
            locus_tags = get_locus_tags(outfile, gbff)
            recs = filter_gbk_refs(gbff, locus_tags)
            prune_gbk(recs, gbff.replace('gbff', 'gbk'))
            all_locus_tags += locus_tags


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
"bin67-3.gbk":''
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
})

t4ss_seqs = [line.strip() for line in open('T4SS_seqs.txt')]

for o, i in i2o.items():
    recs = filter_gbk_mags(i, gb2p[o], gb2r[o], t4ss_seqs)
    prune_gbk(recs, o)
