from Bio import SeqIO
import glob
import ete3

sub2col = {'virB1':'darksalmon',
           'virB2':'darkseagreen',
           'virB3':'greenyellow',
           'virB4':'hotpink',
           'virB5':'firebrick',
           'virB6':'lightseagreen',
           'virB8':'goldenrod',
           'virB9':'orangered',
           'virB10':'saddlebrown',
           'virB11':'seagreen',
           'virD4':'darkviolet'}

sp2ass = {'BfRick':'Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26',
'AnaPha':'Anaplasma_phagocytophilum_HZ',
'bin125':'bin_125',
'bin671':'bin_67_1',
'bin673':'bin_67_3',
'MidMit':'Candidatus_Midichloria_mitochondrii_IricVA',
'CauCre':'Caulobacter_crescentus_CB15',
'DeiVas':'Deianiraea_vastatrix',
'EhrCha':'Ehrlichia_chaffeensis_str_Arkansas',
'NeoRis':'Neorickettsia_risticii_Illinois',
'OriTsu':'Orientia_tsutsugamushi_Ikeda',
'RhiEtl':'Rhizobium_etli_CFN',
'RickBe':'Rickettsia_bellii_RML369',
'TARANE':'TARA_ANE_MAG_00011',
'TOBGRS':'TOBG_RS_372',
'UBA2645':'UBA2645',
'UBA6149':'UBA6149',
'UBA6178':'UBA6178',
'UBA6187':'UBA6187',
'UBA6189':'UBA6189',
'WolCul':'Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel'}

magfix = {'L671':'L671..L671',
          'L673':'L673..L673',
          'L125':'L125..L125',
          'TOBGRS':'TOBG_RS-372',
          'TARANE':'TARA_ANE_MAG_00011'}

prots2seqid = defaultdict(str)
for f in glob.glob("../*.out"):
    for line in open(f):
        prots2seqid[line.split()[0].split('..')[1]] = line.split()[1]
for f in glob.glob('../faa_subunits/*'):
    for rec in SeqIO.parse(f, 'fasta'):
        if 'BfRick..' in rec.id:
            fixed_id = rec.id.split('..')[1]
            prots2seqid[fixed_id] = fixed_id.split('_')[0]

def translate_prot_to_tag(gbff):
    for rec in SeqIO.parse(gbff, 'genbank'):
        for f in list(rec.features):
            if 'protein_id' in f.qualifiers:
                for s, p in prots2seqid.items():
                    if f.qualifiers['protein_id'][0] == p:
                        prots2seqid[s] = f.qualifiers['locus_tag'][0]

for f in glob.glob("genbank_files/*"):
    translate_prot_to_tag(f)


def get_correct_id(seqid):
    if seqid in prots2seqid:
        return prots2seqid[seqid]
    else:
        if not 'unlab' in seqid:
            seqid = seqid.replace('_' + seqid.split('_')[-1], '')
        else:
            seqid = seqid.replace('_unlab' + seqid.split('_unlab')[-1], '')
        if any([seqid.startswith(x) for x in magfix]):
            for p,r in magfix.items():
                seqid = seqid.replace(p, r)
            return seqid
        else:
            return seqid

with open('genoplotR_ids.tsv', 'w') as out:
    for sub, col in sub2col.items():
        for rec in SeqIO.parse("../faa_subunits/{}.faa".format(sub), 'fasta'):
            recid = rec.id.split('..')
            seqid = recid[1]
            sp = recid[0]
            if sp in sp2ass:
                seqid = get_correct_id(seqid)
                print(seqid, sp2ass[sp], col, 'plot', sep='\t', file=out)

t = ete3.PhyloTree('../../34_ALE_updated_clusters/species_trees/marineRick_b5000g20100_wDeia_clean.tree')
t.prune(sp2ass.keys())
t.ladderize(direction=1)
for sp, ass in sp2ass.items():
    n = t & sp
    n.name = ass
t.write(outfile='tree_selection_genoplotR.new', format=9)
