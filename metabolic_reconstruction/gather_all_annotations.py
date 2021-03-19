from collections import defaultdict
from Bio.KEGG import REST
from urllib.error import HTTPError
from Bio import SeqIO
import ete3

header = ['Species', 'protein length', 'KO', 'KEGG description', \
          'KEGG pathway', 'alphaNOG', 'eggNOG annot', 'eggNOG cat', \
          'PROKKA annotation', 'PROKKA gene', 'Pfam', 'Pfam description', \
          'CAZy', 'CAZy description', 'TCDB', 'TCDB description', 'VFDB', \
          'VFDB genes', 'VFDB description', 'VFDB keywords', 'SecReT4', 'TIGRFAM', \
          'TIGRFAM description', 'IPR domain', 'IPR description', \
          'IPR pathway', 'NR hits', 'NR taxonomy', 'NR taxID']

def def_annot():
    return {'Species':'',
           'CAZy':[],
           'CAZy description':set(),
           'NR hits':[],
           'NR taxonomy':"",
           'NR taxID':"",
           'KO':[],
           'KEGG description':[],
           'KEGG pathway':set(),
           'alphaNOG':'',
           'eggNOG annot':'',
           'eggNOG cat':'',
           'TCDB':[],
           'TCDB description':set(),
           'VFDB':set(),
           'VFDB genes':set(),
           'VFDB description':set(),
           'VFDB keywords':set(),
           'SecReT4':"",
           'protein length':'',
           'Pfam':[],
           'Pfam description':[],
           'TIGRFAM':[],
           'TIGRFAM description':[],
           'Gene3D':[],
           'Gene3D description':[],
           'SUPERFAMILY':[],
           'SUPERFAMILY description':[],
           'Hamap':[],
           'Hamap description':[],
           'PIRSF':[],
           'PIRSF description':[],
           'SMART':[],
           'SMART description':[],
           'CDD':[],
           'CDD description':[],
           'Coils':[],
           'Coils description':[],
           'SFLD':[],
           'SFLD description':[],
           'IPR domain':set(),
           'IPR description':set(),
           'IPR pathway':set(),
           'PROKKA annotation':'',
           'PROKKA gene':''}

species = ['Acidiphilium_cryptum_JF-5',
'Alphaproteobacteria_bacterium_MarineAlpha10_Bin2',
'Alphaproteobacteria_bacterium_MarineAlpha11_Bin1',
'Alphaproteobacteria_bacterium_MarineAlpha12_Bin1',
'Alphaproteobacteria_bacterium_MarineAlpha3_Bin2',
'Alphaproteobacteria_bacterium_MarineAlpha3_Bin5',
'Alphaproteobacteria_bacterium_MarineAlpha9_Bin6',
'Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26',
'Anaplasma_centrale_Israel',
'Anaplasma_marginale_Florida',
'Anaplasma_phagocytophilum_HZ',
'Azospirillum_brasilense_Sp245',
'bin_125',
'bin_67_1',
'bin_67_3',
'Candidatus_Arcanobacter_lacustris',
'Candidatus_Fokinia_solitaria',
'Candidatus_Jidaibacter_acanthamoeba',
'Candidatus_Midichloria_mitochondrii_IricVA',
'Candidatus_Neoehrlichia_lotoris_str_RAC413',
'Candidatus_Phycoricketsia_trachydisa',
'Candidatus_Puniceispirillum_marinum_IMCC1322',
'Candidatus_Xenolissoclinum_pacificiensis',
'Caulobacter_crescentus_CB15',
'Deianiraea_vastatrix',
'Ehrlichia_canis_Jake',
'Ehrlichia_chaffeensis_Arkansas',
'Ehrlichia_ruminantium_Gardel',
'endosymbiont_of_Acanthamoeba_sp_UWC8',
'Neorickettsia_helminthoeca_str_Oregon',
'Neorickettsia_risticii_Illinois',
'Neorickettsia_sennetsu_Miyayama',
'Orientia_tsutsugamushi_Ikeda',
'Rhizobium_etli_CFN',
'Rhodobacter_sphaeroides_2_4_1',
'Rhodospirillum_rubrum_ATCC_11170',
'Rickettsia_bellii_RML369',
'Rickettsiaceae_bacterium_Os18',
'Rickettsia_felis_URRWXCal2',
'Rickettsiales_bacterium_Ac37b',
'Rickettsiales_endosymbiont_of_Peranema_trichophorum',
'Rickettsiales_endosymbiont_of_Stachyamoeba_lipophora',
'Rickettsia_prowazekii_Madrid',
'Rickettsia_rickettsii_Sheila_Smith',
'TARA_ANE_MAG_00011',
'TOBG_RS-372',
'UBA2645',
'UBA6149',
'UBA6178',
'UBA6187',
'UBA6189',
'Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel',
'Wolbachia_endosymbiont_of_Drosophila_melanogaster',
'Wolbachia_endosymbiont_of_Onchocerla_melanogaster',
'Wolbachia_endosymbiont_TRS_of_Brugia_malayi',
'Zymomonas_mobilis_ATCC_10988']


def parse_eggnog():
    for sp in species:
        for line in open('EggNOG/{}.emapper.annotations'.format(sp)):
            if not line.startswith('#'):
                line = line.strip().split('\t')
                seq = line[0].split('..')[1]
                annotations[seq]['KO'] = line[6].split(',') if line[6] else []
                annotations[seq]['alphaNOG'] = line[10].split('|')[0]
                annotations[seq]['eggNOG cat'] = line[11]
                annotations[seq]['eggNOG annot'] = line[12]
                OGs = {o.split('@')[1]:o.split('@')[0] for o in line[9].split(',')}

def parse_KEGG():
    ko_desc = {}
    ko_path = {}
    for line in open('KEGG_KO_Pathways.csv'):
        line = line.split('\t')
        ko_desc[line[0]] = line[1]
        ko_path[line[0]] = line[2].strip().split(';')

    for k,v in annotations.items():
        if annotations[k]['KO']:
            for ko in annotations[k]['KO']:
                annotations[k]['KEGG description'].append(ko_desc[ko])
                for p in ko_path[ko]:
                    annotations[k]['KEGG pathway'].add(p)

def parse_CAZy():
    cazy_annot = defaultdict(str)
    for line in open("CAZy/fam_activities.txt"):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            if len(line) > 1:
                cazy_annot[line[0]] = line[1].strip()

    cazy_seq = defaultdict(lambda: {sp:set() for sp in species})
    for sp in species:
        for line in open("CAZy/{}.out".format(sp)):
            if not line.startswith('#'):
                line = line.split()
                seq = line[0].split('..')[1]
                cazy = line[2].replace('.hmm', '')
                annotations[seq]['CAZy'].append(cazy)
                annotations[seq]['CAZy description'].add(cazy_annot[cazy.split('_')[0]])

def parse_TCDB():
    tcdb_annot = {}
    for line in open("TCDB/families.tsv"):
        line = line.strip().split('\t')
        tcdb_annot[line[0]] = line[1].replace('</sup>', '').replace('<sup>', '')

    for sp in species:
        for line in open("TCDB/{}.out".format(sp)):
            line = line.split()
            seq = line[0].split('..')[1]
            if len(line[1].split('|')) == 4:
                fam = line[1].split('|')[3]
                superfam = ".".join(fam.split('.')[0:3])
                annotations[seq]['TCDB'].append(fam)
                annotations[seq]['TCDB description'].add(tcdb_annot[superfam])

def parse_VFDB():
    VFs = defaultdict(tuple)
    for line in open('VFDB/VFs.csv'):
        line = line.strip().split('\t')
        if line[1]:
            VFs[line[8]] = ("{} ({})".format(line[0], line[1]), line[3], line[5], line[7])
        else:
            VFs[line[8]] = (line[0], line[3], line[5], line[7])
    seq2VF = {}
    seq2desc = {}
    for rec in SeqIO.parse('VFDB/VFDB_setB_pro.fas', 'fasta'):
        for elem in rec.description.split('('):
            if elem.split(')')[0] in VFs:
                seq2VF[rec.id] = elem.split(')')[0]
        if not rec.id in seq2VF:
                seq2desc[rec.id] = rec.description.replace(rec.id, '')
    for sp in species:
        for line in open("VFDB/{}.out".format(sp)):
            if not line.startswith('#'):
                line = line.strip().split('\t')
                seq = line[0].split('..')[1]
                if line[1] in seq2VF:
                    VF = seq2VF[line[1]]
                    annotations[seq]['VFDB'].add(VF)
                    annotations[seq]['VFDB genes'].add(VFs[VF][0])
                    annotations[seq]['VFDB description'].add(VFs[VF][1])
                    annotations[seq]['VFDB description'].add(VFs[VF][2])
                    for kw in VFs[VF][3].split(';'):
                        annotations[seq]['VFDB keywords'].add(kw)
                else:
                    annotations[seq]['VFDB description'].add(seq2desc[line[1]])
                    
def parse_interproscan():
    for sp in species:
        for line in open("INTERPROSCAN/{}.faa.tsv".format(sp)):
            line = line.strip().split('\t')
            seq = line[0].split('..')[1]
            analysis = line[3]
            annotations[seq][analysis].append(line[4])
            annotations[seq]["{} description".format(analysis)].append(line[5])
            for i, k in zip([11, 12, 13], ['IPR domain', 'IPR description', 'IPR pathway']):
                if len(line) >= i + 1:
                    annotations[seq][k].add(line[i])

def parse_prokka():
    for sp in species:
        for line in open("prokka/{}.gff".format(sp)):
            line = line.strip().split('\t')
            if len(line) == 9 and not line[2] == 'gene':
                line_dict = {}
                for elem in line[8].split(';'):
                    elem = elem.split('=')
                    line_dict[elem[0]] = elem[1]
                try:
                    annotations[line_dict['ID'].split('..')[1]]['PROKKA annotation'] = line_dict['product']
                    annotations[line_dict['ID'].split('..')[1]]["Species"] = sp
                except:
                    print(line)
                if 'gene' in line_dict:
                    annotations[line_dict['ID'].split('..')[1]]['PROKKA gene'] = line_dict['gene']

def ncbi_get_common_ancestor(taxids, ncbi):
    taxids = [t for t in taxids if t]
    if not taxids:
        return '', ''
    t = ncbi.get_topology(taxids, intermediate_nodes=False)
    tax = ncbi.get_taxid_translator([int(t.name)])[int(t.name)]
    return tax, t.name

def parse_diamond():
    ncbi = ete3.ncbi_taxonomy.NCBITaxa()
    seq2taxids = defaultdict(set)
    for sp in species:
        for line in open('DIAMOND_NR/{}.faa.out'.format(sp)): 
            line = line.strip().split('\t')
            sequence = line[0].split('..')[1]
            annotations[sequence]["NR hits"].append("{}:{}".format(line[1], line[3]))
            for t in line[2].split(';'):
                seq2taxids[sequence].add(t)
    for seq, taxids in seq2taxids.items():
        tax, taxid = ncbi_get_common_ancestor(taxids, ncbi)
        annotations[seq]['NR taxonomy'] = tax
        annotations[seq]['NR taxID'] = taxid

def parse_secret4():
    for sp in species:
        for line in open('SecReT4/{}.out'.format(sp)):
            line = line.strip().split('\t')
            seq = line[0].split('..')[1]
            hit = line[12].split('|')[4]
            if not annotations[seq]['SecReT4']:
                annotations[seq]['SecReT4'] = hit
            
def parse_sequences():
    for sp in species:
        for rec in SeqIO.parse("proteomes/{}.faa".format(sp), 'fasta'):
            seq = rec.id.split('..')[1]
            annotations[seq]["protein length"] = str(len(rec.seq))
            annotations[seq]["Species"] = sp

def print_annotations():
    with open("all_annotations_marineRick.csv", 'w') as out:
        print('\t'.join(['Sequence ID'] + header), file=out)
        for k, v in annotations.items():
            cols = [k]
            for h in header:
                if isinstance(v[h], list) or isinstance(v[h], set):
                    cols.append(';'.join(v[h]))
                else:
                    cols.append(v[h])
            print('\t'.join(cols), file=out)

def main():
    annotations = defaultdict(lambda: def_annot())
    parse_eggnog()
    parse_KEGG()
    parse_CAZy()
    parse_TCDB()
    parse_VFDB()
    parse_prokka()
    parse_interproscan()
    parse_secret4()
    parse_sequences()
    # parse_diamond()
    print_annotations()

if __name__ == '__main__':
    main()

# def prepare_KEGG_lookup():
#     KOs = set()
#     for sp in species:
#         for line in open('{}.emapper.annotations'.format(sp)):
#             if not line.startswith('#'):
#                 line = line.strip().split('\t')
#                 if line[6]:
#                     for ko in line[6].split(','):
#                         KOs.add(ko)
#
    general_pathways = ['path:map01110', 'path:map01100', 'path:map01230', 'path:map01120']
    KO_path = defaultdict(lambda: {'description':"", 'pathways':set()})
    for ko in KOs:
        try:
            KO_path[ko]['description'] = REST.kegg_list(ko).read().strip().split('\t')[1]
            paths = REST.kegg_link('pathway', ko).read().strip()
            if paths:
                for path in paths.split('\n'):
                    path = path.split('\t')
                    if path[1].startswith('path:map') and not any([path[1] == p for p in general_pathways]):
                        path_str = ':'.join(REST.kegg_list(path[1]).read().strip().split('\t'))
                        KO_path[ko]['pathways'].add(path_str)
        except:
            print("{} did not finish, errors occured".format(ko))
            KO_path[ko]

#     with open('KEGG_KO_Pathways.csv', 'a') as out:
#         for k,v in KO_path.items():
#             print("{}\t{}\t{}".format(k, v['description'], ';'.join(v['pathways'])), file=out)
