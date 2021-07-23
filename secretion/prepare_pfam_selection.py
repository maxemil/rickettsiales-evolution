from collections import defaultdict
import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

order = [
'Anaplasmataceae',
'Anaplasma_marginale_Florida',
'Anaplasma_centrale_Israel',
'Anaplasma_phagocytophilum_HZ',
'Candidatus_Neoehrlichia_lotoris_str_RAC413',
'Ehrlichia_chaffeensis_Arkansas',
'Ehrlichia_canis_Jake', 
'Ehrlichia_ruminantium_Gardel',
'Wolbachia_endosymbiont_of_Drosophila_melanogaster',
'Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel',
'Wolbachia_endosymbiont_of_Onchocerla_melanogaster',
'Wolbachia_endosymbiont_TRS_of_Brugia_malayi',
'Neorickettsia_sennetsu_Miyayama',
'Neorickettsia_risticii_Illinois',
'Neorickettsia_helminthoeca_str_Oregon',
'Midichloriaceae',
'endosymbiont_of_Acanthamoeba_sp_UWC8',
'Candidatus_Jidaibacter_acanthamoeba',
'Rickettsiales_endosymbiont_of_Peranema_trichophorum', 
'Candidatus_Midichloria_mitochondrii_IricVA',
'Gamibaceraceae',
'UBA2645',
'Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26', 
'bin_67_1',
'bin_67_3', 
'Deianiraeaceae',
'Deianiraea_vastatrix',
'Rickettsiaceae', 
'Rickettsia_rickettsii_Sheila_Smith',
'Rickettsia_felis_URRWXCal2', 
'Rickettsia_prowazekii_Madrid',
'Rickettsia_bellii_RML369',
'Rickettsiaceae_bacterium_Os18',
'Orientia_tsutsugamushi_Ikeda', 
'Candidatus_Phycoricketsia_trachydisa',
'Rickettsiales_bacterium_Ac37b', 
'Candidatus_Arcanobacter_lacustris',
'Athabascaceae',
'UBA6187',
'TARA_ANE_MAG_00011',
'Mitibacteraceae',
'bin_125',
'UBA6178',
'TOBG_RS-372',
'UBA6189',
'UBA6149',
'Outgroup',
'Alphaproteobacteria_bacterium_MarineAlpha3_Bin5',
'Alphaproteobacteria_bacterium_MarineAlpha3_Bin2',
'Alphaproteobacteria_bacterium_MarineAlpha12_Bin1',
'Alphaproteobacteria_bacterium_MarineAlpha11_Bin1',
'Alphaproteobacteria_bacterium_MarineAlpha9_Bin6',
'Alphaproteobacteria_bacterium_MarineAlpha10_Bin2',
'Candidatus_Puniceispirillum_marinum_IMCC1322',
'Azospirillum_brasilense_Sp245',
'Acidiphilium_cryptum_JF-5',
'Rhodospirillum_rubrum_ATCC_11170',
'Rhodobacter_sphaeroides_2_4_1',
'Rhizobium_etli_CFN',
'Caulobacter_crescentus_CB15',
'Zymomonas_mobilis_ATCC_10988']

clans_ids = []
clan_labels = {}
for line in open('effector_clans.txt'):
    line = line.split()
    clans_ids.append(line[0])
    clan_labels[line[0]] = "{} ({})".format(line[0], line[1])

profiles = {}
for line in open('Pfam-A.clans.tsv'):
    if any([c in line for c in clans_ids]):
        line = line.split()
        if line[1] in clans_ids:
            profiles[line[0]] = line[1]
        else:
            profiles[line[0]] = line[0]

with open('Pfam_effector_selection.txt', 'w') as out:
    for p in profiles.keys():
        print(p, file=out)

# for d in $(cat Pfam_effector_selection.txt); do wget http://pfam.xfam.org/family/$d/hmm -O $d.hmm; cat $d.hmm >> Pfam_effector_selection.hmm; rm $d.hmm; done
# for f in ../33_metabolic_reconstruction/proteomes/*.faa; do hmmsearch --tblout outfiles/$(basename ${f%%.faa}).out --cpu 5 -E 1e-05 Pfam_effector_selection.hmm $f;done

specs = [os.path.basename(f).replace('.out', '') for f in glob.glob('outfiles/*.out')]

df = pd.DataFrame(0, columns=specs, index=clan_labels.values())

for f in glob.glob('outfiles/*.out'):
    spec = os.path.basename(f).replace('.out', '')
    seqs = defaultdict(set)
    for line in open(f):
        #print(line)
        if line.startswith('#'):
            continue
        line = line.split()
        seqs[profiles[line[3].split('.')[0]]].add(line[0])
    for k, v in seqs.items():
        df.loc[clan_labels[k], spec] = len(v)

df[['Anaplasmataceae', 'Rickettsiaceae', 'Outgroup', 'Midichloriaceae', 
    'Athabascaceae', 'Deianiraeaceae', 'Gamibaceraceae', 'Mitibacteraceae']] = 0
df = df[order]
fig, ax = plt.subplots(figsize=(25, 8))
sns.heatmap(df, annot=True, 
        cmap=sns.color_palette("viridis_r", as_cmap=True), 
        ax=ax, fmt='d', norm=LogNorm())
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=45)
plt.tight_layout()
plt.savefig("effector_domains_heatmap.pdf")
plt.cla()
