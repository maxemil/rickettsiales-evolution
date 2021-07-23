import glob
import pandas as pd

df = pd.read_csv('../33_metabolic_reconstruction/all_annotations_marineRick_secret4.csv', sep='\t')

tat = []
sp = []
lipo = []
for f in glob.glob("*.signalp5"):
    fdf = pd.read_csv(f, sep='\t', skiprows=1, header=0)
    fdf['# ID'] = fdf['# ID'].apply(lambda x: x.replace('&', '..'))
    fdf['spec'] = fdf['# ID'].apply(lambda x: x.split('..')[0])
    fdf['seq'] = fdf['# ID'].apply(lambda x: x.split('..')[1])
    tat.append(fdf[fdf['Prediction'] == 'TAT(Tat/SPI)'])
    sp.append(fdf[fdf['Prediction'] == 'SP(Sec/SPI)'])
    lipo.append(fdf[fdf['Prediction'] == 'LIPO(Sec/SPII)'])

signals = pd.concat(tat + sp + lipo)
tat = pd.concat(tat)
sp = pd.concat(sp)
lipo = pd.concat(lipo)

signals = df.merge(signals, how='inner', left_on='Sequence ID', right_on='seq')
tat = df.merge(tat, how='inner', left_on='Sequence ID', right_on='seq')
sp = df.merge(sp, how='inner', left_on='Sequence ID', right_on='seq')
lipo = df.merge(lipo, how='inner', left_on='Sequence ID', right_on='seq')

del tat['seq']
del tat['spec']
del tat['# ID']

signals = signals[['Sequence ID', 'Species', 'Prediction', 'protein length', 'KO', 'KEGG description',
       'KEGG pathway', 'alphaNOG', 'eggNOG annot', 'eggNOG cat',
       'PROKKA annotation', 'PROKKA gene', 'Pfam', 'Pfam description', 'CAZy',
       'CAZy description', 'TCDB', 'TCDB description', 'VFDB', 'VFDB genes',
       'VFDB description', 'VFDB keywords', 'SecReT4', 'TIGRFAM',
       'TIGRFAM description', 'IPR domain', 'IPR description', 'IPR pathway',
       'NR hits', 'NR taxonomy', 'NR taxID', 'SP(Sec/SPI)',
       'TAT(Tat/SPI)', 'LIPO(Sec/SPII)', 'OTHER', 'CS Position']]
signals.to_csv('secretion_signals.csv', sep='\t', header=True, index=False)


for gs in combinations(groups.keys(), r=2):
    in_sp = [sp for g in gs for sp in groups[g]]
    anogs_in = tat.loc[tat['spec'].isin(in_sp), 'alphaNOG'].value_counts()
    anogs_out = tat.loc[~tat['spec'].isin(in_sp), 'alphaNOG'].value_counts()
    anogs_in = set(anogs_in[anogs_in > 2].index)
    anogs_out = set(anogs_out[anogs_out > 2].index)
    print(gs)
    print(anogs_in.difference(anogs_out))

groups = {'Alpha':['Acidiphilium_cryptum_JF-5','Alpha10', 'Alpha11', 'Alpha12',
                   'Alpha32', 'Alpha35', 'Alpha96', 'Azospirillum_brasilense_Sp245',
                   'Candidatus_Puniceispirillum_marinum_IMCC1322','Caulobacter_crescentus_CB15',
                   'Rhizobium_etli_CFN', 'Rhodobacter_sphaeroides_2_4_1',
                   'Rhodospirillum_rubrum_ATCC_11170','Zymomonas_mobilis_ATCC_10988'],
          'Rickettsiaceae':['Candidatus_Arcanobacter_lacustris','Orientia_tsutsugamushi_Ikeda',
                   'PhyTra', 'Rickettsia_bellii_RML369','Rickettsia_felis_URRWXCal2',
                   'Rickettsia_prowazekii_Madrid', 'Rickettsia_rickettsii_Sheila_Smith',
                   'Rickettsiaceae_bacterium_Os18','Rickettsiales_bacterium_Ac37b'],
          'Anaplasmataceae':['Anaplasma_centrale_Israel', 'Anaplasma_marginale_Florida',
                   'Anaplasma_phagocytophilum_HZ', 'Ehrlichia_canis_Jake',
                   'Ehrlichia_chaffeensis_Arkansas', 'Ehrlichia_ruminantium_Gardel',
                   'Neorickettsia_helminthoeca_str_Oregon', 'Neorickettsia_risticii_Illinois',
                   'Neorickettsia_sennetsu_Miyayama', 'Wolbachia_endosymbiont_TRS_of_Brugia_malayi',
                   'Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel',
                   'Wolbachia_endosymbiont_of_Drosophila_melanogaster',
                   'Wolbachia_endosymbiont_of_Onchocerla_melanogaster',
                   'Candidatus_Neoehrlichia_lotoris_str_RAC413'],
          'Midichloriaceae':[ 'Candidatus_Jidaibacter_acanthamoeba',
                   'Candidatus_Midichloria_mitochondrii_IricVA',
                   'endPer', 'endosymbiont_of_Acanthamoeba_sp_UWC8'],
          'Deianiraeaceae':['DeiVas'],
          'Gamibacteraceae':['BfRick', 'UBA2645', 'bin671', 'bin673'],
          'Athabascaceae':['TARANE', 'UBA6187'],
          'Mitibacteraceae':['TOBGRS', 'UBA6149', 'UBA6178', 'UBA6189', 'bin125']}
