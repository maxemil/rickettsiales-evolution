import pandas as pd

df = pd.read_csv('all_annotations_marineRick_secret4.csv', sep='\t')

effectors = set()
for ss in df[-df['SecReT4'].isna()]['SecReT4']:
    for s in ss.split(';'):
        effectors.add(s.strip())

ef = pd.DataFrame(index=effectors, columns=set(df['Species'].dropna()))
ef = ef.fillna(0)


for row in df[-df['SecReT4'].isna()].iterrows():
    for s in row[1]['SecReT4'].split(';'):
        s = s.strip()
        ef.loc[s, row[1]['Species']] += 1
        
ef = ef.loc[ef.sum(axis=1).sort_values(ascending=False).index]
ef = ef[[sp for k,g in groups.items() for sp in g]]
ef.to_csv('effectors_SecReT4.csv', sep='\t', header=True, index=True)





groups = {'Rickettsiaceae':['Candidatus_Arcanobacter_lacustris','Orientia_tsutsugamushi_Ikeda',
                   'Candidatus_Phycoricketsia_trachydisa', 'Rickettsia_bellii_RML369','Rickettsia_felis_URRWXCal2',
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
                   'Rickettsiales_endosymbiont_of_Peranema_trichophorum', 'endosymbiont_of_Acanthamoeba_sp_UWC8'],
          'Deianiraeaceae':['Deianiraea_vastatrix'],
          'Gamibacteraceae':['Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26', 'UBA2645', 'bin_67_1', 'bin_67_3'],
          'Athabascaceae':['TARA_ANE_MAG_00011', 'UBA6187'],
          'Mitibacteraceae':['TOBG_RS-372', 'UBA6149', 'UBA6178', 'UBA6189', 'bin_125'],
          'Alpha':['Acidiphilium_cryptum_JF-5','Alphaproteobacteria_bacterium_MarineAlpha10_Bin2', 
                             'Alphaproteobacteria_bacterium_MarineAlpha11_Bin1', 'Alphaproteobacteria_bacterium_MarineAlpha12_Bin1',
                             'Alphaproteobacteria_bacterium_MarineAlpha3_Bin2', 'Alphaproteobacteria_bacterium_MarineAlpha3_Bin5', 'Alphaproteobacteria_bacterium_MarineAlpha9_Bin6', 'Azospirillum_brasilense_Sp245',
                             'Candidatus_Puniceispirillum_marinum_IMCC1322','Caulobacter_crescentus_CB15',
                             'Rhizobium_etli_CFN', 'Rhodobacter_sphaeroides_2_4_1',
                             'Rhodospirillum_rubrum_ATCC_11170','Zymomonas_mobilis_ATCC_10988']}
