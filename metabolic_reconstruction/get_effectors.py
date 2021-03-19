import pandas as pd

species = ['Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26',
'Anaplasma_centrale_Israel',
'Anaplasma_marginale_Florida',
'Anaplasma_phagocytophilum_HZ',
'bin_125',
'bin_67_1',
'bin_67_3',
'Candidatus_Arcanobacter_lacustris',
'Candidatus_Fokinia_solitaria',
'Candidatus_Jidaibacter_acanthamoeba',
'Candidatus_Midichloria_mitochondrii_IricVA',
'Candidatus_Neoehrlichia_lotoris_str_RAC413',
'Candidatus_Phycoricketsia_trachydisa',
'Deianiraea_vastatrix',
'Ehrlichia_canis_Jake',
'Ehrlichia_chaffeensis_Arkansas',
'Ehrlichia_ruminantium_Gardel',
'endosymbiont_of_Acanthamoeba_sp_UWC8',
'Neorickettsia_helminthoeca_str_Oregon',
'Neorickettsia_risticii_Illinois',
'Neorickettsia_sennetsu_Miyayama',
'Orientia_tsutsugamushi_Ikeda',
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
'Wolbachia_endosymbiont_TRS_of_Brugia_malayi']

df = pd.read_csv('all_annotations_marineRick_secret4.csv', sep='\t')
apronogs = [l.strip() for l in open('effector_aproNOGs.txt')]

for a in apronogs:
    dfa = df[df['alphaNOG'] == a]
    if set(dfa['Species']).issubset(set(species)):
        print(a, ";".join(set(dfa['PROKKA gene'].dropna())), ";".join(set(dfa['KO'].dropna())), ";".join(set(dfa['SecReT4'].dropna())), sep='\t')
