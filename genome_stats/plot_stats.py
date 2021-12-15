import pandas as pd
import seaborn as sns
import seaborn as sns
import matplotlib.pyplot as plt
import colorsys
import matplotlib.colors

clade2color = {'Rickettsiaceae':'goldenrod', 'Midichloriaceae':'olive', 
                  'Anaplasmataceae':'chocolate', 'Deianiraeaceae':'firebrick', 
                  'Alphaproteobacteria':'grey', 'Mitibacteraceae':'cyan', 
                  'UBA6187-group':'teal'}
order = ['Anaplasma_marginale_str_Florida',
'Anaplasma_centrale_str_Israel',
'Anaplasma_phagocytophilum_str_HZ',
'Candidatus_Neoehrlichia_lotoris_str_RAC413',
'Ehrlichia_chaffeensis_str_Arkansas',
'Ehrlichia_canis_str_Jake', 
'Ehrlichia_ruminantium_str_Gardel',
'Wolbachia_endosymbiont_of_Drosophila_melanogaster',
'Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel',
'Wolbachia_endosymbiont_of_Onchocerca_ochengi',
'Wolbachia_endosymbiont_strain_TRS_of_Brugia_malayi',
'Neorickettsia_sennetsu_str_Miyayama',
'Neorickettsia_risticii_str_Illinois',
'Neorickettsia_helminthoeca_str_Oregon',
'Endosymbiont_of_Acanthamoeba_sp_UWC8',
'Candidatus_Jidaibacter_acanthamoeba_strain_UWC36_NF27_EY',
'Rickettsiales_endosymbiont_of_Peranema_trichophorum', 
'Candidatus_Midichloria_mitochondrii_IricVA',
'UBA2645',
'Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26', 
'bin_67-1',
'bin_67-3', 
'Deianiraea_vastatrix', 
'Rickettsia_rickettsii_str_Sheila_Smith',
'Rickettsia_felis_URRWXCal2', 
'Rickettsia_prowazekii_str_Madrid_E',
'Rickettsia_bellii_RML369-C',
'Occidentia_massiliensis_strain_Os18',
'Orientia_tsutsugamushi_str_Ikeda', 
'Candidatus_Phycoricketsia_trachydisa',
'Rickettsiales_bacterium_Ac37b', 
'Candidatus_Arcanobacter_lacustris_strain_SCGC',
'UBA6187',
'TARA_ANE_MAG_00011',
'bin_125',
'UBA6178',
'TOBG_RS-372',
'UBA6189',
'UBA6149',
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


def scale_lightness(rgb, scale_l):
    h, l, s = colorsys.rgb_to_hls(*rgb)
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

df = pd.read_csv('all_genomes.csv', sep='\t', header=0, index_col=0)
df['Estimated length'] = df['Length'] * 1/df['Completeness']
df['Color'] = df['Clade'].apply(lambda x: matplotlib.colors.to_rgb(clade2color[x]))
df = df.reindex(order)
df['coding density'] = df['coding density'] * 100


f, ax = plt.subplots(1,3, figsize=(20, 10), sharey=True)

sns.barplot(x='Estimated length', y=df.index, data=df, 
            palette=sns.color_palette(df['Color'].apply(lambda x: scale_lightness(x, 1.2)), as_cmap=True),
            ax=ax[0])
sns.barplot(x='Length', y=df.index, data=df, 
            palette=sns.color_palette(df['Color'].apply(lambda x: scale_lightness(x, 0.7)), as_cmap=True), 
            ax=ax[0])

ax[0].yaxis.grid(False)
ax[0].xaxis.grid(True)

sns.barplot(x='GC-content', y=df.index, data=df, color='grey', ax=ax[1],
            palette=sns.color_palette(df['Color'], as_cmap=True))
ax[1].yaxis.grid(False)
ax[1].xaxis.grid(True)

sns.barplot(x='coding density', y=df.index, data=df, color='grey', ax=ax[2],
            palette=sns.color_palette(df['Color'], as_cmap=True))
ax[2].yaxis.grid(False)
ax[2].xaxis.grid(True)

# plt.show()
plt.tight_layout()
plt.savefig('length_GC.pdf')
plt.close()

f, ax = plt.subplots(1,2, figsize=(15, 4), sharey=True)

sns.boxplot(x='Estimated length', y='Clade', data=df, 
            palette=sns.color_palette([df.loc[df['Clade'] == c, 'Color'].iloc[0] for c in sorted(set(df['Clade']))], as_cmap=True), ax=ax[0])

sns.boxplot(x='GC-content', y='Clade', data=df, 
            palette=sns.color_palette([df.loc[df['Clade'] == c, 'Color'].iloc[0] for c in sorted(set(df['Clade']))], as_cmap=True), ax=ax[1])
ax[0].yaxis.grid(False)
ax[0].xaxis.grid(True)
ax[1].yaxis.grid(False)
ax[1].xaxis.grid(True)
plt.tight_layout()
plt.savefig('length_GC_box.pdf')


from scipy.stats import ttest_ind
cat1 = df[df['Clade'].apply(lambda x: x in ['Mitibacteraceae', 'UBA6187-group'])]
cat2 = df[df['Clade'].apply(lambda x: x in ['Rickettsiaceae','Anaplasmataceae','Midichloriaceae','Deianiraeaceae','Gamibacteraceae'])]
cat3 = df[df['Clade'].apply(lambda x: x in ['Alphaproteobacteria'])]

for val in ['GC-content', 'Estimated length', 'coding density']:
    print(ttest_ind(cat1[val], cat2[val]))
    print(ttest_ind(cat3[val], cat2[val]))

    print(cat1[val].describe())
    print(cat2[val].describe())
    print(cat3[val].describe())
