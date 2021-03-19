import pandas as pd

df = pd.read_csv('all_annotations_marineRick_secret4.csv', sep='\t')
df = df[df['Species'].isin(['bin_125', 'TOBG_RS-372', 'UBA6149', 'UBA6178', 
                    'UBA6189', 'TARA_ANE_MAG_00011', 'UBA6187',
                    'Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26', 
                    'bin_67-1', 'bin_67-3', 'UBA2645'])]

pdf = pd.DataFrame(index=set(df['alphaNOG'].dropna()), columns=set(df['Species'].dropna()))

pdf_dummy = pdf.loc[['01S0K', '01V6A']]
pdf_dummy.stack(dropna=False).to_frame().apply(lambda x: len(df.loc[(df['alphaNOG'] == x.name[0]) & (df['Species'] == x.name[1])]), axis=1).unstack()

pdf = pdf.stack(dropna=False).to_frame().apply(lambda x: len(df.loc[(df['alphaNOG'] == x.name[0]) & (df['Species'] == x.name[1])]), axis=1).unstack()

pdf['eggNOG annot'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'eggNOG annot'].iloc[0], axis=1)
pdf['eggNOG cat'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'eggNOG cat'].iloc[0], axis=1)
pdf['KEGG description'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'KEGG description'].iloc[0], axis=1)
pdf['KEGG pathway'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'KEGG pathway'].iloc[0], axis=1)
pdf['KO'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'KO'].iloc[0], axis=1)

pdf.to_csv('STx_clusters_mitibacter_UBA6187_Paradeianiraeaceae.csv', header=True, sep='\t', index=True)


df = pd.read_csv('all_annotations_marineRick_secret4.csv', sep='\t')
df = df[df['Species'].isin(['bin_125', 'TOBG_RS-372', 'UBA6149', 'UBA6178', 
                    'UBA6189', 'TARA_ANE_MAG_00011', 'UBA6187',
                    'Alphaproteobacteria_bacterium_RIFCSPLOWO2_01_FULL_40_26', 
                    'bin_67-1', 'bin_67-3', 'UBA2645'])]
df = df[~df['CAZy'].isna()]

pdf = pd.DataFrame(index=set(df['alphaNOG'].dropna()), columns=set(df['Species'].dropna()))
pdf = pdf.stack(dropna=False).to_frame().apply(lambda x: len(df.loc[(df['alphaNOG'] == x.name[0]) & (df['Species'] == x.name[1])]), axis=1).unstack()
pdf['CAZy'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'CAZy'].iloc[0], axis=1)
pdf['CAZy description'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'CAZy description'].iloc[0], axis=1)
pdf.to_csv('clst_mitibacter_UBA6187_Paradeianiraeaceae_CAZy.csv', header=True, sep='\t', index=True)



#############################################
df = pd.read_csv('all_annotations_marineRick_secret4.csv', sep='\t')

alphanogs = []
for line in open('STx_highlight_pathways.csv'):
    line = line.strip().split()
    if line[0].startswith('0'):
        alphanogs.append(line[0])

df = df[df['alphaNOG'].isin(alphanogs)]

pdf = pd.DataFrame(index=set(df['alphaNOG'].dropna()), columns=set(df['Species'].dropna()))
pdf = pdf.stack(dropna=False).to_frame().apply(lambda x: len(df.loc[(df['alphaNOG'] == x.name[0]) & (df['Species'] == x.name[1])]), axis=1).unstack()

pdf['eggNOG annot'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'eggNOG annot'].iloc[0], axis=1)
pdf['eggNOG cat'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'eggNOG cat'].iloc[0], axis=1)
pdf['KEGG description'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'KEGG description'].iloc[0], axis=1)
pdf['KEGG pathway'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'KEGG pathway'].iloc[0], axis=1)
pdf['KO'] = pdf.apply(lambda x: df.loc[df['alphaNOG'] == x.name, 'KO'].iloc[0], axis=1)

with open('STx_highlight_pathways_annot.csv', 'w') as out:
    print("\t".join(['cluster', 'gene'] + list(pdf.columns)), file=out)
    for line in open('STx_highlight_pathways.csv'):
        if line.startswith('0'):
            line = line.strip('\n').split()
            if len(line) == 1:
                line += ['']
            print("\t".join(line + list(pdf.astype(str).loc[line[0]])), file=out)
        else:
            print(line.strip('\n'), file=out)
