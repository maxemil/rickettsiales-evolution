import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('paraller_sim_query_2274/IMNGS_99/IMNGS_99_counts_matrix.tab', 
                    sep='\t', index_col=0)

df = df[df['1'] > 10]

df.loc[df['Description'] == 'seawater metagenome', 'Description'] = 'marine metagenome'
df.loc[df['Description'] == '16S rRNA gene amplicons', 'Description'] = 'marine metagenome'
df.loc[df['Description'] == 'coral metagenome', 'Description'] = 'filter feeder metagenome'
df.loc[df['Description'] == 'Fungia', 'Description'] = 'filter feeder metagenome'
df.loc[df['Description'] == 'gut metagenome', 'Description'] = 'filter feeder metagenome'

hist = df['Description'].value_counts()

sns.barplot(y=hist, x=hist.index)
