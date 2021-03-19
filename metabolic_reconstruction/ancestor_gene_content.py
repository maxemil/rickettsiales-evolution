import pandas as pd
from math import modf


header = ['species_tree', 'cluster', 'node', 'duplications',
          'transfers', 'losses', 'originations', 'copies']
events = pd.read_csv('../34_ALE_updated_clusters/events_wDeia.txt', 
                            sep='\t', names=header)
del events['species_tree']
events['cluster'] = events['cluster'].str.replace('.clean.ufboot.ale', '')

for e in ['duplications', 'transfers', 'losses', 'originations', 'copies']:
    event = pd.pivot_table(events, values=e, index=['cluster'], columns=['node'])
    event.loc[event.sum(axis=1).sort_values(ascending=False).index]
    event.fillna(0, inplace=True)
    event.to_csv('{}_all.csv'.format(e), sep='\t', header=True, index=True)


df = pd.read_csv('all_annotations_marineRick_secret4.csv', sep='\t')
def get_annot(clst, col):
    try:
        return df.loc[df['alphaNOG'] == clst, col].iloc[0]
    except:
        return ''
        
nodes = ['90']
ancestors = ['Anaplasmataceae']

for node, anc in zip(nodes, ancestors):
    sel = events[events['node'] == node]
    del sel['node']
    for e in ['duplications', 'transfers', 'losses', 'originations', 'copies']:
        sel[e] = sel[e].apply(lambda x: modf(x)[1] + int(modf(x)[0] >= 0.3))
    sel = sel[(sel['losses'] > 0.3) | (sel['copies'] > 0.3)]
    sel['eggNOG annot'] = sel['cluster'].apply(lambda x: get_annot(x, 'eggNOG annot'))
    sel['eggNOG cat'] = sel['cluster'].apply(lambda x: get_annot(x, 'eggNOG cat'))
    sel.index = sel['cluster']
    del sel['cluster']
    sel.to_csv('{}.csv'.format(anc), header=True, sep='\t', index=True)
