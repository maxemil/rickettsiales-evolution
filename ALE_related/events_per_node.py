import pandas as pd
from math import modf
import numpy as np

header = ['species_tree', 'cluster', 'node', 'duplications',
          'transfers', 'losses', 'originations', 'copies']
events = pd.read_csv('marineRick-aphaNOG-wDeia/events.txt', 
                            sep='\t', names=header)
del events['species_tree']
events['cluster'] = events['cluster'].str.replace('.clean.ufboot.ale', '')
events['cluster'] = events['cluster'].str.replace('unlab_60_90', 'silix')

for e in ['duplications', 'transfers', 'losses', 'originations', 'copies']:
    events[e] = events[e].apply(lambda x: int(modf(x)[1] + int(modf(x)[0] >= 0.3)))

names = {}
names['104'] = 'root'
names['103'] = 'LRCA'
names['101'] = 'LhRCA'

events['node'] = events['node'].apply(lambda x: names[x] if x in names else x)

events.to_csv("marineRick-aphaNOG-wDeia-0.3.csv", sep='\t', header=True, index=False)

names.update({line.split()[1]:line.split()[0] for line in open('species_map.txt')})

events['node'] = events['node'].apply(lambda x: names[x] if x in names else x)

events_summary = pd.pivot_table(events, values=['duplications', 'transfers', 'losses', 'originations', 'copies'],
    index=['node'], 
    aggfunc=np.sum)

events_summary.to_csv('events_summary.tsv', sep='\t', index=True, header=True)        
events_summary.to_excel('events_summary.xlsx', engine='openpyxl', index=True, header=True)        


# event.to_excel('{}_{}_raw.xlsx'.format(e, tree_suffix), engine='openpyxl', 
#                 header=True, index=True)
