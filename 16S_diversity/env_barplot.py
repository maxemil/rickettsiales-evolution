import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.DataFrame({
        'marine':8, 
        'lake':7, 
        'groundwater':5, 
        'host associated':2,
        'glacier ice':1, 
        'fish culture water':1}, index=['para'])
df.loc['para'] = df.loc['para']/24 * 100
df = df.append({'marine':100}, ignore_index=True)
df.index = ['para', 'miti']
df.fillna(0, inplace=True)

fig, ax = plt.subplots(figsize=(10, 4))
axs = df.T.plot(kind='pie', stacked=True, subplots=True, ax=ax)
my_circle=plt.Circle( (0,0), 0.7, color='white')
axs[0].add_artist(my_circle)
my_circle=plt.Circle( (0,0), 0.7, color='white')
axs[1].add_artist(my_circle)
plt.tight_layout()
plt.savefig("16S_env_donut.pdf")
