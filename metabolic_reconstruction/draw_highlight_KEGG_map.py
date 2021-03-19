from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
import os
import time
import pandas as pd

def draw_pathway_map(pathway, outfile):
    pathway = KGML_parser.read(kegg_get(pathway, "kgml"))
    present = False
    for p in pathway.orthologs:
        p.graphics[0].bgcolor = "#FFFFFF"
        if len(p.graphics) > 1:
            print("more than one graphics item for {}".format(p.name))
        ks = p.name.replace("ko:","").split()
        if any([k in k_list for k in ks]):
                p.graphics[0].bgcolor = "#4585D7"
                present = True
        # elif any([k in k_list for k in ks]):
        #         p.graphics[0].bgcolor = "#C7DAF3"
        #         present = True
    if present:
        canvas = KGMLCanvas(pathway, import_imagemap=True)
        canvas.draw(outfile)

df = pd.read_csv('clst_125_annotations.csv', sep='\t', index_col=0)
# df = df[df['bin_125'] > 0]
k_list = [ko for kos in set(df['KO'].dropna()) for ko in kos.split(';')]

maps = ['ko00010', 'ko00030', 'ko00020', 'ko00190', 'ko00071', 'ko00061',
        'ko00900', 'ko00630', 'ko00230', 'ko00240', 'ko00740', 'ko00770',
        'ko00790', 'ko00780', 'ko00730', 'ko00760', 'ko00860']
pathways = ['Glycolysis', 'Pentose_Phosphate', 'Citrate_cycle', 
            'Oxidative_phosphorylation', 'Fatty_acid_degradation', 
            'Fatty_acid_biosynthesis', 'Terpenoid_biosynthesis',
            'Glyoxylate_dicarboxylate', 'Purine_metabolism', 
            'Pyrimidine_metabolism', 'Riboflavin_metabolism', 
            'Pantothenate_CoA', 'Folate', 'Biotin', 'Thiamine', 
            'Nicotinate_nicotinamide', 'Porphyrin']
for map, path in zip(maps, pathways):
        draw_pathway_map(map, "{}.pdf".format(path))
