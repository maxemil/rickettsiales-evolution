
spec = ["UBA6178", "UBA6189", "UBA6149", "TOBG_RS-372",  "bin_125", "UBA6187", "TARA_ANE_MAG_00011"]

PG_clust = set()
for sp in spec:
    f_handle = "42_effectors_repeats/outfiles/{}.out".format(sp)
    PG_prots = [p.split()[0] for p in open(f_handle) if 'PG_binding' in p]

    f_handle = '28_eggnog_annotation_updated/annotations_aproNOG/{}.emapper.annotations'.format(sp)
    for line in open(f_handle):
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[0] in PG_prots:
            PG_clust.add(line[10].split('|')[0])


for line in open('34_ALE_updated_clusters/marineRick-aphaNOG-wDeia-0.3.csv'):
    line = line.split('\t')
    if line[0] in PG_clust and line[1] == 'LRCA':
        print(line)
