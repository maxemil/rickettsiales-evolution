for faa in *.faa 
do
   base=${faa%%.*}
   mafft-einsi --thread 20 $faa > $base.aln
   trimal -in $base.aln -out $base.99.aln -gt 0.01
   iqtree2 -s $base.99.aln -pre $base.99 -nt 20 -m LG+R8 -bb 1000
 done


for gene cog in $(cat gene2cog.csv );
# for gene cog in "PilB" "COG2804";
do
  wget http://eggnogapi5.embl.de/nog_data/file/fasta/$cog -O $cog.faa.gz
  gunzip $cog.faa.gz
  python remove_Rickettsiales.py $cog.faa "$cog"_noRick.faa
  
  diamond makedb --in "$cog"_noRick.faa --db $cog.dmnd --threads 10
  diamond blastp -d $cog.dmnd \
                 -q ../faa_subunits/$gene.faa \
                 --quiet --threads 20 \
                 -o "$gene"_EggNOG.out \
                 --outfmt 5 --top 50 --ultra-sensitive 
  /local/two/Software/misc-scripts/extract_blast_hits.py -i "$gene"_EggNOG.out -o "$gene"_EggNOG.hits -f $cog.faa -m blastp

  cd-hit -i "$gene"_EggNOG.hits -o "$gene"_EggNOG.80.hits -c 0.80 -T 1 -M 8000
  cat "$gene"_EggNOG.80.hits ../faa_subunits/$gene.faa > "$gene"_eggnog_ref.faa

  mafft-einsi --thread 20 "$gene"_eggnog_ref.faa > "$gene"_eggnog_ref.aln
  trimal -in "$gene"_eggnog_ref.aln -out "$gene"_eggnog_ref.99.aln -gt 0.01

  iqtree2 -s "$gene"_eggnog_ref.99.aln -pre "$gene"_eggnog_ref.99 -nt 20 -m MFP -mset LG -bb 1000
  python /local/two/Software/misc-scripts/color_tree.py -i "$gene"_eggnog_ref.99.treefile \
                                                        -o "$gene"_eggnog_ref.99.nex \
                                                        -c taxa_colors.txt --taxonomy
done
