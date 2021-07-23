for f in genomes/*.fna
do 
  base=$(basename ${f%%.fna})
  python /local/two/Software/misc-scripts/parse_barrnap.py -l <(singularity exec -B /local ../../picozoa-workflow/picozoa.sif barrnap --reject 0.2 $f) -r $f -o 23S_sequences/$base.23S.fasta -f 23S
done
