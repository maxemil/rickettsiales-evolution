nextflow run /local/two/Software/misc-scripts/rp15-on-contigs.nf --assemblies 'assemblies/*' \
        --proteomes 'proteomes/*' \
        --genomes 'genomes/*' \
        --references 'rp15_bacteria_alpha_mito_db_tax/*' \
        --phylo_method iqtree \
        -with-singularity rp15.sif \
        --prune_percentage 0.4 \
        --cpus 40 \
        -resume
