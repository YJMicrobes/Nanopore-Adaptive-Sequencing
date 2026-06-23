##############################################################################
## 9. binning: metacoag (per group)
##############################################################################

mkdir -p 9_metacoag_out


source ~/.bashrc
conda activate metacoag_env


groups=(colony_no_subtraction colony_subtraction H G T)


for g in "${groups[@]}"; do

    echo "Running MetaCoAG for ${g}"


    mkdir -p 9_metacoag_out/${g}


    ##########################################################################
    # Inputs 
    ##########################################################################

    CONTIGS=6_flye_output/${g}/assembly.fasta
    GRAPH=6_flye_output/${g}/assembly_graph.gfa
    PATHS=6_flye_output/${g}/assembly_info.txt
    ABUNDANCE=8_metabat_out/${g}.depth.txt



    ##########################################################################
    # Run MetaCoAG
    ##########################################################################

    metacoag \
        --assembler flye \
        --graph "$GRAPH" \
        --contigs "$CONTIGS" \
        --paths "$PATHS" \
        --abundance "$ABUNDANCE" \
        --output 9_metacoag_out/${g}


done