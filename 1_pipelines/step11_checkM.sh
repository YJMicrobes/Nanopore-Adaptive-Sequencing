##############################################################################
## 11. CheckM (per group MAG quality assessment)
##############################################################################

source ~/.bashrc
conda activate checkm_env


mkdir -p 11_checkm_out



groups=(colony_no_subtraction colony_subtraction H G T)



for g in "${groups[@]}"; do

    echo "Running CheckM for ${g}"


    ##########################################################################
    # Input bins from DAS Tool
    ##########################################################################

    BIN_DIR=10_dastool/${g}/output/${g}_dastool_DASTool_bins
    OUT_DIR=11_checkm_out/${g}

    mkdir -p "$OUT_DIR"
    
    ##########################################################################
    # Run CheckM
    ##########################################################################

    checkm lineage_wf \
        -x fa \
        -t 48 \
        "$BIN_DIR" \
        "$OUT_DIR"


done