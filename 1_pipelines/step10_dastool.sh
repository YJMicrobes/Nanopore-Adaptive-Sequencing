##############################################################################
## 10. DAS Tool (per group integration)
##############################################################################

source ~/.bashrc
conda activate das_tool_env


groups=(colony_no_subtraction colony_subtraction H G T)


for g in "${groups[@]}"; do

    echo "Running DAS Tool for ${g}"


    mkdir -p 10_dastool/${g}/bins
    mkdir -p 10_dastool/${g}/output



    ##########################################################################
    # Step 1: link bins
    ##########################################################################

    rm -f 10_dastool/${g}/bins/*


    # ---- MaxBin bins ----
    for f in 7_maxbin_out/${g}/${g}_maxbin.*.fasta; do
        [ -e "$f" ] || continue
        base=$(basename "$f" .fasta)
        ln -s "$(realpath "$f")" 10_dastool/${g}/bins/maxbin_${base}.fa
    done



    # ---- MetaBAT bins ----
    for f in 8_metabat_out/${g}_bins/bin.*.fa; do
        [ -e "$f" ] || continue
        base=$(basename "$f" .fa)
        ln -s "$(realpath "$f")" 10_dastool/${g}/bins/metabat_${base}.fa
    done



    # ---- MetaCoAG bins ----
    for f in 9_metacoag_out/${g}/*.fa; do
        [ -e "$f" ] || continue
        base=$(basename "$f" .fa)
        ln -s "$(realpath "$f")" 10_dastool/${g}/bins/metacoag_${base}.fa
    done



    ##########################################################################
    # Step 2: build contig-to-bin table
    ##########################################################################

    cd 10_dastool/${g}/bins || exit 1

    > contig_to_bin.tsv

    for f in *.fa; do

        bin=$(basename "$f" .fa)

        awk -v bin="$bin" '
        /^>/ {
            gsub(/^>/, "", $1);
            print $1 "\t" bin
        }' "$f" >> contig_to_bin.tsv

    done


    cd ../../..



    ##########################################################################
    # Step 3: run DAS Tool
    ##########################################################################

    DAS_Tool \
        -i 10_dastool/${g}/bins/contig_to_bin.tsv \
        -c 6_flye_output/${g}/assembly.fasta \
        -o 10_dastool/${g}/output/${g}_dastool \
        --write_bin_evals \
        --write_bins \
        --write_unbinned \
        --threads 8 \
        --score_threshold 0.5


done