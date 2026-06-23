##############################################################################
## 8. binning: MetaBAT2 (per group)
##############################################################################

mkdir -p 8_metabat_out


#!/bin/bash
#SBATCH --job-name=metabat
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mem=500G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


module load minimap2/2.28
module load samtools/1.19.2
module load metabat/2.12.1



groups=(colony_no_subtraction colony_subtraction H G T)



for g in "${groups[@]}"; do

    echo "Processing ${g}"


    ##########################################################################
    # Step 1: mapping reads to assembly
    ##########################################################################

    minimap2 \
      --split-prefix tmp \
      -t 48 \
      -a \
      -x map-ont \
      6_flye_output/${g}/assembly.fasta \
      5_nanofilt/${g}/*_filt.fastq | \
    samtools sort -@ 48 -o 8_metabat_out/${g}.sorted.bam


    samtools index 8_metabat_out/${g}.sorted.bam



    ##########################################################################
    # Step 2: compute coverage depth
    ##########################################################################

    jgi_summarize_bam_contig_depths \
      --outputDepth 8_metabat_out/${g}.depth.txt \
      8_metabat_out/${g}.sorted.bam



    ##########################################################################
    # Step 3: MetaBAT2 binning
    ##########################################################################

    mkdir -p 8_metabat_out/${g}_bins


    metabat2 \
      -i 6_flye_output/${g}/assembly.fasta \
      -a 8_metabat_out/${g}.depth.txt \
      -o 8_metabat_out/${g}_bins/bin \
      --minContig 1500 \
      --unbinned

done