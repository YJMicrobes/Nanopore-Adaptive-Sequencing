##############################################################################
## 7. binning: maxbin
##############################################################################

mkdir -p 7_maxbin_out


#!/bin/bash
#SBATCH --job-name=maxbin
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=150G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


source ~/.bashrc
conda activate maxbin_env


groups=(colony_no_subtraction colony_subtraction H G T)


for g in "${groups[@]}"; do

    echo "Running MaxBin for ${g}"

    mkdir -p "7_maxbin_out/${g}"


    run_MaxBin.pl \
        -contig 6_flye_output/${g}/assembly.fasta \
        -reads 5_nanofilt/${g}/*_filt.fastq \
        -out 7_maxbin_out/${g}/${g}_maxbin \
        -thread 48

done

mkdir 7_maxbin_out

#check bins
ls -lhA 7_maxbin_out
