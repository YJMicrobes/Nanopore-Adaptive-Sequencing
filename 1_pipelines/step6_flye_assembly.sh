##############################################################################
## 6. Contig assembly - Flye (5 groups separately)
##############################################################################

mkdir -p 6_flye_output


##############################################################################
## Step 1. Merge reads within each group
##############################################################################

groups=(
colony_no_subtraction
colony_subtraction
H
G
T
)


for group in "${groups[@]}"; do

    echo "Processing ${group}"

    cat 5_nanofilt/${group}/*_filt.fastq \
        > 5_nanofilt/${group}_merged.fastq

    echo "Created: 5_nanofilt/${group}_merged.fastq"

done



##############################################################################
## Step 2. Run Flye separately for each group
##############################################################################

#!/bin/bash
#SBATCH --job-name=flye_5groups
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mem=500G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


module load flye/2.9.5


for group in colony_no_subtraction colony_subtraction H G T

do

    echo "Starting Flye assembly for ${group}"


    flye \
      --nano-raw 5_nanofilt/${group}_merged.fastq \
      --meta \
      --out-dir 6_flye_output/${group} \
      -t 48


    echo "Finished ${group}"

done