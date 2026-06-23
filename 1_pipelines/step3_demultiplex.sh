##############################################################################
## 3. demultiplexing
##############################################################################


mkdir 2_demultiplex


#!/bin/bash
#SBATCH --job-name=basecalling
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=150G
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


#No need to use GPU

module load Dorado/0.9.0
dorado demux --kit-name SQK-NBD114-96 --emit-fastq --output-dir 2_demultiplex/ 1_basecall/hgt.bam -t 30