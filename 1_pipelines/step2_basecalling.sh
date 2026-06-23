
##############################################################################
## 2. basecalling
##############################################################################


mkdir 1_basecall

#check gpu
sinfo -o "%P %G %D %N"


#!/bin/bash
#SBATCH --job-name=basecalling
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --partition=general
#SBATCH --gres=gpu:A100:1
#SBATCH --qos=general
#SBATCH --mem=150G
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


#Pre-work 1:
#download the module first
#dorado download 

#Pre-work 2:
#check your nanopore library kit name

#Pre-work 3:
#check your module Dorado version

module load Dorado/0.9.0
dorado basecaller ~/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/ --kit-name SQK-NBD111-96 --no-trim > 1_basecall/hgt.bam