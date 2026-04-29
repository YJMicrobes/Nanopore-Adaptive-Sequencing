#!/bin/bash
set -euo pipefail

##############################################################################
## 0. USER CONFIGURATION (EDIT HERE ONLY)
##############################################################################
#Name: Jing Yuan
#Date: 2026-4-29
#Project: Nanopore adaptive sequencing of mosquito microbiome
#Institute: Connecticut Agricultural Experiment Station

########################
# INPUTS
########################
RAW_DIR="1_pod5"
MOSQUITO_GENOME="mosquito_genome.fa"


########################
# OUTPUT DIRECTORIES
########################
DEMUX_DIR="2_demultiplex"
NPLOT_DIR="3_nanoplot"
FILTER_DIR="4_nanofilt"
KAIJU_OUT="5_kaiju_reads"
ASSEMBLY_DIR="6_assembly"
KRAKEN_DIR="7_kraken"
PROKKA_DIR="8_prokka"
EGGNOG_DIR="9_eggnog"
BIN_DIR="10_binning"
DASTOOL_DIR="11_dastool"
CHECKM_DIR="12_checkm"
GTDB_DIR="13_gtdbtk"
BLAST_DIR="14_host_check"
SPLIT_DIR="15_contig_split"
ADAPT_DIR="16_adaptive_targets"

########################
# THREADS
########################
THREADS=30

########################
# DATABASE PATHS
########################
KAIJU_DB="$HOME/kaijudb"
KRAKEN_DB="$HOME/standard_database"
EGGNOG_DB="$HOME/eggnogdb"
GTDBTK_DATA_PATH="$HOME/gtdb_database"


##############################################################################
## 1. HPC MODULES (VERSION DEPENDS ON YOUR SERVER)
##############################################################################
# NOTE: versions may differ per cluster → check via `module avail`

module load Dorado
module load kraken
module load flye
module load prokka
module load eggnog-mapper
module load bowtie2
module load samtools
module load metabat
module load maxbin
module load metacoag
module load das_tool
module load CheckM
module load blast
module load edirect
module load seqtk
module load kaiju
module load nanoplot
module load NanoFilt


##############################################################################
## 2. BASECALLING + DEMULTIPLEXING
##############################################################################
mkdir -p ${DEMUX_DIR}

dorado basecaller $HOME/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    ${RAW_DIR}/ \
    --kit-name SQK-NBD111-24 \
    --no-trim \
    > ${DEMUX_DIR}/reads.bam

dorado demux \
    --kit-name SQK-NBD114-24 \
    --emit-fastq \
    --output-dir ${DEMUX_DIR} \
    ${DEMUX_DIR}/reads.bam \
    -t ${THREADS}

##############################################################################
## 3. QC + FILTERING
##############################################################################
mkdir -p ${NPLOT_DIR} ${FILTER_DIR}

NanoPlot --fastq ${DEMUX_DIR}/*.fastq -o ${NPLOT_DIR}/

for fq in ${DEMUX_DIR}/b*.fastq; do
    base=$(basename ${fq} .fastq)
    NanoFilt -q 10 -l 1000 ${fq} > ${FILTER_DIR}/${base}_filt.fastq
done

##############################################################################
## 4. READ TAXONOMY (KAIJU)
##############################################################################
mkdir -p ${KAIJU_OUT}

for fq in ${FILTER_DIR}/*.fastq; do
    base=$(basename ${fq} .fastq)

    kaiju \
        -t ${KAIJU_DB}/nodes.dmp \
        -f ${KAIJU_DB}/refseq_nr/kaiju_db_refseq_nr.fmi \
        -i ${fq} \
        -o ${KAIJU_OUT}/${base}.out \
        -a greedy -e 3 -m 11 -s 65 -l 7 -E 0.01
done

##############################################################################
## 5. ASSEMBLY (metaFlye)
##############################################################################
mkdir -p ${ASSEMBLY_DIR}

cat ${FILTER_DIR}/*.fastq > ${FILTER_DIR}/merged.fastq

flye \
    --nano-raw ${FILTER_DIR}/merged.fastq \
    --meta \
    --out-dir ${ASSEMBLY_DIR} \
    -t ${THREADS}

##############################################################################
## 6. CONTIG TAXONOMY (Kraken2)
##############################################################################
mkdir -p ${KRAKEN_DIR}

asm="${ASSEMBLY_DIR}/assembly.fasta"
base=$(basename "${asm}" .fasta)

kraken2 \
    --db ${KRAKEN_DB} \
    --output ${KRAKEN_DIR}/${base}.out \
    --report ${KRAKEN_DIR}/${base}.report \
    --use-names \
    --confidence 0.1 \
    --threads ${THREADS} \
    "${asm}"

python scripts/kraken_summary.py

## KRAKEN SUMMARY TABLE (Python)
#import pandas as pd
#import glob
#import os
#files = glob.glob("*_report.txt")
#dfs = []
#for f in files:
#    sample = os.path.basename(f).replace("_report.txt", "")   
#    rows = []
#    with open(f) as fh:
#        for line in fh:
#            parts = line.strip().split(None, 5)
#            if len(parts) == 6:
#                percent, _, _, _, _, name = parts
#                rows.append([name, float(percent)])
#    df = pd.DataFrame(rows, columns=["name", sample]).set_index("name")
#    dfs.append(df)
#summary = pd.concat(dfs, axis=1).fillna(0)
#summary.to_csv("kraken_summary_table.csv")


##############################################################################
## 7. FUNCTIONAL ANNOTATION (Prokka)
##############################################################################
mkdir -p ${PROKKA_DIR}

prokka \
    --cpus ${THREADS} \
    --outdir ${PROKKA_DIR}/${base} \
    --prefix ${base} \
    ${asm}

##############################################################################
## 8. FUNCTIONAL ANNOTATION (eggNOG)
##############################################################################
mkdir -p ${EGGNOG_DIR}

export OMP_NUM_THREADS=10

find ${PROKKA_DIR} -name "*.faa" | \
parallel -j 3 emapper.py \
    -i {} \
    -o ${EGGNOG_DIR}/{/.}_eggnog \
    --cpu 10 \
    -m diamond \
    --dmnd_algo ctg

##############################################################################
## 9. READ MAPPING
##############################################################################
mkdir -p ${BIN_DIR}

bowtie2-build ${asm} ${BIN_DIR}/index

bowtie2 \
    -x ${BIN_DIR}/index \
    -U ${FILTER_DIR}/merged.fastq \
    -S ${BIN_DIR}/mapped.sam \
    -p ${THREADS}

samtools view -bS ${BIN_DIR}/mapped.sam | \
samtools sort -o ${BIN_DIR}/mapped.bam

samtools index ${BIN_DIR}/mapped.bam

##############################################################################
## 10. BINNING
##############################################################################
jgi_summarize_bam_contig_depths \
    --outputDepth ${BIN_DIR}/depth.txt \
    ${BIN_DIR}/mapped.bam

metabat2 -i ${asm} -a ${BIN_DIR}/depth.txt -o ${BIN_DIR}/metabat/bin -t ${THREADS}

run_MaxBin.pl -contig ${asm} -abund ${BIN_DIR}/depth.txt -out ${BIN_DIR}/maxbin/bin -thread ${THREADS}

metacoag \
    --assembler flye \
    --graph ${ASSEMBLY_DIR}/assembly_graph.gfa \
    --contigs ${asm} \
    --paths ${ASSEMBLY_DIR}/assembly_info.txt \
    --abundance ${BIN_DIR}/depth.txt \
    --output ${BIN_DIR}/metacoag

##############################################################################
## 11. BIN REFINEMENT (DAS Tool)
##############################################################################
mkdir -p ${DASTOOL_DIR}

DAS_Tool \
    -i bins.tsv \
    -l metabat,maxbin,metacoag \
    -c ${asm} \
    -o ${DASTOOL_DIR}/output \
    --write_bins \
    --threads ${THREADS}


##############################################################################
## 12. BIN QUALITY (CheckM)
##############################################################################
mkdir -p ${CHECKM_DIR}

checkm lineage_wf \
    -x fa \
    ${DASTOOL_DIR}/output_DASTool_bins/ \
    ${CHECKM_DIR} \
    -t ${THREADS}


##############################################################################
## 13. TAXONOMY (GTDB-Tk)
##############################################################################
mkdir -p ${GTDB_DIR}

export GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}

gtdbtk classify_wf \
    --genome_dir ${DASTOOL_DIR}/output_DASTool_bins \
    --out_dir ${GTDB_DIR} \
    --extension fa \
    --cpus ${THREADS} \
    --skip_ani_screen

##############################################################################
## 14. HOST CONTAMINATION CHECK (BLAST)
##############################################################################
mkdir -p ${BLAST_DIR}

blastn \
    -query ${asm} \
    -db nt \
    -remote \
    -out ${BLAST_DIR}/host.out \
    -outfmt 6 \
    -max_target_seqs 10 \
    -evalue 1e-5


##############################################################################
## 15. ADAPTIVE SEQUENCING TARGETS
##############################################################################
mkdir -p ${ADAPT_DIR}

cat mosquito_matched_contigs.fa ${MOSQUITO_GENOME} > ${ADAPT_DIR}/depletion.fa
cat non_mosquito_microbial_contigs.fa > ${ADAPT_DIR}/enrichment.fa

echo "Pipeline completed."


