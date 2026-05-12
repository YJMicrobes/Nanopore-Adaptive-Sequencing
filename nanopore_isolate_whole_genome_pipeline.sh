#!/bin/bash
set -euo pipefail

##############################################################################
## 0. USER CONFIGURATION (EDIT HERE ONLY)
##############################################################################
# Name: Jing Yuan
# Date: 2026-05-12
# Project: Nanopore whole-genome sequencing of microbial isolates
# Institute: Connecticut Agricultural Experiment Station
#
# Description:
# End-to-end Oxford Nanopore pipeline for bacterial isolate whole-genome
# sequencing, including:
#   1. Basecalling (Dorado)
#   2. Demultiplexing
#   3. Quality control and read filtering
#   4. Genome assembly (Flye)
#   5. Genome quality assessment (CheckM)
#   6. Optional assembly validation (BUSCO)
#   7. Genome annotation (Prokka)
#   8. Taxonomic classification (GTDB-Tk)
#   9. Genome comparison (FastANI)

########################
# INPUTS
########################
RAW_DIR="1_pod5"                       # Directory containing POD5 files
MODEL="$HOME/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

########################
# BARCODING KITS
########################
BASECALL_KIT="SQK-NBD111-24"
DEMUX_KIT="SQK-NBD114-24"

########################
# OUTPUT DIRECTORIES
########################
DEMUX_DIR="2_demultiplex"
NPLOT_DIR="3_nanoplot"
FILTER_DIR="4_nanofilt"
ASSEMBLY_DIR="5_flye"
CHECKM_DIR="6_checkm"
BUSCO_DIR="7_busco"
PROKKA_DIR="8_prokka"
GTDB_DIR="9_gtdbtk"
FASTANI_DIR="10_fastani"

########################
# PARAMETERS
########################
THREADS=30
GENOME_SIZE="5m"       # Estimated bacterial genome size
MIN_QUALITY=10
MIN_LENGTH=150

########################
# DATABASE PATHS
########################
GTDBTK_DATA_PATH="$HOME/gtdb_database"

########################
# OPTIONAL QUERY GENOMES
########################
# Place genomes to compare against the isolate assemblies in this directory.
# Example: previously assembled metagenome bins.
QUERY_DIR="query_genomes"
QUERY_EXTENSION="fa"

##############################################################################
## 1. HPC MODULES / CONDA ENVIRONMENTS
##############################################################################
# Adjust based on your cluster environment.

# Example module loading:
# module load Dorado
# module load flye/2.9.5
# module load busco/5.7.0
# module load parallel
#
# Example conda environments:
# source ~/.bashrc
# conda activate nanoplot_env
# conda activate checkm_env
# conda activate prokka_env
# conda activate gtdbtk_env
# conda activate fastANI_env

##############################################################################
## 2. CREATE OUTPUT DIRECTORIES
##############################################################################
mkdir -p \
    "${DEMUX_DIR}" \
    "${NPLOT_DIR}" \
    "${FILTER_DIR}" \
    "${ASSEMBLY_DIR}" \
    "${CHECKM_DIR}" \
    "${BUSCO_DIR}" \
    "${PROKKA_DIR}" \
    "${GTDB_DIR}" \
    "${FASTANI_DIR}"

##############################################################################
## 3. BASECALLING
##############################################################################
echo "=================================================================="
echo "STEP 1: Basecalling with Dorado"
echo "=================================================================="

dorado basecaller "${MODEL}" "${RAW_DIR}/" \
    --kit-name "${BASECALL_KIT}" \
    --no-trim \
    > "${DEMUX_DIR}/reads.bam"

##############################################################################
## 4. DEMULTIPLEXING
##############################################################################
echo "=================================================================="
echo "STEP 2: Demultiplexing"
echo "=================================================================="

dorado demux \
    --kit-name "${DEMUX_KIT}" \
    --emit-fastq \
    --output-dir "${DEMUX_DIR}" \
    "${DEMUX_DIR}/reads.bam" \
    -t "${THREADS}"

##############################################################################
## 5. QUALITY CONTROL (NanoPlot)
##############################################################################
echo "=================================================================="
echo "STEP 3: Quality Control with NanoPlot"
echo "=================================================================="

# Activate NanoPlot environment if needed:
# source ~/.bashrc
# conda activate nanoplot_env

NanoPlot \
    --fastq "${DEMUX_DIR}"/*.fastq \
    -o "${NPLOT_DIR}"

##############################################################################
## 6. READ FILTERING (NanoFilt)
##############################################################################
echo "=================================================================="
echo "STEP 4: Filtering reads with NanoFilt"
echo "=================================================================="

for fq in "${DEMUX_DIR}"/*.fastq; do
    sample=$(basename "${fq}" .fastq)

    echo "Filtering ${sample}..."

    NanoFilt \
        -q "${MIN_QUALITY}" \
        -l "${MIN_LENGTH}" \
        "${fq}" \
        > "${FILTER_DIR}/${sample}_filt.fastq"
done

##############################################################################
## 7. GENOME ASSEMBLY (Flye)
##############################################################################
echo "=================================================================="
echo "STEP 5: Genome Assembly with Flye"
echo "=================================================================="

# module load flye/2.9.5

for fq in "${FILTER_DIR}"/*_filt.fastq; do
    sample=$(basename "${fq}" _filt.fastq)
    outdir="${ASSEMBLY_DIR}/${sample}"

    echo "Assembling ${sample}..."

    flye \
        --nano-raw "${fq}" \
        --genome-size "${GENOME_SIZE}" \
        --threads "${THREADS}" \
        --out-dir "${outdir}"

    echo "Completed assembly for ${sample}"
done

##############################################################################
## 8. GENOME QUALITY ASSESSMENT (CheckM)
##############################################################################
echo "=================================================================="
echo "STEP 6: Genome Quality Assessment with CheckM"
echo "=================================================================="

# source ~/.bashrc
# conda activate checkm_env

for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")
    outdir="${CHECKM_DIR}/${sample}"

    echo "Running CheckM for ${sample}..."

    mkdir -p "${outdir}"

    checkm lineage_wf \
        -x fasta \
        "$(dirname "${assembly}")" \
        "${outdir}" \
        -t "${THREADS}"

    checkm qa \
        "${outdir}/lineage.ms" \
        "${outdir}" \
        -f "${outdir}/checkm_summary.tsv"

    echo "Completed CheckM for ${sample}"
done

##############################################################################
## 9. OPTIONAL ASSEMBLY VALIDATION (BUSCO)
##############################################################################
echo "=================================================================="
echo "STEP 7: Optional BUSCO Analysis"
echo "=================================================================="

# Uncomment to run.
# module load busco/5.7.0
#
# for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
#     sample=$(basename "$(dirname "${assembly}")")
#
#     busco \
#         -i "${assembly}" \
#         -l bacteria_odb10 \
#         -o "${sample}" \
#         -m genome \
#         --out_path "${BUSCO_DIR}" \
#         -f \
#         -c "${THREADS}"
# done

##############################################################################
## 10. GENOME ANNOTATION (Prokka)
##############################################################################
echo "=================================================================="
echo "STEP 8: Genome Annotation with Prokka"
echo "=================================================================="

# source ~/.bashrc
# conda activate prokka_env

for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")
    outdir="${PROKKA_DIR}/${sample}"

    echo "Annotating ${sample}..."

    prokka \
        --cpus "${THREADS}" \
        --outdir "${outdir}" \
        --prefix "${sample}" \
        "${assembly}"

    echo "Completed Prokka for ${sample}"
done

##############################################################################
## 11. TAXONOMIC CLASSIFICATION (GTDB-Tk)
##############################################################################
echo "=================================================================="
echo "STEP 9: Taxonomic Classification with GTDB-Tk"
echo "=================================================================="

# source ~/.bashrc
# conda activate gtdbtk_env

mkdir -p "${GTDB_DIR}/genomes"

# Copy assemblies and rename to .fa
for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")
    cp "${assembly}" "${GTDB_DIR}/genomes/${sample}.fa"
done

export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH}"

gtdbtk classify_wf \
    --genome_dir "${GTDB_DIR}/genomes" \
    --out_dir "${GTDB_DIR}/results" \
    --extension fa \
    --cpus "${THREADS}" \
    --skip_ani_screen \
    --force

##############################################################################
## 12. GENOME COMPARISON (FastANI)
##############################################################################
echo "=================================================================="
echo "STEP 10: Genome Comparison with FastANI"
echo "=================================================================="

# source ~/.bashrc
# conda activate fastANI_env

mkdir -p "${FASTANI_DIR}/references"

# Prepare reference genome list from assemblies
for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")
    cp "${assembly}" "${FASTANI_DIR}/references/${sample}.fa"
done

cd "${FASTANI_DIR}"

ls references/*.fa > reference_list.txt

# Compare external query genomes (if provided)
if [[ -d "${QUERY_DIR}" ]] && ls "${QUERY_DIR}"/*.${QUERY_EXTENSION} >/dev/null 2>&1; then
    echo "Running FastANI comparisons..."

    ls "${QUERY_DIR}"/*.${QUERY_EXTENSION} > query_list.txt

    fastANI \
        --ql query_list.txt \
        --rl reference_list.txt \
        -o fastani_output.tsv

    echo "FastANI results written to ${FASTANI_DIR}/fastani_output.tsv"
else
    echo "No query genomes found in ${QUERY_DIR}/"
    echo "Skipping FastANI comparison."
fi

cd - >/dev/null

##############################################################################
## 13. OUTPUT SUMMARY
##############################################################################
echo "=================================================================="
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "=================================================================="
echo ""
echo "Output directories:"
echo "  ${DEMUX_DIR}     - Demultiplexed FASTQ files"
echo "  ${NPLOT_DIR}     - NanoPlot QC reports"
echo "  ${FILTER_DIR}    - Filtered FASTQ files"
echo "  ${ASSEMBLY_DIR}  - Flye assemblies"
echo "  ${CHECKM_DIR}    - CheckM quality reports"
echo "  ${BUSCO_DIR}     - BUSCO completeness reports (optional)"
echo "  ${PROKKA_DIR}    - Genome annotations"
echo "  ${GTDB_DIR}      - GTDB-Tk taxonomy assignments"
echo "  ${FASTANI_DIR}   - ANI comparison results"
echo ""
echo "Recommended downstream analyses:"
echo "  - Functional annotation (eggNOG-mapper)"
echo "  - AMR detection (ABRicate)"
echo "  - Virulence gene screening"
echo "  - Pan-genome analysis (Roary or Panaroo)"
echo "  - Phylogenetic tree construction"
echo "  - NCBI GenBank submission"
echo ""
echo "Done."