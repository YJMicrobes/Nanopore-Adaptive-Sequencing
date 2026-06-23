##############################################################################
## Nanopore Whole-Genome Sequencing Pipeline (Isolates)
## Author: Jing Yuan
##############################################################################

#!/bin/bash
set -euo pipefail


##############################################################################
## 0. USER CONFIGURATION
##############################################################################

RAW_DIR="1_pod5"
MODEL="$HOME/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

BASECALL_KIT="SQK-NBD111-24"
DEMUX_KIT="SQK-NBD114-24"

THREADS=30
GENOME_SIZE="5m"
MIN_QUALITY=10
MIN_LENGTH=150

##############################################################################
## 1. OUTPUT DIRECTORIES
##############################################################################

DEMUX_DIR="2_demultiplex"
NPLOT_DIR="3_nanoplot"
FILTER_DIR="4_nanofilt"
ASSEMBLY_DIR="5_flye"
CHECKM_DIR="6_checkm"
BUSCO_DIR="7_busco"
PROKKA_DIR="8_prokka"
GTDB_DIR="9_gtdbtk"

mkdir -p \
    "${DEMUX_DIR}" \
    "${NPLOT_DIR}" \
    "${FILTER_DIR}" \
    "${ASSEMBLY_DIR}" \
    "${CHECKM_DIR}" \
    "${BUSCO_DIR}" \
    "${PROKKA_DIR}" \
    "${GTDB_DIR}/genomes" \
    "${GTDB_DIR}/results"


##############################################################################
## 2. BASECALLING (Dorado)
##############################################################################

echo "Step 1: Basecalling"

dorado basecaller "${MODEL}" "${RAW_DIR}/" \
    --kit-name "${BASECALL_KIT}" \
    --no-trim \
    --emit-fastq \
    > "${DEMUX_DIR}/reads.fastq"


##############################################################################
## 3. DEMULTIPLEXING
##############################################################################

echo "Step 2: Demultiplexing"

dorado demux \
    --kit-name "${DEMUX_KIT}" \
    --output-dir "${DEMUX_DIR}" \
    "${DEMUX_DIR}/reads.fastq" \
    -t "${THREADS}"


##############################################################################
## 4. QUALITY CONTROL
##############################################################################

echo "Step 3: NanoPlot QC"

NanoPlot \
    --fastq "${DEMUX_DIR}"/*.fastq \
    --outdir "${NPLOT_DIR}"


##############################################################################
## 5. READ FILTERING
##############################################################################

echo "Step 4: Filtering reads"

for fq in "${DEMUX_DIR}"/*.fastq; do
    sample=$(basename "${fq}" .fastq)

    NanoFilt \
        -q "${MIN_QUALITY}" \
        -l "${MIN_LENGTH}" \
        "${fq}" \
        > "${FILTER_DIR}/${sample}_filt.fastq"
done


##############################################################################
## 6. GENOME ASSEMBLY (Flye)
##############################################################################

echo "Step 5: Assembly"

for fq in "${FILTER_DIR}"/*_filt.fastq; do
    sample=$(basename "${fq}" _filt.fastq)
    outdir="${ASSEMBLY_DIR}/${sample}"

    flye \
        --nano-raw "${fq}" \
        --genome-size "${GENOME_SIZE}" \
        --threads "${THREADS}" \
        --out-dir "${outdir}"
done


##############################################################################
## 7. CHECKM QUALITY ASSESSMENT
##############################################################################

echo "Step 6: CheckM"

for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")
    outdir="${CHECKM_DIR}/${sample}"

    mkdir -p "${outdir}"

    checkm lineage_wf \
        -x fasta \
        -t "${THREADS}" \
        "$(dirname "${assembly}")" \
        "${outdir}"

    checkm qa \
        "${outdir}/lineage.ms" \
        "${outdir}" \
        -f "${outdir}/checkm_summary.tsv"
done


##############################################################################
## 8. OPTIONAL BUSCO
##############################################################################

# Uncomment if needed
# for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
#     sample=$(basename "$(dirname "${assembly}")")
#     busco \
#         -i "${assembly}" \
#         -l bacteria_odb10 \
#         -o "${sample}" \
#         -m genome \
#         --out_path "${BUSCO_DIR}" \
#         -c "${THREADS}" \
#         -f
# done


##############################################################################
## 9. PROKKA ANNOTATION
##############################################################################

echo "Step 7: Prokka annotation"

for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")
    outdir="${PROKKA_DIR}/${sample}"

    prokka \
        --cpus "${THREADS}" \
        --outdir "${outdir}" \
        --prefix "${sample}" \
        "${assembly}"
done


##############################################################################
## 10. GTDB-Tk CLASSIFICATION
##############################################################################

echo "Step 8: GTDB-Tk"

export GTDBTK_DATA_PATH="$HOME/gtdb_database"

for assembly in "${ASSEMBLY_DIR}"/*/assembly.fasta; do
    sample=$(basename "$(dirname "${assembly}")")

    cp -n "${assembly}" "${GTDB_DIR}/genomes/${sample}.fa"
done


gtdbtk classify_wf \
    --genome_dir "${GTDB_DIR}/genomes" \
    --out_dir "${GTDB_DIR}/results" \
    --extension fa \
    --cpus "${THREADS}" \
    --skip_ani_screen \
    --force


##############################################################################
## PIPELINE COMPLETE
##############################################################################

echo "Pipeline completed successfully"