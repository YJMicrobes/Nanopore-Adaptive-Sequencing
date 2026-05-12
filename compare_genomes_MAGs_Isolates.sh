#!/bin/bash
set -euo pipefail

##############################################################################
## 0. PROJECT INFO
##############################################################################
# Name: Jing Yuan
# Date: 2026-05-12
# Project: Comparative genomics of MAGs and isolate genomes
# Institute: Connecticut Agricultural Experiment Station

##############################################################################
## 1. INPUTS
##############################################################################
MAG_DIR="MAGs"
ISOLATE_DIR="Isolates"
EXT="fa"

##############################################################################
## 2. OUTPUT STRUCTURE
##############################################################################
OUT="comparative_genomics"

ANI="${OUT}/01_ani"
TREE="${OUT}/02_tree"
PRODIGAL="${OUT}/03_prodigal"
MMSEQS="${OUT}/04_mmseqs"
EGGNOG="${OUT}/05_eggnog"
PROKKA="${OUT}/06_prokka"
PANAROO="${OUT}/07_panaroo"
MAUVE="${OUT}/08_mauve"
CIRCOS="${OUT}/09_circos"

THREADS=30

mkdir -p ${OUT} ${ANI} ${TREE} ${PRODIGAL} ${MMSEQS} ${EGGNOG} ${PROKKA} ${PANAROO} ${MAUVE} ${CIRCOS}

##############################################################################
## 3. COLLECT GENOMES
##############################################################################
echo "STEP 1: Collect genomes"

ALL="${OUT}/all_genomes"
mkdir -p ${ALL}

cp ${MAG_DIR}/*.${EXT} ${ALL}/ 2>/dev/null || true
cp ${ISOLATE_DIR}/*.${EXT} ${ALL}/ 2>/dev/null || true

ls ${ALL}/*.${EXT} > ${OUT}/genomes.txt

##############################################################################
## 4. ANI (FastANI)
##############################################################################
echo "STEP 2: FastANI"

fastANI \
  --ql ${OUT}/genomes.txt \
  --rl ${OUT}/genomes.txt \
  -o ${ANI}/fastani.tsv

##############################################################################
## 5. ANI (skani)
##############################################################################
echo "STEP 3: skani"

skani triangle ${ALL}/*.${EXT} -o ${ANI}/skani.tsv

##############################################################################
## 6. GENE PREDICTION (Prodigal)
##############################################################################
echo "STEP 4: Prodigal"

for g in ${ALL}/*.${EXT}; do
    base=$(basename $g .${EXT})

    prodigal \
        -i $g \
        -a ${PRODIGAL}/${base}.faa \
        -d ${PRODIGAL}/${base}.fna \
        -p single
done

##############################################################################
## 7. PROTEIN CLUSTERING (MMseqs2)
##############################################################################
echo "STEP 5: MMseqs clustering"

cat ${PRODIGAL}/*.faa > ${PRODIGAL}/all.faa

mmseqs createdb ${PRODIGAL}/all.faa ${MMSEQS}/db

mmseqs cluster ${MMSEQS}/db ${MMSEQS}/cluster ${MMSEQS}/tmp \
    --min-seq-id 0.85 -c 0.8 --cov-mode 1

mmseqs createtsv \
    ${MMSEQS}/db ${MMSEQS}/db ${MMSEQS}/cluster \
    ${MMSEQS}/clusters.tsv

##############################################################################
## 8. FUNCTIONAL ANNOTATION (eggNOG)
##############################################################################
echo "STEP 6: eggNOG annotation"

emapper.py \
  -i ${PRODIGAL}/all.faa \
  -o eggnog \
  --output_dir ${EGGNOG} \
  --cpu 10 \
  -m diamond

##############################################################################
## 9. PROKKA + PAN-GENOME (Panaroo)
##############################################################################
echo "STEP 7: Prokka + Panaroo"

for g in ${ALL}/*.${EXT}; do
    base=$(basename $g .${EXT})

    prokka \
      --cpus ${THREADS} \
      --outdir ${PROKKA}/${base} \
      --prefix ${base} \
      $g
done

find ${PROKKA} -name "*.gff" > ${PANAROO}/gff.txt

panaroo \
  -i $(cat ${PANAROO}/gff.txt) \
  -o ${PANAROO} \
  -t ${THREADS} \
  -c 0.9 \
  --clean-mode moderate

##############################################################################
## 10. GTDB-Tk TREE - FastTree
##############################################################################
echo "STEP 8: GTDB-Tk tree (optional)"

# Requires GTDB-Tk already run separately

if [ -f ${TREE}/gtdbtk.bac120.user_msa.fasta ]; then
    FastTree -nt \
      ${TREE}/gtdbtk.bac120.user_msa.fasta \
      > ${TREE}/gtdb_tree.nwk
fi

##############################################################################
## 11. MULTIPLE GENOME ALIGNMENT (Mauve)
##############################################################################
echo "STEP 9: Mauve alignment"

if [ $(ls ${ALL}/*.${EXT} | wc -l) -ge 2 ]; then
    progressiveMauve \
      --output=${MAUVE}/genomes.xmfa \
      ${ALL}/*.${EXT}
fi

##############################################################################
## 12. XMFA → CIRCOS LINKS
##############################################################################
echo "STEP 10: Circos links"

if [ -f ${MAUVE}/genomes.xmfa ]; then
    python scripts/xmfa_to_links.py \
      ${MAUVE}/genomes.xmfa \
      ${CIRCOS}/links.txt
fi

##############################################################################
## 13. DONE
##############################################################################
echo "PIPELINE COMPLETE"
echo "Results in: ${OUT}"