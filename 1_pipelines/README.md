# 1_pipelines — Nanopore Metagenomics & Genome-Resolved Analysis Pipeline

This directory contains a modular, stepwise workflow for processing Oxford Nanopore sequencing data, from raw reads to genome-resolved analyses, including assembly, binning, quality control, taxonomic classification, genome comparison, and pangenome analysis.

Each step is implemented as an independent shell script to ensure reproducibility and flexibility on HPC environments.

---

## Overview of Workflow

The pipeline is organized into sequential stages:

1. Data upload & preprocessing  
2. Nanopore basecalling and demultiplexing  
3. Read quality control  
4. Metagenome assembly  
5. Genome binning (multiple methods)  
6. Bin refinement and integration  
7. Quality assessment (CheckM)  
8. Taxonomic classification (GTDB-Tk)  
9. Genome comparison & pangenome analysis  
10. Isolate / whole-genome Nanopore pipeline (extended workflow)

---

## Pipeline Steps

### 1. Data Upload
- `step1_upload_data.sh`  
Uploads raw sequencing data to the computing environment or HPC storage system.

---

### 2. Basecalling
- `step2_basecalling.sh`  
Performs basecalling of raw Nanopore signal data to generate FASTQ files.

---

### 3. Demultiplexing
- `step3_demultiplex.sh`  
Separates reads by barcode/sample identifiers.

---

### 4. Quality Assessment
- `step4_nanoplot_quality_check.sh`  
Generates sequencing quality summaries using NanoPlot.

---

### 5. Read Filtering & Subsampling
- `step5_nanoFilt_quality_control.sh`  
Filters reads based on quality and length thresholds.  
- `step5_sub_reads_based_analysis.sh`  
Optional downstream analysis on subsetted reads.

---

### 6. Assembly
- `step6_flye_assembly.sh`  
Performs long-read assembly using Flye.  
- `step6_sub_contigs_analyses.sh`  
Performs downstream contig-level analysis and summaries.

---

### 7. Genome Binning (Multiple Methods)
- `step7_maxbin.sh` — MaxBin2 binning  
- `step8_metabat.sh` — MetaBAT2 binning  
- `step9_metacoag.sh` — MetaCoAG binning  

---

### 8. Bin Refinement & Integration
- `step10_dastool.sh`  
Integrates multiple binning results using DAS Tool to generate high-quality consensus bins.

---

### 9. Genome Quality Assessment
- `step11_checkM.sh`  
Evaluates completeness and contamination of MAGs using CheckM.

---

### 10. Taxonomic Classification
- `step12_GTDB.sh`  
Assigns taxonomy to MAGs using GTDB-Tk.

---

### 11. Nanopore Isolate / Whole Genome Pipeline
- `step13_nanopore_isolate_whole_genome_pipeline.sh`  
End-to-end pipeline for isolate genomes or high-quality single-genome assemblies.

---

### 12. Genome Comparison
- `step14_genome_comparison.sh`  
Performs pairwise or groupwise genome comparisons (e.g., ANI, alignment-based metrics).

---

### 13. Pangenome Analysis
- `step15_Pangenome_comparison.sh`  
Constructs and analyzes pangenomes across multiple genomes or MAGs.

---

## Output Structure (Typical)

results/
├── raw_reads/
├── basecalled/
├── demultiplexed/
├── qc_reports/
├── filtered_reads/
├── assemblies/
├── bins/
│ ├── maxbin/
│ ├── metabat/
│ ├── metacoag/
│ └── dastool/
├── checkm/
├── gtdbtk/
├── genome_comparison/
└── pangenome/


---

## Dependencies (Typical)

- Nanopore Tools: Guppy / Dorado
- Quality Control: NanoPlot, NanoFilt
- Assembly: Flye
- Binning: MaxBin2, MetaBAT2, MetaCoAG
- Bin refinement: DAS Tool
- Quality assessment: CheckM
- Taxonomy: GTDB-Tk
- Comparative genomics: ANI tools, MUMmer / FastANI
- Pangenome: Roary / Panaroo / Anvi’o (depending on implementation)

---

## Notes

- Each script is designed for HPC execution (SLURM-compatible).
- Steps can be run independently if intermediate outputs already exist.
- Recommended to maintain consistent directory structure for reproducibility.
- Ensure database paths (e.g., GTDB-Tk, CheckM) are correctly configured in each script.

---

## Author
Jing Yuan, Ph.D.

---