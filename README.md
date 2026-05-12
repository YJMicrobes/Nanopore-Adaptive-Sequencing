🧬 Nanopore Microbiome & Comparative Genomics Toolkit

A modular suite of Oxford Nanopore long-read bioinformatics pipelines for:

🦟 Microbiome metagenomics (mosquito-associated systems)
🧬 Single-isolate whole genome assembly
🔬 MAG vs isolate comparative genomics

Developed for HPC environments with SLURM/module support.

📦 Repository Overview

This repository contains three main pipelines:

Pipeline	Purpose	Script
🦟 Metagenomics	Mosquito microbiome analysis	nanopore_adaptive_mosquito_microbiome.sh
🧬 Isolates	Whole genome assembly + annotation	nanopore_isolate_whole_genome_pipeline.sh
🔬 Comparative genomics	MAG vs isolate comparison	compare_genomes_MAGs_Isolates.sh
🦟 1. Mosquito Microbiome Pipeline
🧭 Purpose

Metagenomic analysis of Nanopore sequencing data from mosquito-associated microbiomes.

⚙️ Workflow
Basecalling (Dorado)
Read QC (NanoPlot, NanoFilt)
Read taxonomy (Kaiju)
Metagenome assembly (Flye)
Contig taxonomy (Kraken2)
Functional annotation (Prokka, eggNOG)
Genome binning (MetaBAT2, MaxBin2, MetaCoAG)
Bin refinement (DAS Tool)
Quality assessment (CheckM)
Host filtering (BLAST)
Adaptive sequencing target design
📁 Script
nanopore_adaptive_mosquito_microbiome.sh
🧬 2. Nanopore Isolate Genome Pipeline
🧭 Purpose

High-quality bacterial genome assembly and annotation from isolate Nanopore sequencing.

⚙️ Workflow
Basecalling & QC
Genome assembly (Flye or equivalent)
Gene prediction (optional depending on setup)
Functional annotation (Prokka, eggNOG-mapper)
Basic genome statistics
📁 Script
nanopore_isolate_whole_genome_pipeline.sh
📊 Output
Assembled genomes (FASTA)
Annotated genes (GFF, FAA, FNA)
Functional annotations
🔬 3. Comparative Genomics (MAGs vs Isolates)
🧭 Purpose

Genome-scale comparison between metagenome-assembled genomes (MAGs) and cultured isolates.

⚙️ Workflow
🧬 Genome similarity
FastANI (pairwise ANI)
skani (fast genome distance estimation)
🧪 Gene-level analysis
Prodigal gene prediction
MMseqs2 clustering
eggNOG functional annotation
🧬 Pan-genome analysis
Prokka annotation
Panaroo pan-genome construction
Gene cluster comparison
🌳 Phylogenomics (optional)
GTDB-Tk marker alignment
FastTree phylogenetic tree
🌐 Genome structure comparison
Mauve whole-genome alignment
Circos synteny visualization
📁 Script
compare_genomes_MAGs_Isolates.sh
📊 Outputs
ANI similarity matrices
Gene clusters (core/accessory)
Functional profiles (KO/COG)
Pan-genome structure
Phylogenetic trees
Synteny plots (Circos)
⚙️ Software Requirements

Designed for HPC systems (SLURM/module or conda).

Core tools
Dorado
Flye
FastANI / skani
Kraken2 / Kaiju
Prokka
eggNOG-mapper
Bowtie2 / Samtools
MetaBAT2 / MaxBin2 / MetaCoAG
DAS Tool
CheckM
GTDB-Tk
MMseqs2
Panaroo
FastTree
Mauve / MUMmer
Circos
🧠 Input Data
Nanopore POD5 / FASTQ files
MAG assemblies (FASTA)
Isolate genomes (FASTA)
Optional host genome (e.g., mosquito reference)
🚀 How to Run
🦟 Microbiome pipeline
sbatch nanopore_adaptive_mosquito_microbiome.sh
🧬 Isolate pipeline
sbatch nanopore_isolate_whole_genome_pipeline.sh
🔬 Comparative genomics pipeline
sbatch compare_genomes_MAGs_Isolates.sh
📁 Output Directory Structure
Mosquito pipeline
1_pod5/
2_demultiplex/
3_nanofilt/
4_kaiju_reads/
5_assembly/
6_kraken/
7_prokka/
8_eggnog/
9_binning/
10_dastool/
11_checkm/
12_gtdbtk/
13_host_check/
15_adaptive_targets/
Comparative genomics
comparative_genomics/
├── 01_ani
├── 02_tree
├── 03_prodigal
├── 04_mmseqs
├── 05_eggnog
├── 06_prokka
├── 07_panaroo
├── 08_mauve
└── 09_circos
📊 Applications
Mosquito microbiome profiling
Environmental metagenomics
Genome-resolved microbiology
MAG vs isolate evolutionary comparison
Horizontal gene transfer detection
Adaptive sequencing design (ONT)
⚠️ Notes
Designed for HPC environments
Requires prebuilt databases:
Kraken2
GTDB-Tk
Kaiju
Large intermediate files expected (BAM, assemblies)
BLAST remote search may be slow
📜 Citation

If you use this repository, please cite:

Flye
Kraken2
GTDB-Tk
MetaBAT2
DAS Tool
eggNOG-mapper
Prokka
Panaroo
🧩 Support

If you need additional:

📊 Python scripts (plots, ANI heatmaps, trees)
📈 R scripts (statistical analysis, visualization)
🧬 Custom comparative genomics analyses

feel free to contact me — I’m happy to share helper scripts and templates.

👤 Author

Jing Yuan, Ph.D.
Microbial Ecology | Metagenomics | Microbiome Bioinformatics
