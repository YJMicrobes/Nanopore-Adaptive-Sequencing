# 🧬 Nanopore Adaptive Sequencing Pipelines

A modular suite of Oxford Nanopore long-read bioinformatics pipelines for:

- 🦟 Microbiome metagenomics (mosquito-associated systems)
- 🧬 Single-isolate whole genome assembly
- 🔬 MAG vs isolate comparative genomics

Developed for HPC environments with SLURM/module support.

------------------
# 📦 Repository Overview

This repository contains three main pipelines:

- 🦟 Metagenomics	Mosquito microbiome analysis: nanopore_adaptive_mosquito_microbiome.sh
- 🧬 Isolates	Whole genome assembly + annotation: nanopore_isolate_whole_genome_pipeline.sh
- 🔬 Comparative genomics	MAG vs isolate comparison: compare_genomes_MAGs_Isolates.sh

------------------
## 🦟 1. Nanopore Adaptive Microbiome Pipeline
### 🧭 Purpose
Metagenomic analysis of Nanopore sequencing data from mosquito-associated microbiomes.

### ⚙️ Workflow
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

### 📁 Script
nanopore_adaptive_mosquito_microbiome.sh

How to run:
```
sbatch nanopore_adaptive_mosquito_microbiome.sh
```

------------------
## 🧬 2. Nanopore Isolate Genome Pipeline
### 🧭 Purpose
High-quality bacterial genome assembly and annotation from isolate Nanopore sequencing.

### ⚙️ Workflow
Basecalling & QC
Genome assembly (Flye or equivalent)
Gene prediction (optional depending on setup)
Functional annotation (Prokka, eggNOG-mapper)
Basic genome statistics

### 📁 Script
nanopore_isolate_whole_genome_pipeline.sh

How to run:
```
sbatch nanopore_isolate_whole_genome_pipeline.sh
```

------------------
## 🔬 3. Comparative Genomics (MAGs vs Isolates)
### 🧭 Purpose
Genome-scale comparison between metagenome-assembled genomes (MAGs) and cultured isolates.

### ⚙️ Workflow
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

### 📁 Script
compare_genomes_MAGs_Isolates.sh

How to run:
```
sbatch compare_genomes_MAGs_Isolates.sh
```

------------------
## ⚙️ Software Requirements
Designed for HPC systems (SLURM/module or conda).

## Core tools
- Dorado
- Flye
- FastANI / skani
- Kraken2 / Kaiju
- Prokka
- eggNOG-mapper
- Bowtie2 / Samtools
- MetaBAT2 / MaxBin2 / MetaCoAG
- DAS Tool
- CheckM
- GTDB-Tk
- MMseqs2
- Panaroo
- FastTree
- Mauve / MUMmer
- Circos

## 🧠 Input Data
- Nanopore POD5 / FASTQ files
- MAG assemblies (FASTA)
- Isolate genomes (FASTA)
- Host genome (e.g., mosquito reference)


## Module system dependencies

This pipeline assumes an HPC environment using the module system.

Before running, always check available software versions:

```
module avail
```

Then load appropriate modules:

```
module load Dorado
module load samtools
module load cutadapt
```
⚠️ Module names and versions may differ between systems.



### 3. No module system available?

If your system does not support module load, you can create a Conda environment instead:

```
conda create -n nanopore_env dorado samtools cutadapt seqkit mothur -c bioconda -c conda-forge
conda activate nanopore_env
```

Or install tools manually if required.



------------------
📜 Citation
If you use this pipeline, please cite our paper:

> TBD (under review) 

---

For the Nanopore metagenomic sequencing workflow, please see:

> Yuan, J., LaReau, J., Lawrence, B., Meadows-McDonnell, M., Steven, B., & Shabtai, I. (2026). Thresholds within a soil moisture gradient drive abrupt transitions in microbial community structure resulting in distinct carbon utilization patterns. Soil Ecology Letters. https://doi.org/10.1016/j.apsoil.2026.106966



If you need additional:

📊 Python scripts (plots, ANI heatmaps, trees)
📈 R scripts (statistical analysis, visualization)
🧬 Custom comparative genomics analyses

feel free to contact me — I’m happy to share helper scripts and templates.

👤 Author
> Jing Yuan, Ph.D. \
> Microbial Ecology | Metagenomics | Microbiome Bioinformatics
