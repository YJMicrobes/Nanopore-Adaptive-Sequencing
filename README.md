# Nanopore Adaptive Sequencing Resolves Bacterial Genomes from the Mosquito Microbiome

This repository contains bioinformatic pipelines and analysis workflows developed for the study:

**"Adaptive Sequencing Resolves Bacterial Genomes from the Mosquito Microbiome: Linking Function to Larval Development"**

This project applies Oxford Nanopore long-read adaptive sequencing approaches to resolve bacterial genomes from complex mosquito microbiome communities. The workflow integrates raw read processing, genome assembly, metagenome-assembled genome (MAG) reconstruction, genome quality assessment, taxonomic classification, comparative genomics, and functional analysis.

The complete analysis workflow is organized in the `1_pipelines/` directory, including:

- Nanopore raw signal processing, basecalling, demultiplexing, and read quality control
- Long-read genome assembly using Flye
- MAG recovery using multiple genome binning approaches (MetaBAT2, MaxBin2, and MetaCoAG) followed by bin refinement with DAS Tool
- Genome quality evaluation and taxonomy assignment using CheckM and GTDB-Tk
- Nanopore isolate genome assembly and annotation
- Genome comparison based on ANI, read mapping, phylogenomics, and functional annotation
- Pangenome analysis, gene presence/absence comparison, and whole-genome visualization

Repository structure:

1_pipelines/ # Analysis scripts
2_reads/ # Sequencing reads
3_contigs/ # Genome assemblies and contigs
4_MAGs/ # Metagenome-assembled genomes
5_Isolates/ # Isolate genome analyses
6_Comparison/ # Comparative genomic analyses


The scripts are organized sequentially (`step1`–`step15`) to represent the complete computational workflow from raw Nanopore sequencing data to genome-resolved microbiome analysis.
Please adapt file paths, database locations, and computational resources according to your local environment before running the pipeline.
For questions or collaboration regarding this workflow, please open an issue in this repository.


step1_upload_data.sh
step2_basecalling.sh
step3_demultiplex.sh
step4_nanoplot_quality_check.sh
step5_nanoFilt_quality_control.sh
step5_sub_reads_based_analysis.sh
step6_flye_assembly.sh
step6_sub_contigs_analyses.sh
step7_maxbin.sh
step8_metabat.sh
step9_metacoag.sh
step10_dastool.sh
step11_checkM.sh
step12_GTDB.sh
step13_nanopore_isolate_whole_genome_pipeline.sh
step14_genome_comparison.sh
step15_Pangenome_comparison.sh