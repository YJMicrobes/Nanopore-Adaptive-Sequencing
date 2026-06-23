# 4_MAGs — Metagenome-Assembled Genomes (MAGs) Analysis

This directory contains downstream analyses of metagenome-assembled genomes (MAGs), including quality summaries, read recruitment, ANI comparisons, and phylogenetic placement.

---

## Overview

This module represents the genome-resolved stage of the metagenomics pipeline, focusing on:

- MAG quality and summary statistics  
- Read recruitment to MAGs  
- Genome-wide similarity (ANI) analysis  
- Phylogenetic tree placement  
- Horizontal gene transfer / bin-level comparisons  

---

## Files

### MAG Summary

- `MAGs_summary.xlsx`  
  Main summary workbook containing multiple MAG-level analysis tables:

  - **Table S3 Bins** — general MAG statistics (completeness, contamination, size, etc.)  
  - **MAGs_all** — full list of recovered MAGs and metadata  
  - **MAGs_reads** — read recruitment / coverage information per MAG  
  - **MAGs_ANI** — Average Nucleotide Identity (ANI) comparisons between MAGs and isolates 
  - **Bins_colony_HGT_all** — full bin list for all the colony test and HGT sequencing  

---

### Phylogenetic Tree

- `bac120_user_tree.nwk`  
  Phylogenetic tree of bacterial MAGs constructed using the GTDB-Tk bacterial marker set (bac120).  
  This tree places recovered MAGs in a reference genome phylogeny for taxonomic context.



---