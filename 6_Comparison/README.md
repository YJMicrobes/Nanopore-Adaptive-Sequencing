# 6_Comparison — Comparative Genomics & Genome Alignment Analysis

This directory contains comparative genomics analyses of MAGs and isolate genomes, including genome alignment, synteny visualization, gene presence/absence, functional comparison, and ANI-based similarity assessment.

The focus is on **Bacillus** and **Comamonas** genomes derived from MAGs and isolate-based experiments.

---

## Overview

This module integrates multiple comparative approaches:

- Whole-genome alignment (MUMmer / nucmer / XMFA)
- Circos-based genome visualization
- Functional comparison using EggNOG annotations
- Gene presence/absence and KO-level comparisons
- ANI-based genome similarity (skani)
- Read recruitment and abundance mapping
- Structural comparison of genomic organization

---

## Files

### Comparative Summary

- `Comparison_summary.xlsx`  
  Main summary workbook containing integrated comparative genomic results:

  - Unique and core KOs across genomes  
  - Comamonas gene-level and structural comparisons  
  - Bacillus gene-level and structural comparisons  
  - ANI-based similarity matrix (skani_ANI)  

---

### Functional Comparison (EggNOG-based)

- `compare_eggnog_all.xlsx`  
  Functional annotation comparison across genomes, including:

  - *Comamonas_g2_b16*  
  - *Comamonas_h3_b2*  
  - *Comamonas_t5_g20*  
  - *Comamonas_mag7*  
  - *Bacillus_h1_b5*  
  - *Bacillus_t1_b19*  
  - *Bacillus_t2_b23*  
  - *Bacillus_t3_b24*  
  - *Bacillus_mag5*  

Used to compare:
- gene function distributions  
- pathway-level differences  
- orthologous group variation across strains and MAGs  

---

### Genome Alignment & Circos Visualization

#### Bacillus Comparison
Files:
- `Bacilus_nucmer2circos.svg`

Contains Circos-based visualization of genome alignment and structural variation among Bacillus genomes.

---

#### Comamonas Comparison
Files:
- `Comamonas_nucmer2circos.svg`

Contains Circos-based visualization of genome alignment and synteny conservation among Comamonas genomes.

---

### Alignment & Visualization Configuration

- `circos.conf` — Circos configuration file  
- `karyotype.txt` — genome structure definition for Circos plots  
- `links.txt` / `links.filtered.txt` — alignment links used for visualization  

---

### Gene Presence / Functional Analysis Scripts

- `step15.1_xmfa_to_links.py`  
  Converts XMFA alignments into Circos-compatible link format  

- `step15.1_xmfa_to_links_v2.py`  
  Updated XMFA-to-Circos conversion pipeline  

- `step15.2_gene_presence.py`  
  Computes gene presence/absence across genomes  

- `step15.3_gene_presence_heatmap.py`  
  Generates heatmap visualization of gene presence/absence patterns  

---

### Phylogenetic Tree

- `tree_Isolates_MAG.nwk`  
  Phylogenetic tree combining isolate genomes and MAGs for evolutionary comparison  

---

## Removed / Deprecated Content

The following previously used outputs were removed during cleanup to reduce repository size and eliminate redundant or intermediate files:

- Large XMFA alignment files  
- Intermediate Circos temporary directories  
- Raw `.delta`, `.coords`, and alignment caches  
- Additional BRIG/HGT visualization outputs  

This cleanup improves repository portability and ensures only analysis-ready outputs are retained.

---

## Workflow Context

This stage integrates outputs from earlier pipeline steps:

1. MAG reconstruction and isolate genome assembly (4_MAGs, 5_Isolates)  
2. Functional annotation (EggNOG)  
3. Genome alignment (MUMmer / nucmer / XMFA)  
4. Structural comparison (Circos visualization)  
5. Gene presence/absence analysis  
6. ANI-based similarity estimation (skani)  
7. Read recruitment and abundance validation  

---

## Notes

- All analyses are reproducible using scripts in `1_pipelines/`
- Heavy intermediate files were removed to maintain repository efficiency
- This directory retains only **final or analysis-ready outputs**