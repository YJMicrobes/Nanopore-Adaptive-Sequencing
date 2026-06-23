# 3_contigs — Contig-Level Assembly Annotation & Taxonomic Profiling

This directory contains contigs-level downstream analyses from metagenome assembly, including functional annotation (EggNOG), taxonomic classification (Kraken2), and summary statistics across multiple taxonomic resolutions.

---

## Overview

This stage focuses on assembled contigs and includes:

- Functional annotation using EggNOG
- Taxonomic classification using Kraken2
- Multi-level contig summary statistics
- Functional gene (KO) distribution and comparison across samples

---

## Files

### Contig Summary

- `Contigs_summary.xlsx`  
  Main summary file containing multiple sheets of contig-level analyses:

  - **Table S2 Contigs** — general contigs statistics and summary table  
  - **Contigs_kingdom** — contigs taxonomic classification at kingdom level  
  - **Contigs_genus** — contigs genus-level taxonomic profiles  
  - **Contigs_KO_venn** — shared/unique KEGG Ortholog (KO) profiles across samples  
  - **Contigs_KO_heatmap** — KO abundance heatmap data for visualization  
  - **assembly_file_size** — assembly statistics including file sizes and contig metrics  

---

### Functional Annotation

- `EggNOG_functions.xlsx`  
  Functional annotation of contigs using EggNOG database, including gene ontology, KEGG pathways, and orthologous group assignments.

---

### Kraken2 Taxonomic Classification

Directory: `contigs_kraken/`

Contains contig-level taxonomic assignments generated using Kraken2 for multiple sample groups.

#### Files:

- `g_kraken2_output.txt` / `g_kraken2_report.txt`  
  Taxonomic classification for group G samples

- `h_kraken2_output.txt` / `h_kraken2_report.txt`  
  Taxonomic classification for group H samples

- `t_kraken2_output.txt` / `t_kraken2_report.txt`  
  Taxonomic classification for group T samples

**File types:**
- `*_output.txt` → raw Kraken2 contig assignments  
- `*_report.txt` → summarized taxonomic abundance profiles

---

### Scripts

- `step6.1.2_kraken_summary.py`  
  Python script used to parse Kraken2 outputs and generate summarized tables.

---

## Workflow Context

This module follows metagenome assembly and serves to:

1. Annotate contigs functionally (EggNOG)
2. Assign taxonomic labels using Kraken2
3. Compare functional gene (KO) distributions across samples
4. Summarize assembly-level statistics
5. Support downstream MAG binning and ecological interpretation

---

## Key Outputs

- **Taxonomy:** Kingdom → Genus-level contig classification  
- **Function:** KEGG Ortholog (KO) profiles and heatmaps  
- **Comparisons:** Venn-based KO overlap across samples  
- **Assembly metrics:** Contig size and assembly summary statistics  

---



