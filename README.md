# Tcruzi_codes
# Trypanosoma cruzi Genetic Data Workflow

This repository provides tools and documentation for downloading, processing, and searching *Trypanosoma cruzi* genomic sequences using NCBI's Entrez system (`rentrez` package in R) and local BLAST+ databases created with `makeblastdb`.

## Overview

This workflow is designed for bioinformatics analysis of *Trypanosoma cruzi*, the etiological agent of Chagas disease. It allows researchers to:

- Search and download *T. cruzi* nucleotide or protein sequences from NCBI using R.
- Build custom BLAST databases from downloaded FASTA files.
- Perform local BLAST searches for comparative analysis or annotation.

---

## Requirements

### Software

- [R](https://cran.r-project.org/)
  - `rentrez` R package (`install.packages("rentrez")`)
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - Includes `makeblastdb`, `blastn`, `blastp`, etc.
- Python 3 (optional, for processing or integration scripts)



