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

## Script: `download_trypanosoma_18S.R`

This script performs a batch download of 18S rRNA sequences (SSU) for the genus *Trypanosoma* from the NCBI Nucleotide database using the `rentrez` R package.

### Features

- Retrieves up to 10,000 sequences matching a specific query (length: 100–5000 bp).
- Exports sequences in FASTA format.
- Extracts metadata (accession, organism name, sequence definition, and length) from GenBank records.
- Saves metadata to a CSV file.

### Usage

Make sure you have R and the required packages (`rentrez`, `tibble`, `stringr`) installed.

Run the script from R or RStudio:

```r
source("download_trypanosoma_18S.R")
```

## Script: `rename_fasta_headers.sh`

This Bash script processes a FASTA file by shortening long header lines and replacing them with unique, simplified identifiers. It also generates a mapping file to preserve the link between the original and shortened IDs.

### Features

- Truncates complex FASTA headers into short IDs (max 50 characters)
- Ensures uniqueness with an internal counter
- Generates a CSV map linking `original_id` → `short_id`

### Usage

```bash
chmod +x rename_fasta_headers.sh
./rename_fasta_headers.sh
```


