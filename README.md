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
# Trypanosoma 18S BLAST Pipeline (Shell Script Version)

This repository provides a complete workflow for downloading, processing, and searching 18S rRNA sequences from *Trypanosoma* using BLAST. All processing is done through R and Bash scripts, without requiring a Makefile.

## Overview

The pipeline performs:

1. **Download Trypanosoma 18S sequences** from NCBI using R and the `rentrez` package.
2. **Format FASTA headers** to standardized short names using a Bash script (`trypanosoma.sh`).
3. **Create a BLAST database** with `makeblastdb`.
4. **Run BLASTn** with `blastn` and parse results.
5. **Extract best hits and retrieve matching sequences** using `awk` and `seqtk`.

## Requirements

* R with:

  * `rentrez`
  * `tibble`
  * `stringr`
* Bash shell
* NCBI BLAST+ tools (`makeblastdb`, `blastn`)
* `awk`
* `seqtk`

## Usage

### Step 1: Download Trypanosoma 18S sequences

Run the R script to retrieve sequences from NCBI:

```R
source("download_trypanosoma.R")
```

This script will generate two files:

* `trypanosoma_18S.fasta`
* `trypanosoma_18S_metadata.csv`

### Step 2: Format FASTA headers

Run the Bash script to clean and shorten FASTA headers:

```bash
chmod +x rename_fasta_headers.sh
./rename_fasta_headers.sh
```

Output files:

* `trypanosoma_18S_short.fasta`
* `trypanosoma_18S_id_map.csv`

### Step 3: Create a BLAST database

```bash
makeblastdb -in trypanosoma_18S_short.fasta -dbtype nucl -out trypanosoma_18S_db
```

### Step 4: Run BLAST search

Ensure your query file is named `all.fasta`, then run:

```bash
blastn -query all.fasta -db trypanosoma_18S_db -out blast_all.tsv -outfmt 6 -num_threads 4
```

### Step 5: Process BLAST output

Split output per query:

```bash
mkdir -p blast_por_query
awk '{filename = "blast_por_query/" $1 ".blast"; print >> filename; close(filename)}' blast_all.tsv
```

Extract best hits per query:

```bash
awk '
BEGIN { FS="\t"; OFS="\t" }
{
  query=$1;
  score=$12;
  if (score > max_score[query]) {
    max_score[query] = score;
    max_line[query] = $0;
  }
}
END {
  for (q in max_line) print max_line[q];
}' blast_all.tsv > blast_best_hits.tsv
```

Add header:

```bash
(echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"; \
 cat blast_best_hits.tsv) > blast_best_hits_with_header.tsv
```

### Step 6: Extract matched sequences

```bash
cut -f2 blast_all.tsv | sort | uniq > hit_ids_all.txt
seqtk subseq trypanosoma_18S_short.fasta hit_ids_all.txt > matched_all_seqs.fasta
```

### Optional: Install `seqtk`

```bash
# On macOS with Homebrew:
brew install seqtk
```

## Output Files

* `blast_all.tsv`: full BLAST output
* `blast_best_hits_with_header.tsv`: best-scoring hit per query
* `matched_all_seqs.fasta`: extracted subject sequences from BLAST

## License

MIT License

## Author

\[Julio C. Carrion-Olmedo] – \[INABIO / juliocesarcaol]




