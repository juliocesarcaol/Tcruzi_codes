# Install and load required packages
install.packages("rentrez")  # Only run once
library(rentrez)
library(tibble)
library(stringr)

# Set your email (required by NCBI Entrez API)
Sys.setenv(ENTREZ_EMAIL = "jceszar@gmail.com")

# Define the search query: 18S ribosomal RNA (SSU) in Trypanosoma, between 100–5000 bp in length
query <- "Trypanosoma[Organism] AND (18S ribosomal RNA[Title] OR SSU[Title]) AND 100:5000[SLEN]"

# Perform the search on NCBI Nucleotide database
search_results <- entrez_search(db = "nucleotide", term = query, retmax = 10000)
cat("Total sequences found:", search_results$count, "\n")

# Extract the list of unique IDs (UIDs)
uids <- search_results$ids
batch_size <- 200  # Batch size for downloading sequences

# Output files
fasta_file <- "trypanosoma_18S.fasta"
csv_file <- "trypanosoma_18S_metadata.csv"
file.create(fasta_file)  # Create the FASTA output file

# Initialize list to store metadata
metadata_list <- list()

# Loop through UIDs in batches to download data
for (start in seq(1, length(uids), by = batch_size)) {
  end <- min(start + batch_size - 1, length(uids))
  cat("Downloading sequences", start, "to", end, "...\n")
  
  # Download sequences in FASTA format
  fasta_data <- entrez_fetch(
    db = "nucleotide",
    id = uids[start:end],
    rettype = "fasta",
    retmode = "text"
  )
  write(fasta_data, file = fasta_file, append = TRUE)
  
  # Download entries in GenBank format to extract metadata
  gb_data <- entrez_fetch(
    db = "nucleotide",
    id = uids[start:end],
    rettype = "gb",
    retmode = "text"
  )
  
  # Split GenBank file into individual entries
  entries <- strsplit(gb_data, "//\n")[[1]]
  
  # Extract metadata from each entry
  for (entry in entries) {
    acc <- str_match(entry, "ACCESSION\\s+([A-Z0-9_]+)")[,2]
    organism <- str_match(entry, "  ORGANISM\\s+(.+)")[,2]
    definition <- str_match(entry, "DEFINITION\\s+(.+)")[,2]
    length <- str_match(entry, "LOCUS\\s+\\S+\\s+(\\d+)")[,2]
    
    if (!is.na(acc)) {
      metadata_list[[length(metadata_list) + 1]] <- tibble(
        accession = acc,
        organism = organism,
        definition = definition,
        length = as.integer(length)
      )
    }
  }
}

# Combine metadata into a data frame and export to CSV
metadata_df <- bind_rows(metadata_list)
write.csv(metadata_df, csv_file, row.names = FALSE)

# Final message
cat("Download complete.\nFASTA file:", fasta_file, "\nMetadata CSV:", csv_file, "\n")
