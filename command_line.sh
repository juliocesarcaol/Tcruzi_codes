# Create and make the shell script executable
code trypanosoma.sh                    # Open the script in VS Code (or your default code editor)
chmod +x trypanosoma.sh               # Make the script executable

# Basic checks
ls                                     # List files in the current directory
ls -lhs                                # List files with sizes and human-readable format
clear                                  # Clear the terminal for a clean view

# Run the script to rename FASTA headers and generate the ID map
./trypanosoma.sh

# (Optional) Manually create empty files if needed (only if script fails)
touch trypanosoma_18S_short.fasta
touch trypanosoma_18S_id_map.csv

# Run the script again after ensuring input/output files exist
./trypanosoma.sh

# Check that makeblastdb is installed
makeblastdb

# Create a BLAST nucleotide database from the shortened FASTA file
makeblastdb -in trypanosoma_18S_short.fasta -dbtype nucl -out trypanosoma_18S_db

# Run BLASTn with multiple threads and custom output format (format 6 = tabular)
blastn -query all.fasta -db trypanosoma_18S_db -out blast_results.txt -outfmt 6 -num_threads 4

# Overwrite previous results with a new output filename for clarity
blastn -query all.fasta -db trypanosoma_18S_db -out blast_all.tsv -outfmt 6 -num_threads 4

# Split BLAST output into one file per query
mkdir -p blast_por_query              # Create a directory to store individual query results
awk '{ print > "blast_por_query/"$1".blast" }' blast_all.tsv

# Alternative version that closes each file to avoid open file limit issues
awk '{filename = "blast_por_query/" $1 ".blast"; print >> filename; close(filename)}' blast_all.tsv

# Extract the best hit per query based on the highest bit score (column 12)
awk '\
BEGIN { FS="\t"; OFS="\t" }\
{\
  query=$1;\
  score=$12;\
  if (score > max_score[query]) {\
    max_score[query] = score;\
    max_line[query] = $0;\
  }\
}\
END {\
  for (q in max_line) print max_line[q];\
}\
' blast_all.tsv > blast_best_hits.tsv

# Add a header line to the best hits file for easier parsing in downstream analysis
(echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"; \
 cat blast_best_hits.tsv) > blast_best_hits_with_header.tsv

# Extract all unique subject sequence IDs (column 2) from the BLAST results
cut -f2 blast_all.tsv | sort | uniq > hit_ids_all.txt

# (Optional) Install seqtk if it's not already available
brew install seqtk                   # macOS users can use Homebrew to install seqtk

# Extract matching sequences from the original FASTA using the list of IDs
seqtk subseq trypanosoma_18S_short.fasta hit_ids_all.txt > matched_all_seqs.fasta
