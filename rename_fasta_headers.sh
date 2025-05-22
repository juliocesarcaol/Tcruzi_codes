#!/bin/bash

# -------------------------------
# Script: rename_fasta_headers.sh
# Description:
#   This script reads a FASTA file, renames long headers with shorter unique IDs,
#   and creates a CSV mapping of original IDs to shortened ones.
#   Useful for tools that have issues with long FASTA headers.
# -------------------------------

# Input and output files
INPUT="trypanosoma_18S.fasta"             # Original FASTA file
OUTPUT="trypanosoma_18S_short.fasta"      # Output FASTA with short IDs
MAP="trypanosoma_18S_id_map.csv"          # CSV map of original IDs to new IDs

# Check if the input file exists
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file '$INPUT' not found."
    exit 1
fi

# Initialize output files
> "$OUTPUT"
echo "original_id,short_id" > "$MAP"      # Add header to CSV map

# Counter to ensure ID uniqueness
count=1

# Read the FASTA file line by line
while IFS= read -r line; do
    if [[ "$line" == \>* ]]; then
        # Extract the original FASTA header (remove the '>' character)
        original_id="${line:1}"

        # Remove problematic characters, keep alphanumerics, underscores, and hyphens
        # Then truncate to max 40 characters
        base_id=$(echo "$original_id" | tr -cd '[:alnum:]_-' | cut -c1-40)

        # Create short unique ID using a counter and base_id
        short_id="T18S_${count}_${base_id}"
        short_id=$(echo "$short_id" | cut -c1-50)  # Ensure final ID is max 50 characters

        # Write the new FASTA header and update the mapping CSV
        echo ">$short_id" >> "$OUTPUT"
        echo "${original_id},${short_id}" >> "$MAP"

        ((count++))  # Increment counter
    else
        # If the line is a sequence, write it unchanged
        echo "$line" >> "$OUTPUT"
    fi
done < "$INPUT"

# Final message
echo "Short FASTA written to: $OUTPUT"
echo "ID mapping written to: $MAP"
