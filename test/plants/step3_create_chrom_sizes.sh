#!/bin/bash
# The prerequisite for the successful execution of this command is the installation of samtools.

# Usage: bash step3_create_chrom_sizes.sh Oryza_sativa.IRGSP-1.0.dna.toplevel.fa

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi

fasta_file=$1

if [ ! -f "$fasta_file" ]; then
    echo "Error: File '$fasta_file' not found!"
    exit 1
fi

# Samtools faidx creates an index for fasta.
samtools faidx "$fasta_file"

# Remove the .fa suffix with basename and generate a chromosome length file
base_name=$(basename "$fasta_file" .fa)
chrom_sizes_file="${base_name}.chrom.sizes"

cut -f1,2 "${fasta_file}.fai" > "$chrom_sizes_file"

echo "Generated chrom.sizes file: $chrom_sizes_file"

cut -f1,2 "${fasta_file}.fai" | awk '{print "chr"$1"\t"$2}' > "$chrom_sizes_file"

echo "Generated chrom.sizes file: $chrom_sizes_file"
