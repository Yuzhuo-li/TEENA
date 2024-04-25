#!/bin/bash
# Renzhe Lyu, 2024/4/22
# This is a script to process TE annotation file of plants
# If there are two strings separated by "/" in the fourth column, but do not have ":", add ":NA"

# Usage: bash step2-2_makeTEfile_attention.sh irgsp1.TE.bed

# Check input bed file
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bed_file>"
    exit 1
fi

bed_file=$1

if [ ! -f "$bed_file" ]; then
    echo "Error: File '$bed_file' not found!"
    exit 1
fi

# Create a new name of output bed file
output_file="${bed_file%.bed}.repbase.bed"

# process the bed file
awk 'BEGIN {OFS="\t"} {
    if (split($4, arr, ":") == 1 && split($4, arr, "/") == 2) {
        $4 = $4 ":NA"  # If there are two strings separated by "/" in the fourth column, but do not have ":", add ":NA"
    }
    print $0
}' "$bed_file" > "$output_file"

echo "Modified BED file created: $output_file"
