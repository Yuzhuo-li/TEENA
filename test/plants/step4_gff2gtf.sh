#!/bin/bash
# The prerequisite for the successful execution of this command is the installation of gffread.

for gff_file in *.gff3; do
    
    gtf_file="${gff_file%.gff}.gtf"
    
    gffread "$gff_file" -T -o "$gtf_file"
done

# Rename and remove
mv Oryza_sativa.IRGSP-1.0.dna.toplevel.chrom.sizes irgsp1.chrom.sizes
mv irgsp1.TE.repbase.bed irgsp1.repbase.bed
mv Oryza_sativa.IRGSP-1.0.58.gff3.gtf Oryza_sativa.IRGSP-1.0.58.gtf
rm -rf *.gff *TE.bed *.gff3 
