#!/bin/sh
# This is the code to run teena with default parameters.
# It is recommended to have a version of Python above 3.8 to facilitate drawing.
# GATA3_hg38.bed is only a test file for code running, it is recommended to use a bed file that matches the species genome for analysis.
python ../../teena.py -q GATA3_hg38.bed -d irgsp1.repbase.bed -fa Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -ch irgsp1.chrom.sizes -a Oryza_sativa.IRGSP-1.0.58.gtf -o rice
