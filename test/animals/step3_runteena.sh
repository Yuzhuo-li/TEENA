#!/bin/sh
# This is the code to run teena with default parameters.
# It is recommended to have a version of Python above 3.8 to facilitate drawing.
python ../teena.py -q GATA3_hg38.bed -d hg38.repbase.bed -fa hg38.fa -ch hg38.chrom.sizes -a Homo_sapiens.GRCh38.110.gtf -o test
