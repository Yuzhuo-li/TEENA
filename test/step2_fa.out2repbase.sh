#!/bin/sh
# The prerequisite for the successful execution of this command is the installation of bedops.
gunzip *.gz && mv ./hg38.fa.out hg38.repbase.txt;

# To prevent the existence of bugs, it is recommended to use anabsolute path of rmsk2bed (for example: /home/download/bedops/bin/rmsk2bed)       
rmsk2bed <hg38.repbase.txt >hg38.repbase.full.bed
cat hg38.repbase.full.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 ":" $4 "\t" $5 "\t" $6}' >hg38.repbase.bed
   

rm -rf ./hg38.repbase.full.bed; 
rm -rf ./hg38.repbase.txt
