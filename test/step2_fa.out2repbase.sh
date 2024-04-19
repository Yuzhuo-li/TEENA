#!/bin/sh
# The prerequisite for the successful execution of this command is the installation of bedops.
gunzip *.gz && mv ./hg38.fa.out hg38.repbase.txt;

for i in *.repbase.txt; do
    j=${i/.txt/}
    if [ -e ${j}.bed ]; then
        echo "Skip: $i"
    else
        echo "Processing $i ..."
 # To prevent the existence of bugs, it is recommended to use anabsolute path of convert2bed (for example: /home/download//bedops/bin/convert2bed)       
        convert2bed <$i >${j}.full.bed
        cat ${j}.full.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 ":" $4 "\t"
$5 "\t" $6}' >${j}.bed
    fi
done

rm -rf ./hg38.repbase.full.bed; 
rm -rf ./hg38.repbase.txt
