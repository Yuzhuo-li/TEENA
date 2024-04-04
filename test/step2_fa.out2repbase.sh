#!/bin/sh
gunzip *.gz && mv ./hg38.fa.out hg38.repbase.txt;

for i in *.repbase.txt; do
    j=${i/.txt/}
    if [ -e ${j}.bed ]; then
        echo "Skip: $i"
    else
        echo "Processing $i ..."
        rmsk2bed <$i >${j}.full.bed
        cat ${j}.full.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $11 ":" $4 "\t"
$5 "\t" $6}' >${j}.bed
    fi
done

rm -rf ./hg38.repbase.full.bed; 
rm -rf ./hg38.repbase.txt