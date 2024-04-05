#!/usr/bin/env python
#-*- coding:utf-8 -*-
# author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# datetime: 2023-11-23
# software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-09: set a function to extract the promoter regions from genome annotation file (like gtf file)

def get_sub_result(fasta_file, sub_area_file, usize: int, dsize: int):
    # print('Start extracting promoter regions!')
    # output_file = "GRCh38.promoter.bed"
    output_file = sub_area_file
    with open(fasta_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.replace('"', '\t')
            parts = line.split("\t")
            if len(parts) > 7 and parts[2] == "gene":
                if parts[6] == "+":
                    start = int(parts[3]) - usize
                    end = int(parts[3]) + dsize
                elif parts[6] == "-":
                    start = int(parts[4]) - dsize
                    end = int(parts[4]) + usize
                if start < 0:
                    start = 0
                if(parts[0].startswith('chr')):
                    outfile.write(f"{parts[0]}\t{start}\t{end}\t{parts[13]}\t{parts[9]}\t{parts[6]}\n")
                else:
                    outfile.write(f"chr{parts[0]}\t{start}\t{end}\t{parts[13]}\t{parts[9]}\t{parts[6]}\n")
    # print('Extracting the promoter regions have been completed!')
