#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime:2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1 
# Update 2023-12-14: we have preprocessed TE-derived-peak with string operations to facilitate annotatePeaks.pl of homer

def fix_winbed_content(file_name, query_prefix):
    file = open(file_name, 'r')
    result_fixed = open(f"{query_prefix}_homer.bed", 'w')
    lines = file.readlines()
    for line in lines:
        splited = line.split('\t')
        origin_split = splited[3].split(":")
        newline = f"{origin_split[0]}\t{origin_split[1]}\t{origin_split[2]}\t{splited[-1].strip()}&{splited[-4].strip()}&{splited[-3].strip()}&{splited[-2].strip()}\t.\n"
        result_fixed.write(newline)
    result_fixed.flush()
    result_fixed.close()
    file.close()
    return f"{query_prefix}_homer.bed"

# if __name__ == '__main__':
#     fix_content('result2.xlsx_winbed')
