#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime: 2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1

import argparse
import re
import sys
from Bio import SeqIO

# 1. Extract the gaps from genome sequence

def process_origin_file(fasta_file, out_n_prefix):
    print('Prepare to generate gap file')
    # Open FASTA, search for masked regions, print in GFF3 format
    n_bed = open(f'{out_n_prefix}_N.bed', 'w')
   
    with open(fasta_file) as handle:
        i = 0
        for record in SeqIO.parse(handle, "fasta"):
            for match in re.finditer('[Nn]+', str(record.seq)):
                if record.id.startswith('chr'):
                    n_bed.write(f"{record.id}\t{match.start() + 1}\t{match.end()}\n")
                else:
                    n_bed.write(f"chr{record.id}\t{match.start() + 1}\t{match.end()}\n")
    print('The gap file has been generated from teena!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--db_file',  required=True)
    parser.add_argument('-fa', '--fasta', required=True)
    args = parser.parse_args()
    db_prefix = args.db_file[0: args.db_file.find('.')]
    process_origin_file(args.fasta, db_prefix)
