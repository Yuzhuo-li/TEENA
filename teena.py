#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime:2023-11-23
# Software: TEENA:determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-11: We have added  homer annotation function to TE-derived-peak, and integrated the neighboring genes information of overlaps from query file and TE in Excel sheet.
# Update 2024-1-17: Retrieved the last string from the query file for renaming
# Update 2024-01-20: Read species files which named species.txt
# Update 2024-02-08: Solved the problem of divisor being zero in Fisher's exact test
# Update 2024-02-08: Fixed pattern issues in string columns
# Update 2024-02-09: If the number of species annotation is zero, it means null, and directly returns without drawing a graph
# Update 2024-02-18: Fixed the bug of extracting database file name during TE annotation

import argparse
# from metrics import calculate_metrics
import os
import platform
import pandas as pd
from fix_winbed_result import fix_winbed_content
from homer_to_excel import convert_homer_result_to_excel
from draw_piechart import draw_pie
import re
from find_n import process_origin_file
from sub_start_area import get_sub_result
from fisher import run
import subprocess
from plot_img import draw_plot

def fix_file_content(file):
    file_open = open(file, 'r')
    file_fixed = open(f"{file}.fixed", 'w')
    lines = file_open.readlines()
    col = 0
    for line in lines:
        if line.startswith('#'):
            file_fixed.write(line)
        if not line.startswith('chr'):
            file_fixed.write('chr' + line)
        else:
            file_fixed.write(line)
        col = len(line.split('\t'))
    file_fixed.flush()
    file_fixed.close()
    return col, f"{file}.fixed"



def get_top_col(col, file, query_prefix):
    file_open = open(file, 'r')
    top_file_fixed = open(f"{query_prefix}_top_{col}.fixed", 'w')
    lines = file_open.readlines()
    for line in lines:
        splits = re.split(r'\s+', line)
        middle = round((int(splits[1]) + int(splits[2])) / 2)
        m_left = middle - 1
        m_right = middle + 1
        newline = f"{splits[0]}\t{m_left}\t{m_right}\t{':'.join(splits[0:3])}\n"
        top_file_fixed.write(newline)
    top_file_fixed.flush()
    top_file_fixed.close()
    return f"{query_prefix}_top_{col}.fixed"


def fix_repbase_content(file):
    file_open = open(file, 'r')
    file_fixed = open(f"{file}.fixed", 'w')
    lines = file_open.readlines()
    for line in lines:
        line_split = line.split('\t')
        if line_split[3].startswith('DNA') or line_split[3].startswith('LINE') or line_split[3].startswith('SINE') or line_split[3].startswith('LTR'):
            file_fixed.write('\t'.join(line_split[0: 4]) + '\n')
    file_fixed.flush()
    file_fixed.close()
    return f"{file}.fixed"
            

def execute_homer(query_file,  homer_gtf_file, out_put, homer_param='hg38'):
    # run_bedtools_subtract(["annotatePeaks.pl", f"{query_file}", f"{homer_param}", ">", f"{out_put}_homer.txt"])
    run_bedtools_subtract(["annotatePeaks.pl", f"{query_file}", f"{homer_param.split('/')[-1]}", ">", f"{out_put}_homer.txt"])
    return f"{out_put}_homer.txt"

def run_bedtools_subtract(cmd):
        print(" ".join(cmd))
        subprocess.run(" ".join(cmd), shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Name\n

    TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test\n

Description\n

    This script performs fisher's exact test on the number of overlaps/unique intervals between 2 files. By removing the gap and promoter files, pipelines or analyses can be run on the disturbed and original files, and then compared.This script provides function for data filtering, TEENA analysis, scatter plot, bar plot, pie plot, TE-derived-peak annotation and result summarizing.\n

Version\n
    Author: Yuzhuo Li (DX120230210@stu.yzu.edu.cn) \n
    Version: 1.1    Date: 2023-11-23 \n

Usage\n
    python teena.py [parameters] -q query_bed_file -d te_bed_file -a gtf_file -n gap_file -ch chrom_size_file -o result_file \n

Exmple\n
    python teena.py -q GATA3_hg38.bed -d hg38.repbase.bed -a Homo_sapiens.GRCh38.110.gtf -n hg38_N.bed -ch hg38.chrom.sizes -o GATA3 """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-q', '--query_file', required=True, help="The query bed file you want to analyze. Your query bed file only needs to include the chromosome name and its starting and ending positions in three columns(required=True).")
    parser.add_argument('-d', '--db_file',  required=True, help="The repbase annotation bed file we have already downloaded from UCSC, the first three columns are the same as the query file, but the 4th column contains the TE family names(required=True).")
    parser.add_argument('-uk', '--usize', type=int, required=False, default=500, help="Upstream from transcription initiation site(default=500, required=False).")
    parser.add_argument('-dk', '--dsize', type=int, required=False, default=500, help="Downstream from transcription initiation site(default=500, required=False).")
    parser.add_argument('-o', '--output_file', required=True, help="Result file after running the process, and just enter the prefix name which you want, and the suffix xlsx will be automatically added(required=True).")
    parser.add_argument('-a', '--gtf', required=True, help="The GTF annotation file we have already downloaded from Ensembl(required=True).")
    parser.add_argument('-fa', '--fasta', required=False, help="The genome sequence file we have already downloaded from UCSC(required=False).")
    parser.add_argument('-ch', '--chrom_file', required=True, help="The chromosome length file we have already downloaded from UCSC(required=True).")
    parser.add_argument('-m', '--midpoint', type=bool, required=False, default=True, help="We believe that the midpoint of the query file is on the interval of the corresponding annotation file, that is, two intervals overlap; If your query file is certain broad-peak file(histone modified like H3K9me3, H3K27me3 etc.), then we consider the intersection of the comment file interval and the corresponding interval of the query file as overlap(required=False, default='True').")
    parser.add_argument('-n', '--nbed', required=False, default=None, help="The bed file of gap from genome sequence which we have preprocessed(required=False, default=None).")
    parser.add_argument('-rn', '--remove_n', type=str, required=False, default='True', help="Choose whether to remove the gaps from genome sequence(required=False, default='True').")
    parser.add_argument('-rp', '--remove_p', type=str, required=False, default='True', help="Choose whether to remove the promoter regions(required=False, default='True').")
    parser.add_argument('-V', '--version', action='version', help="Author: Yuzhuo Li (DX120230210@stu.yzu.edu.cn)   Version: 1.1   Date: 2023-11-23")
    args = parser.parse_args()


    db_prefix = args.db_file[0: args.db_file.find('.')]
    remove_n = args.remove_n
    remove_p = args.remove_p
    ## Generate gap file

    bed_n_file = args.nbed
    if not bed_n_file:
        if args.fasta:
            print('The gap file is generated through the teena.py script')
            process_origin_file(args.fasta, db_prefix)
            bed_n_file = f'{db_prefix}_N.bed'
        else:
            print('Please pass in gap file or fasta file to generate gap file!')
            exit(0)
    
    ## Extract the promoter regions (The promoter region refers to the upstream and downstream 500bp of the transcription start site, which is removed by default. This parameter can choose to remove more or fewer promoter regions)
    _, gtf_fixed = fix_file_content(args.gtf)

    # print('Extract the promoter regions!')

    get_sub_result(gtf_fixed, f'{db_prefix}.promoter.bed', args.usize, args.dsize)
   
    ## Use bedtools subtract for subtraction operation (subtraction operation, subtract the area in B from the area in A, and eliminate all overlapping areas through parameter - A)

    # Perform additional preprocessing on bed files,to prevent extra tab at the end of your line.
    try:
        if platform.system() == 'Darwin':
            run_bedtools_subtract(["sed", "-i", "''", "'s/\t*$//'", f"{args.db_file}"])
            run_bedtools_subtract(["sed", "-i", "''", "'s/\t*$//'", f"{args.query_file}"])
        else:
            run_bedtools_subtract(["sed", "-i", "'s/\t*$//'", f"{args.db_file}"])
            run_bedtools_subtract(["sed", "-i", "'s/\t*$//'", f"{args.query_file}"])
    except Exception as ex:
        print(f"An error occurred while executing the sed command: {ex}")   

    fixed_db_file = fix_repbase_content(args.db_file)
    # Update 2024-01-17: Retrieve the last string from the query file for renaming
    query_prefix = args.query_file.split(".")[-2].split('_') [0]
    # query_prefix = args.query_file.split('_')[0]
    col, query_file_fixed = fix_file_content(args.query_file)
    middle_file = ''
    out_put_file_name = ''
    # Remove gaps and promoter regions from genome sequence
    if remove_n == 'True' and remove_p == 'True':
        run_bedtools_subtract(["bedtools", "subtract", "-a", fixed_db_file, "-b", bed_n_file, ">", f"{db_prefix}_repbese_no_N.bed"])
        run_bedtools_subtract(["bedtools", "subtract",  "-a", f"{db_prefix}_repbese_no_N.bed", "-b", f'{db_prefix}.promoter.bed', ">", f"{db_prefix}_repbese_no_N_no_promoter.bed"])
        # calculate db 
        # print('calculate db ')
        # calculate_metrics(f"{db_prefix}_metrics.txt", f"{db_prefix}_repbese_no_N_no_promoter.bed", args.db_file);

       
        run_bedtools_subtract(["bedtools", "subtract", "-a", query_file_fixed, "-b", bed_n_file, ">", f"{query_prefix}_no_N.bed"])
        run_bedtools_subtract(["bedtools", "subtract",  "-a", f"{query_prefix}_no_N.bed", "-b", f'{db_prefix}.promoter.bed', ">", f"{query_prefix}_no_N_no_promoter.bed"])
        # print('calculate query ')
        # calculate_metrics(f"{query_prefix}_metrics.txt", f"{query_prefix}_no_N_no_promoter.bed", args.query_file);
    
        ## Remove gaps from genome sequence
        middle_file, out_put_file_name = run(args.chrom_file, f"{query_prefix}_no_N_no_promoter.bed", f"{db_prefix}_repbese_no_N_no_promoter.bed", args.output_file, args.midpoint)
    elif remove_n == 'True' and remove_p == 'False':
        run_bedtools_subtract(["bedtools", "subtract", "-a", fixed_db_file, "-b", bed_n_file, ">", f"{db_prefix}_repbese_no_N.bed"])
        # run_bedtools_subtract(["bedtools", "subtract",  "-a", f"{db_prefix}_repbese_no_N.bed", "-b", f'{db_prefix}.promoter.bed', ">", f"{db_prefix}_repbese_no_N_no_promoter.bed"])
        # calculate db 
        # print('calculate db ')
        # calculate_metrics(f"{db_prefix}_metrics.txt", f"{db_prefix}_repbese_no_N.bed", args.db_file);

        run_bedtools_subtract(["bedtools", "subtract", "-a", query_file_fixed, "-b", bed_n_file, ">", f"{query_prefix}_no_N.bed"])
        # run_bedtools_subtract(["bedtools", "subtract",  "-a", f"{query_prefix}_no_N.bed", "-b", f'{db_prefix}.promoter.bed', ">", f"{query_prefix}_no_N_no_promoter.bed"])
        # print('calculate query ')
        # calculate_metrics(f"{query_prefix}_metrics.txt", f"{query_prefix}_no_N.bed", args.query_file);
    

        ## Remove promoter regions from genome sequence
        middle_file, out_put_file_name = run(args.chrom_file, f"{query_prefix}_no_N.bed", f"{db_prefix}_repbese_no_N.bed", args.output_file, args.midpoint)
    elif remove_n == 'False' and remove_p == 'True':
        # run_bedtools_subtract(["bedtools", "subtract", "-a", fixed_db_file, "-b", bed_n_file, ">", f"{db_prefix}_repbese_no_N.bed"])
        run_bedtools_subtract(["bedtools", "subtract",  "-a", fixed_db_file, "-b", f'{db_prefix}.promoter.bed', ">", f"{db_prefix}_repbese_no_promoter.bed"])
        # calculate db 
        # print('calculate db ')
        # calculate_metrics(f"{db_prefix}_metrics.txt", f"{db_prefix}_repbese_no_promoter.bed", args.db_file);

        # run_bedtools_subtract(["bedtools", "subtract", "-a", query_file_fixed, "-b", bed_n_file, ">", f"{query_prefix}_no_N.bed"])
        run_bedtools_subtract(["bedtools", "subtract",  "-a", query_file_fixed, "-b", f'{db_prefix}.promoter.bed', ">", f"{query_prefix}_no_promoter.bed"])
        # print('calculate query ')
        # calculate_metrics(f"{query_prefix}_metrics.txt", f"{query_prefix}_no_promoter.bed", args.query_file);
    
        middle_file, out_put_file_name = run(args.chrom_file, f"{query_prefix}_no_promoter.bed", f"{db_prefix}_repbese_no_promoter.bed", args.output_file, args.midpoint)
    else:
        # run_bedtools_subtract(["bedtools", "subtract", "-a", fixed_db_file, "-b", bed_n_file, ">", f"{db_prefix}_repbese_no_N.bed"])
        # run_bedtools_subtract(["bedtools", "subtract",  "-a", f"{db_prefix}_repbese_no_N.bed", "-b", f'{db_prefix}.promoter.bed', ">", f"{db_prefix}_repbese_no_N_no_promoter.bed"])
        # # calculate db 
        # print('calculate db ')
        # calculate_metrics(f"{db_prefix}_metrics.txt", f"{db_prefix}_repbese_no_N_no_promoter.bed", args.db_file);

        # run_bedtools_subtract(["bedtools", "subtract", "-a", query_file_fixed, "-b", bed_n_file, ">", f"{query_prefix}_no_N.bed"])
        # run_bedtools_subtract(["bedtools", "subtract",  "-a", f"{query_prefix}_no_N.bed", "-b", f'{db_prefix}.promoter.bed', ">", f"{query_prefix}_no_N_no_promoter.bed"])
        # print('calculate query ')
        # calculate_metrics(f"{query_prefix}_metrics.txt", f"{query_prefix}_no_N_no_promoter.bed", args.query_file);
    
        middle_file, out_put_file_name = run(args.chrom_file, query_file_fixed, fixed_db_file, args.output_file, args.midpoint)
    
    # Plot image
    print('start plot image!')
    if os.path.exists(middle_file):
        draw_plot(middle_file, out_put_file_name.split('.')[0])
        os.remove(middle_file)
    
    if col >= 3:
        print(f'Intercept the first 3 columns of the {query_file_fixed} file')
        top_col_file_name = get_top_col(3, query_file_fixed, query_prefix)
        # Executing windowBed
        
        run_bedtools_subtract(["windowBed", "-a", f"{top_col_file_name}", "-b", f"{fixed_db_file}", "-w", "0", ">", f"{query_prefix}_mid.bed"])

        winbed_fixed_file_name = fix_winbed_content(f"{query_prefix}_mid.bed", query_prefix)
        print(f'The file correction before homer annotation has been completed! Your file name is {winbed_fixed_file_name}')

        print("Annotate TE-derived-peak with homer annotatePeaks.pl")

        
        homer_result_txt_file = execute_homer(query_file=winbed_fixed_file_name, homer_gtf_file=gtf_fixed, out_put=out_put_file_name.split('.')[0], homer_param=db_prefix)
        print("")
        # annotation_list = convert_homer_result_to_excel('demo_output_homer.txt', 'demo_output')
        # annotation_pds = convert_homer_result_to_excel(homer_result_txt_file, out_put_file_name.split('.')[0])
        # Update 2024-01-17: Read species files which named species.txt
        species_reader = open('species.txt', 'r')
        line_striped = []
        lines = species_reader.readlines()
        for line in lines:
            line_striped.append(line.strip())
        spec_exist = args.db_file.split('/')[-1][0: args.db_file.split('/')[-1].find('.')] in line_striped
        annotation_pds = convert_homer_result_to_excel(homer_result_txt_file, out_put_file_name.split('.')[0], spec_exist)
        print("Extract and integrate the annotated txt file from homer and generate a new Excel sheet")
        # Pie chart generation
        # draw_pie(annotation_list, 'demo_output')
        print("address lllllogs ....................................................................................")
        if spec_exist:
            draw_pie(annotation_pds, out_put_file_name.split('.')[0])
