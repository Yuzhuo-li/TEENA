#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime: 2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-15: We have fixed the bug that -log10(P) would not show up in barplot when p value is 0 after correction of benjamini_hochberg method
# Update 2024-02-08: Solved the problem of divisor being zero in Fisher's exact test

from scipy.stats import fisher_exact
import sys
import time
import openpyxl
from openpyxl import Workbook
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.hyperlink import Hyperlink
import pandas as pd
import numpy as np
import re

def read_chrom_content(file_name: str):
    chr_size_dict = {}
    chr_file = open(file_name)
    lines = chr_file.readlines()
    for line in lines:
        splits = line.split('\t')
        if 'chr' in splits[0]:
            chr_size_dict[splits[0]] = int(splits[1])
        else:
            chr_size_dict['chr' + splits[0]] = int(splits[1])
    return chr_size_dict



def fisher_valid(chrom_size_file, query_file, db_file, out_put:str, midpoint = True):
    # chr_size_dict = read_chrom_content('GRCh38.chrom.sizes')
    # file_1 = open('human_s_no_N_no_promoter.bed')
    # file_2 = open('hg38_repbese_no_N_no_promoter.bed')

    chr_size_dict = read_chrom_content(chrom_size_file)
    file_1 = open(query_file)
    file_2 = open(db_file)

    f1_lines = file_1.readlines()
    f2_lines = file_2.readlines()

    f_dna = {}
    for f2_line in f2_lines:
        if f2_line.strip() == '':
            continue
        arr = f2_line.split('\t')
        #print(arr)
        if arr[3] in f_dna.keys():
            if arr[0] in f_dna[arr[3]].keys():
                f_dna[arr[3]][arr[0]].append(arr)
            else:
                f_dna[arr[3]][arr[0]] = [arr]
        else:
            f_dna[arr[3]] = {arr[0]: [arr]}

    for f_key in f_dna.keys():
        for chr_key in f_dna[f_key]:
            f_dna[f_key][chr_key].sort(key=lambda x: int(x[1]))

    f_chr = {}
    for f1_line in f1_lines:
        arr = f1_line.split('\t')
        if arr[0] in f_chr.keys():
            f_chr[arr[0]].append(arr)
        else:
            f_chr[arr[0]] = [arr]

    for f_key in f_chr.keys():
        f_chr[f_key].sort(key=lambda x: int(x[1]))

    ## Output to file
    workbook = Workbook()
    sheet = workbook.active
    sheet['A1'] = 'TE_name'
    sheet['B1'] = 'TE_family'
    sheet['C1'] = 'TE_class'
    sheet['D1'] = 'a1_b1'
    sheet['E1'] = 'a1_b0'
    sheet['F1'] = 'a0_b1'
    sheet['G1'] = 'a0_b0'
    sheet['H1'] = 'p_left'
    sheet['I1'] = 'p_value'
    sheet['J1'] = 'p_both'
    sheet['K1'] = 'fold_enrich'
    sheet['L1'] = 'Dfam_link'
    sheet['M1'] = 'Repbase'
    index = 1
    for f_key in f_dna.keys():
        name, n11, n12, n21, n22, p_value_one_tail1, p_value_one_tail2, p_value_two_tail,foldenrich = read_dna_count(f_key, f_dna[f_key], f_chr, chr_size_dict, midpoint)
        keyword = name.split(':')[1].strip()
        te_class = name.split('/')[0].strip()
        if  te_class == name.strip():
            te_class = name.split(':')[0].strip()
            te_family = ''
        else:
            te_family = name.split(':')[0].split('/')[1].strip()
        index = index + 1
        sheet[f'A{index}'] = keyword
        sheet[f'B{index}'] = te_family
        sheet[f'C{index}'] = te_class
        sheet[f'D{index}'] = n11
        sheet[f'E{index}'] = n12
        sheet[f'F{index}'] = n21
        sheet[f'G{index}'] = n22
        sheet[f'H{index}'] = p_value_one_tail1
        sheet[f'I{index}'] = p_value_one_tail2
        sheet[f'J{index}'] = p_value_two_tail
        sheet[f'K{index}'] = foldenrich
        
        link_url = f'https://dfam.org/browse?keywords={keyword}'
        sheet[f'L{index}'] = 'Dfam_link'
        cell = sheet[f'L{index}'] 
        cell.value = link_url
        cell.hyperlink = link_url  

        rep_link_url = f'https://www.girinst.org/repbase/update/search.php?query={keyword}'
        sheet[f'M{index}'] = 'Repbase_link'
        cell1 = sheet[f'M{index}'] 
        cell1.value = rep_link_url
        cell1.hyperlink = rep_link_url  
        
    if out_put.endswith('xlsx'):
        workbook.save(f'{out_put}')
        return f'{out_put}'
    else:
        workbook.save(f'{out_put}.xlsx')
        return f'{out_put}.xlsx'
    


# intervals2 query
# intervals1 db
def count_intersections(intervals1, intervals2,  midpoint = True):
    intervals1.sort(key=lambda x: x[0])  # Sort the intervals of file 1 by starting coordinates
    intervals2.sort(key=lambda x: x[0])  # Sort the intervals of file 2 by starting coordinates

    i, j = 0, 0  # i points to the interval of file 1, j points to the interval of file 2
    intersections = 0

    while i < len(intervals1) and j < len(intervals2):
        start1, end1 = intervals1[i]
        start2, end2 = intervals2[j]

        # If the end coordinate of file 1 is smaller than the start coordinate of file 2, move the pointer of file 1
        if end1 < start2:
            i += 1
        # If the end coordinate of file 2 is smaller than the start coordinate of file 1, move the pointer of file 2
        elif end2 < start1:
            j += 1
        # Otherwise, the two intervals intersect
        else:
            # If the parameter is true, the midpoint needs to be compared
            if midpoint:
                mid = (start2 + end2) / 2
                if mid >= start1 and mid <= end1:
                    intersections += 1
            else:
                # Otherwise, as long as it intersects, it is considered usable
                intersections += 1
            # To avoid duplicate counting, only move one pointer
            if end1 < end2:
                i += 1
            else:
                j += 1

    return intersections


def read_dna_count(name: str, dna: dict, chr: dict, chr_size_dict: dict, midpoint = True):
    n11 = 0
    n12 = 0
    n21 = 0
    n22 = 0
    for key in dna.keys():
        n11_item = 0
        n12_item = 0
        n21_item = 0
        n22_item = 0
        if (key in chr.keys()) and (key in chr_size_dict.keys()):
            intervals1 = []
            intervals2 = []
            size = chr_size_dict[key]
            dna_chr_list = dna[key]
            chr_list = chr[key]
            # Scanning-line algorithm
            db_union = 0
            for dna_chr_item in dna_chr_list:
                start = int(dna_chr_item[1])
                end = int(dna_chr_item[2])
                intervals1.append([start, end])
                db_union += (end - start)

            query_union = 0
            for chr_item in chr_list:
                start = int(chr_item[1])
                end = int(chr_item[2])
                intervals2.append([start, end])
                query_union += (end - start)
            
            intersect_count = count_intersections(intervals1, intervals2,  midpoint)
            
            db_count = len(dna_chr_list)
            query_count = len(chr_list)

            n22_full_base = size
            d_mean = 1 + db_union / db_count
            q_mean = 1 + query_union / query_count
            b_mean = d_mean + q_mean
            

            n11_item = intersect_count
            n21_item = db_count - n11_item
            n12_item = query_count - n11_item
            n22_full = max(n21_item + n12_item + n11_item, n22_full_base / b_mean)
            n22_item = max(0, n22_full - n12_item - n21_item - n11_item)
        n11 += n11_item
        n21 += n21_item
        n12 += n12_item
        n22 += n22_item
    n22 = round(n22)
    

    # Create a 2x2 table
    table = [[n11, n12], [n21, n22]]

    # Two-sided test
    oddsratio, p_value_two_tail = fisher_exact(table, alternative='two-sided')
        
    # One-sided test
    oddsratio, p_value_one_tail1 = fisher_exact(table, alternative='less')

    oddsratio, p_value_one_tail2 = fisher_exact(table, alternative='greater')

    # Foldenrich
    foldenrich = calculate_foldenrich(n11, n21, n12, n22)

    # a1_b1\ta1_b0\ta0_b1\ta0_b0
    return name, n11, n12, n21, n22, p_value_one_tail1, p_value_one_tail2, p_value_two_tail,foldenrich

def calculate_foldenrich(a1_b1, a0_b1, a1_b0, a0_b0):
    obs_diff = (a1_b1 + a0_b1)
    if obs_diff == 0:
        obs_hits = 0
    else:
        obs_hits = (a1_b1) / (a1_b1 + a0_b1)
    if obs_hits == 0:
        return 0
    exp_hits = (a1_b1 + a1_b0) / (a1_b1 + a1_b0 + a0_b1 + a0_b0)
    return obs_hits / exp_hits

def benjamini_hochberg(p_values, threshold=1e-250):
    m = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = np.sort(p_values)
    adjusted_p_values = m * sorted_p_values / (np.arange(m) + 1)
    adjusted_p_values = np.minimum.accumulate(adjusted_p_values[::-1])[::-1]
    adjusted_p_values[adjusted_p_values < threshold] = threshold
    return adjusted_p_values[np.argsort(sorted_indices)]

def update_adj(middle_file, out_put, query_file):
    d = pd.read_excel(middle_file)
    d['p_adj'] = benjamini_hochberg(d['p_value'], threshold=1e-250)
    new_df = pd.DataFrame()
    new_df['TE_name'] = d['TE_name']
    new_df['TE_family'] = d['TE_family']
    new_df['TE_class'] = d['TE_class']
    new_df['overlap_number'] = d['a1_b1']
    new_df['p_value'] = d['p_value']
    new_df['p_adj'] = d['p_adj']
    new_df['fold_enrich'] = d['fold_enrich']
    new_df['Dfam_links'] = d['Dfam_link']
    new_df['Repbase_links'] = d['Repbase']
    new_df_sorted = new_df.sort_values(by='p_value')
    # Update intermediate files
    d.to_excel(middle_file, index=False)
    # Save the new DataFrame as an Excel sheet
    if out_put.endswith('xlsx'):
        new_df_sorted.to_excel(out_put, index=False)
    else:
        new_df_sorted.to_excel(f'{out_put}_TEENA_result.xlsx', index=False)
    # Save significantly enriched TEs as an Excel sheet
    df_sign_sorted = new_df.sort_values(by='p_adj')
    if df_sign_sorted['p_adj'].min() < 0.05:
        df_sign_sorted = df_sign_sorted[df_sign_sorted['p_adj'] < 0.05]
        query_prefix = query_file.split('_')[0]
        # Save
        df_sign_sorted.to_excel(f'{out_put}_TEENA_significant.xlsx', index=False)
    else:
        print('No significant TE identified')
    
    if out_put.endswith('xlsx'):
        return out_put
    else:
        return f'{out_put}.xlsx'


def run(chrom_size_file, query_file, db_file, out_put, midpoint = True):
    print("Start Fisher's exact test")
    start_time = time.time()

    # Calculate raw data of Fisher's exact test
    middle_file =  fisher_valid(chrom_size_file, query_file, db_file, out_put, midpoint)
    
    # Calculate the updated padj data
    out_put_file_name = update_adj(middle_file, out_put, query_file)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"The Fisher's exact test has been completed! Code runtime: {elapsed_time:.6f} seconds")
    return middle_file, out_put_file_name


#print('LTR:MamRep1527'.split('/')[0].strip())
