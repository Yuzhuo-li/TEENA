#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime: 2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-15: We have fixed the bug that -log10(P) would not show up in barplot when p value is 0 after correction of benjamini_hochberg method
# Update 2024-02-08: Solved the problem of divisor being zero in Fisher's exact test
# Update 2024-06-06: Fixed bugs to make results more accurate

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from openpyxl import Workbook
import time
import re
from multiprocessing import Pool, cpu_count

def read_chrom_content1(file_name: str):
    chr_size_dict = {}
    genome_size = 0
    with open(file_name) as chr_file:
        for line in chr_file:
            splits = re.split(r'\s+', line.strip())
            if 'chr' in splits[0]:
                chr_size_dict[splits[0]] = int(splits[1])
            else:
                chr_size_dict['chr' + splits[0]] = int(splits[1])
            genome_size += int(splits[1])
    return chr_size_dict, genome_size

def read_chrom_content(file_name: str):
    chr_size_dict = {}
    genome_size = 0
    with open(file_name) as chr_file:
        lines = chr_file.readlines()
        for line in lines:
            splits = line.split('\t')
            if 'chr' in splits[0]:
                chr_size_dict[splits[0]] = int(splits[1])
            else:
                chr_size_dict['chr' + splits[0]] = int(splits[1])
            genome_size += int(splits[1])
    return chr_size_dict, genome_size

def filter_file(infile, chr_dict):
    filtered_lines = []
    with open(infile) as in_file:
        for line in in_file:
            chrom = re.split(r'\s+', line.strip())[0]
            if chrom in chr_dict:
                filtered_lines.append(line)
    return filtered_lines

def get_midpoint(filtered_lines):
    midpoint_lines = []
    for line in filtered_lines:
        splits = re.split(r'\s+', line.strip())
        if len(splits) < 3:
            continue
        try:
            mid = int((int(splits[1]) + int(splits[2])) / 2)
            splits[1] = str(mid)
            splits[2] = str(mid + 1)
            midpoint_lines.append('\t'.join(splits) + '\n')
        except ValueError:
            continue
    return midpoint_lines

def calculate_fisher(n11, n12, n21, n22):
    table = [[n11, n12], [n21, n22]]
    _, p_value_two_tail = fisher_exact(table, alternative='two-sided')
    _, p_value_one_tail1 = fisher_exact(table, alternative='less')
    _, p_value_one_tail2 = fisher_exact(table, alternative='greater')
    return p_value_one_tail1, p_value_one_tail2, p_value_two_tail

def process_family(args):
    # family, lines, pks_intervals, genome_size = args
    dna, chr, family, lines, pks_intervals, genome_size, midpoint, chr_size_dict = args
    query_counts = len(lines)
    db_counts = len(pks_intervals)

    # Calculate the union lengths
    db_union = sum(int(interval[2]) - int(interval[1]) for interval in pks_intervals)
    query_union = sum(int(re.split(r'\s+', line.strip())[2]) - int(re.split(r'\s+', line.strip())[1]) for line in lines)

    # Convert intervals to a list of tuples for easier processing
    query_intervals = [(int(re.split(r'\s+', line.strip())[1]), int(re.split(r'\s+', line.strip())[2])) for line in lines]
    pks_intervals = [(int(interval[1]), int(interval[2])) for interval in pks_intervals]

    # Sort both lists of intervals
    query_intervals.sort()
    pks_intervals.sort()

    # Use two pointers to find overlapping intervals
    # i, j = 0, 0
    # while i < len(query_intervals) and j < len(pks_intervals):
    #     query_start, query_end = query_intervals[i]
    #     pks_start, pks_end = pks_intervals[j]

    #     if query_end < pks_start:
    #         i += 1
    #     elif pks_end < query_start:
    #         j += 1
    #     else:
    #         # Overlapping intervals found
    #         # Ensure that overlaps are counted accurately
    #         overlap_start = max(query_start, pks_start)
    #         overlap_end = min(query_end, pks_end)
    #         if overlap_start < overlap_end:  # strict less than to ensure proper counting
    #             overlap_counts += 1
    #         if query_end <= pks_end:
    #             i += 1
    #         if pks_end <= query_end:
    #             j += 1
    overlap_counts = read_dna_count(dna, chr, chr_size_dict, midpoint)
    n11 = overlap_counts
    n12 = query_counts - overlap_counts
    n21 = db_counts - overlap_counts
    n22_full_bases = genome_size - query_union - db_union + overlap_counts
    d_mean = 1 + db_union / db_counts
    q_mean = 1 + query_union / query_counts
    b_mean = (q_mean + d_mean)
    n22 = max(0, int(n22_full_bases / b_mean) - n12 - n21 - n11)

    p_value_one_tail1, p_value_one_tail2, p_value_two_tail = calculate_fisher(n11, n12, n21, n22)
    foldenrich = (n11 / (n11 + n12)) / (n21 / (n21 + n22))
    # print('fisher item finished')
    return family, {
        "p_both": p_value_two_tail,
        "FoldEnrich": foldenrich,
        "a1_b1": n11,
        "a1_b0": n21,
        "a0_b1": n12,
        "a0_b0": n22,
        "p_left": p_value_one_tail1,
        "p_right": p_value_one_tail2,
        "Overlap_fraction": n11 / (n11 + n12)
    }

def read_dna_count(dna: dict, chr: dict, chr_size_dict: dict, midpoint=True):
    n11 = 0
    for key in dna.keys():
        n11_item = 0
        if key in chr.keys() and key in chr_size_dict.keys():
            intervals1 = []
            intervals2 = []
            dna_chr_list = dna[key]
            chr_list = chr[key]
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

            intersect_count = count_intersections(intervals1, intervals2, midpoint)

            n11_item = intersect_count
        n11 += n11_item
    return n11

def count_intersections(intervals1, intervals2, midpoint=True):
    intervals1.sort(key=lambda x: x[0])
    intervals2.sort(key=lambda x: x[0])

    i, j = 0, 0
    intersections = 0

    while i < len(intervals1) and j < len(intervals2):
        start1, end1 = intervals1[i]
        start2, end2 = intervals2[j]

        if end1 < start2:
            i += 1
        elif end2 < start1:
            j += 1
        else:
            if midpoint:
                mid = (start2 + end2) / 2
                if mid >= start1 and mid <= end1:
                    intersections += 1
            else:
                intersections += 1
            if end1 < end2:
                i += 1
            else:
                j += 1

    return intersections

def run_fisher(f_dna, f_chr, chr_size_dict, genome_size, rep_lines, pks_lines, midpoint=False, verbose=True):
    filtered_rep_lines = filter_file(rep_lines, chr_size_dict)
    filtered_pks_lines = filter_file(pks_lines, chr_size_dict)
    # print("file filter finished........")
    if midpoint:
        filtered_pks_lines = get_midpoint(filtered_pks_lines)
        # print("mid point  finished........")

    reps = {}
    for line in filtered_rep_lines:
        family = re.split(r'\s+', line.strip())[3]
        if family in reps:
            reps[family].append(line)
        else:
            reps[family] = [line]

    pks_intervals = [re.split(r'\s+', line.strip()) for line in filtered_pks_lines]

    args_list = []
    for family, lines in reps.items():
        if len(lines) >= 10:
            # args_list.append((family, lines, pks_intervals, genome_size))
            args_list.append((f_dna[family], f_chr, family, lines, pks_intervals, genome_size, midpoint, chr_size_dict))

        else:
            print('less 10')
    # args_list = [(family, lines, pks_intervals, genome_size) for family, lines in reps.items() if len(lines) >= 10]
    
    results = {}
    for args in args_list:
        key, value = process_family(args)
        results[key] = value
    return results

def fisher_valid(chrom_size_file, query_file, db_file, out_put, midpoint=False):
    _, genome_size = read_chrom_content1(chrom_size_file)
    chr_size_dict, _ = read_chrom_content(chrom_size_file)
    # chr_size_dict, genome_size = read_chrom_content(chrom_size_file)
    with open(query_file) as file_1, open(db_file) as file_2:
        f1_lines = file_1.readlines()
        f2_lines = file_2.readlines()

    query_file_length = 0

    f_dna = {}
    for f2_line in f2_lines:
        if f2_line.strip() == '':
            continue
        arr = re.split(r'\s+', f2_line.strip())
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
        query_file_length += 1
        arr = re.split(r'\s+', f1_line.strip())
        if arr[0] in f_chr.keys():
            f_chr[arr[0]].append(arr)
        else:
            f_chr[arr[0]] = [arr]

    for f_key in f_chr.keys():
        f_chr[f_key].sort(key=lambda x: int(x[1]))

    fisher_results = run_fisher(f_dna, f_chr,chr_size_dict, genome_size, db_file, query_file, midpoint)

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
        family = f_key
        if family in fisher_results:
            result = fisher_results[family]
            name = f_key
            n11 = result['a1_b1']
            n12 = result['a1_b0']
            n21 = result['a0_b1']
            n22 = result['a0_b0']
            p_value_one_tail1 = result['p_left']
            p_value_one_tail2 = result['p_right']
            p_value_two_tail = result['p_both']
            foldenrich = result['FoldEnrich']
            keyword = name.split(':')[1].strip()
            te_class = name.split('/')[0].strip()
            if te_class == name.strip():
                te_class = name.split(':')[0].strip()
                te_family = ''
            else:
                te_family = name.split(':')[0].split('/')[1].strip()
            index += 1
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
        workbook.save(out_put)
        return out_put
    else:
        workbook.save(f'{out_put}.xlsx')
        return f'{out_put}.xlsx'

def benjamini_hochberg(p_values, threshold=1e-250):
    m = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = np.sort(p_values)
    adjusted_p_values = m * sorted_p_values / (np.arange(m) + 1)
    adjusted_p_values = np.minimum.accumulate(adjusted_p_values[::-1])[::-1]
    adjusted_p_values[adjusted_p_values < threshold] = threshold
    return adjusted_p_values[np.argsort(sorted_indices)]

def update_adj(middle_file, out_put, query_file):
    print(f"Reading file {middle_file}")
    d = pd.read_excel(middle_file)
    print("Columns in file:", d.columns.tolist())
    d['p_adj'] = benjamini_hochberg(d['p_value'].values, threshold=1e-250)
    new_df = pd.DataFrame()
    new_df['TE_name'] = d['TE_name']
    new_df['TE_family'] = d['TE_family']
    new_df['TE_class'] = d['TE_class']
    new_df['overlap_number'] = d['a1_b1']
    # new_df['a1_b0'] = d['a1_b0']
    # new_df['a0_b1'] = d['a0_b1']
    # new_df['a0_b0'] = d['a0_b0']
    new_df['p_value'] = d['p_value']
    new_df['p_adj'] = d['p_adj']
    new_df['fold_enrich'] = d['fold_enrich']
    new_df['Dfam_links'] = d['Dfam_link']
    new_df['Repbase_links'] = d['Repbase']
    new_df_sorted = new_df.sort_values(by='p_value')
    d.to_excel(middle_file, index=False)
    if out_put.endswith('xlsx'):
        # print(f' Output results {out_put}')
        new_df_sorted.to_excel(out_put, index=False)
    else:
        # print(f' Output results {out_put}_result.xlsx')
        new_df_sorted.to_excel(f'{out_put}_TEENA_result.xlsx', index=False)

    df_sign_sorted = new_df.sort_values(by='p_adj')
    if df_sign_sorted['p_adj'].min() < 0.05:
        df_sign_sorted = df_sign_sorted[df_sign_sorted['p_adj'] < 0.05]
        query_prefix = query_file.split('_')[0]
        # print(f' Output results {out_put}_significant.xlsx')
        df_sign_sorted.to_excel(f'{out_put}_TEENA_significant.xlsx', index=False)
    else:
        print('No significant TE identified')

    if out_put.endswith('xlsx'):
        return out_put
    else:
        return f'{out_put}.xlsx'

def run(chrom_size_file, query_file, db_file, out_put, midpoint=False):
    # print("11111111111")
    start_time = time.time()
    middle_file = fisher_valid(chrom_size_file, query_file, db_file, out_put, midpoint)
    # print(f"Fisher test output saved to {middle_file}")

    out_put_file_name = update_adj(middle_file, out_put, query_file)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"The Fisher's exact test has been completed! Code runtime: {elapsed_time:.6f} seconds")
    return middle_file, out_put_file_name

# if __name__ == '__main__':
    # run('mm10.chrom.sizes', 'K27ac_24h.gain.bed','mm10.repbase.bed',  'CTCF.FET.xlsx', midpoint=True)
