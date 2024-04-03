#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime: 2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-15: We have integrated the results obtained from the annotation of homer with TE-derived-peak
# Update 2023-12-19: the header of the generated xlsx file is changed
# Update 2024-01-17: Added a species files which named species.txt
# Update 2024-02-08: Fixed pattern issues in string columns


import pandas as pd

def convert_homer_result_to_excel(homer_txt_file, res_file, spec_exist):
    df_with_header = pd.read_csv(homer_txt_file, sep='\t')
    split_columns = df_with_header[df_with_header.columns[0]].str.split('&', expand=True)
    split_columns.columns = ["TE", "TE_chr", "TE_start", "TE_end"]
    df_with_header = pd.concat([df_with_header, split_columns], axis=1)
    new_df = pd.DataFrame()
    new_df['chr'] = df_with_header['Chr']
    new_df['start'] = df_with_header['Start']
    new_df['end'] = df_with_header['End']
    new_df['TE'] = df_with_header['TE']
    new_df['TE_chr'] = df_with_header['TE_chr']
    new_df['TE_start'] = df_with_header['TE_start']
    new_df['TE_end'] = df_with_header['TE_end']
    print(df_with_header['Annotation'])
    if spec_exist:
        df_with_header['Annotation'] = df_with_header['Annotation'].apply(lambda x: x.split('(')[0] if pd.notnull(x) else x)
        new_df['Gene_name'] = df_with_header['Gene Name']
        new_df['Gene_ID'] = df_with_header['Nearest Ensembl']
        new_df['Distance_to_TSS'] = df_with_header['Distance to TSS']
        new_df['Distribution'] = df_with_header['Annotation']
    new_df = new_df.sort_values(by= ['chr', 'start'], ascending=[True, True])
    new_df.to_excel(f'{res_file}_TEENA_TE-derived-peak_annotation.xlsx', index=False)
    return new_df
    
#if __name__ == '__main__':
#    species_reader = open('species.txt', 'r')
#    lines = species_reader.readlines()
#    line_striped = []
#    for linr in lines:
#        line_striped.append(linr.strip())
#    spec_exist =  'hg38' in line_striped
#    convert_homer_result_to_excel('2_homer.txt', "ttt2", spec_exist)
