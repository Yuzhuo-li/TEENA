#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime: 2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-14: We have added a new function in teena:to calculate and count the genome distribution of neighboring genes after annotationof homer of TE-derived-peak in the query file
# Update 2024-02-09: If the number of species annotation is zero, it means null, and directly returns without drawing a graph


import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd

def draw_pie(data_pds , name):
    four_category = ['DNA', 'LINE', 'SINE', 'LTR']
    fig, axs = plt.subplots(1, 4, figsize=(12, 5))
    for cate in four_category:
        index = four_category.index(cate)
        data = data_pds[data_pds['TE'].str.startswith(cate)]['Distribution'].to_list()
        # print(data)
        counts = Counter(data)
        labels = counts.keys()
        sizes = counts.values()
        # Calculates the proportion and formats it as a string
        total = sum(sizes)
        # Create a new label containing proportion information
        dict_data = {}
        for label, size in zip(labels, sizes):
            dict_data[str(label).strip()] = size
        labels_with_pct = []
        values = []
        lbs = ["promoter-TSS", "5'UTR", "exon", "intron", "3'UTR", "TTS", "non-conding", "Intergenic" ]
        i_list =[]
        for i,lab in enumerate(lbs):
            if lab in dict_data.keys():
                val = dict_data[lab]
                values.append(dict_data[lab])
                if lab == 'promoter-TSS':
                   lab = 'promoter'
                labels_with_pct.append(f"{lab} - {round(val * 100 / total, 2)}%")
                # values.append(dict_data[lab])
                i_list.append(i)

        colors = ['lightsalmon', 'lightskyblue', 'goldenrod', 'mediumpurple', 'y', 'yellowgreen', 'violet',
                  'lightskyblue']
        colors1 = []
        for i in i_list:
            colors1.append(colors[i])
        # print()
        sum_val = 0
        for v in values:
            sum_val = sum_val + v
        if sum_val == 0:
            return 
        else:
            axs[index].pie(values, colors=colors1, startangle=140)
        # axs[index].pie(values, colors=colors, startangle=140)
        axs[index].set_title(cate)
        axs[index].legend(labels_with_pct, loc="best", bbox_to_anchor=(0.90, 0.00))

    plt.tight_layout()
    plt.savefig(f"{name}_TEENA_TE-derived-peak_annotation_pieplot.svg", format='svg')
    plt.close()

# if __name__ == '__main__':
#     draw_pie(pd.read_excel('result2_derived_peak.xlsx'), 'test')

