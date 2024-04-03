#!/usr/bin/env python
#-*- coding:utf-8 -*-
# Author: Yuzhuo Li  (DX120230210@stu.yzu.edu.cn)
# Datetime:2023-11-23
# Software: TEENA: determine the enrichment of TE family within given genomic intervals using Fisher's Exact Test
# Version: 1.1
# Update 2023-12-04: The width of the svg figures is changed from 16 to 18; the x-label and y-label for the bar plots are changed.
# Update 2024-03-29: Added color bars and labels filled with foldenrich values
# Update 2024-03-31: Changed the label color and location

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import seaborn as sns

def read_fet(filename):
    d = pd.read_excel(filename)
    d['Rname'] = d['TE_name']
    d['Rclass'] = d['TE_class']
    pval = d['p_adj']

    df = pd.DataFrame({
        'class': d['Rclass'],
        'name': d['Rname'],
        'obs_hits': d['a1_b1'],
        'obs_freq': d['a1_b1'] / (d['a1_b1'] + d['a1_b0']),
        'exp_hits': d['a0_b1'] * ((d['a1_b1'] + d['a1_b0']) / (d['a0_b1'] + d['a0_b0'])),
        'exp_freq': d['a0_b1'] / (d['a0_b1'] + d['a0_b0']),
        'p': pval,
        'p_adj': pval,
        'ratio': d['fold_enrich']
    })
    df['fold'] = (df['obs_hits'] + 0.001) / (df['exp_hits'] + 0.001)
    df['log2fold'] = np.log2(df['fold'])
    df['minus_log10_padj'] = -np.log10(df['p_adj'])
    return df


def draw_scatter(df_enh_all, prefix):
    four_category = ['DNA', 'LINE', 'SINE', 'LTR']
    # fig, axs = plt.subplots(1, 4, figsize=(16, 4))
    # MAS, 231204
    fig, axs = plt.subplots(1, 4, figsize=(18, 4))

    # Find the maximum value form exp_freq and obs_freq
    max_exp_freq = df_enh_all['exp_freq'].max() * 1.1
    max_obs_freq = df_enh_all['obs_freq'].max() * 1.1
    max_xy = max(max_exp_freq, max_obs_freq)
    #max_xy = 0.04

    for category in four_category:
        index = four_category.index(category)
        df_enh_filter = df_enh_all[df_enh_all['class'] == category]
        sns.scatterplot(ax=axs[index], data=df_enh_filter, x='exp_freq', y='obs_freq',
                        palette='coolwarm', s=50)

        # Set the limits of the x and y axes using the calculated maximum value
        #axs[index].set_xlim(0, max_exp_freq)
        #axs[index].set_ylim(0, max_obs_freq)
        axs[index].set_xlim(0, max_xy)
        axs[index].set_ylim(0, max_xy)

        axs[index].set_xlabel("Expected frequency")
        axs[index].set_ylabel("Observed frequency")
        axs[index].axline((0, 0), slope=1, color="grey", lw=1)
        axs[index].set_title(category)

    plt.tight_layout()
    plt.savefig(f"{prefix}_TEENA_scatterplot.svg", format='svg')
    plt.close()


def draw_bar(df_enh_all, prefix):
    four_category = ['DNA', 'LINE', 'SINE', 'LTR']
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2))

    f_list = []
    for i in range(4):
        df = df_enh_all[df_enh_all['class'] == four_category[i]].sort_values(
            by='minus_log10_padj', ascending=False).head(10)
        f_list.append(df['ratio'].values.tolist())
    vmin = min(min(row) for row in f_list)
    vmax = max(max(row) for row in f_list)

    for category in four_category:
        index = four_category.index(category)
        df_enh_filter = df_enh_all[df_enh_all['class'] == category].sort_values(
            by='minus_log10_padj', ascending=False).head(10)

        names = df_enh_filter['name'].to_list()
        padjs = df_enh_filter['minus_log10_padj'].to_list()

        cmap = plt.cm.Oranges
        # Set the normalized range
        norm = Normalize(vmin=vmin, vmax=vmax)
        # Get the normalized color array
        colors = cmap(norm(df_enh_filter['ratio']))

        axs[index].barh(np.arange(len(names)), padjs,
                        align='center', color=colors, ec='grey', lw=.5)
        axs[index].set_yticks(np.arange(len(names)))
        axs[index].set_yticklabels(names)
        axs[index].set_title(category)
        axs[index].set_xlabel('-log10(P)')
        axs[index].set_ylabel('')
        axs[index].set_xlim([0, None])

        # Create a ScalarMappable object
        sm = cm.ScalarMappable(cmap='Oranges', norm=norm)
        sm.set_array([])  # Set an empty array or list that the range of color mapping is displayed correctly
        cbar = plt.colorbar(sm, ax=axs[index], ticklocation='bottom', shrink=0.95,
                            orientation='horizontal', location='bottom', pad=0.16)

        cbar.ax.xaxis.set_ticks_position('bottom')
        cbar.ax.tick_params(labelsize=8, length=1.5)
        # cbar.ax.set_title('FoldEnrich', fontsize=10, loc='bottom')
        # set the title and move it to the bottom
        # ax = plt.gca()
        cbar.ax.set_xlabel('FoldEnrich', fontsize=10, color='k')
        axs[index].invert_yaxis()  # Reverse the Y-axis
    plt.tight_layout()
    plt.savefig(f"{prefix}_TEENA_barplot.svg", format='svg')
    plt.close()


def draw_plot(file, prefix):
    df_enh_all = read_fet(file)
    if df_enh_all['p_adj'].min() < 0.05:
        draw_scatter(df_enh_all, prefix)
        draw_bar(df_enh_all, prefix)
    else:
        draw_scatter(df_enh_all, prefix)
        print('No significant TE identified')

#draw_plot('result2.xlsx', 'demo')

# df = read_fet('CD14_STAT1_ENA.xlsx')
# draw_bar(df, prefix=None)
