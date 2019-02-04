#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab
import os
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import argparse
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker




def main():
    args = get_args()

    chrom, all_pos, all_cov = read_samtools_depth_file(args.all_cov_file)

    chrom2, uniq_pos, uniq_cov = read_samtools_depth_file(args.uniq_cov_file)

    if chrom != chrom2:
        sys.exit("Error: different reference contig names for samtools depth files... exiting ...\n")

    plot = plot_coverage(chrom, all_pos, all_cov, uniq_pos, uniq_cov, args.title, args.plot_height, args.plot_width, args.normalize_cov, args.add_hline, args.output_file )
    
    plot.savefig(args.output_file, bbox_inches="tight")



def get_args():
    parser = argparse.ArgumentParser(description="Script to plot normalized samtools depth output")

    ## required ##
    parser.add_argument("-a", "--all_cov_file", type=str, help="File of samtools depth output for all coverage (MAPQ >= 0)", required=True)
    parser.add_argument("-u", "--uniq_cov_file", type=str, help="File of samtools depth output for unique coverage (MAPQ >= 1)", required=True)
    parser.add_argument("-o", "--output_file", type=str, help="Output file name", required=True)

    ## optional ##
    parser.add_argument("-n", "--title", help="title of plot", type=str, required=False)
    parser.add_argument("-H", "--plot_height", type=int, help="Plot height", default=3, required=False)
    parser.add_argument("-W", "--plot_width", type=int, help="Plot width", default=10, required=False)
    parser.add_argument("--normalize_cov", type=float, help="Divides raw coverage by given value to plot relative coverage", required=False)
    parser.add_argument("--add_hline", type=float, help="Adds a horizontal line to plot at given y-axis value", required=False)

    args = parser.parse_args()

    return args

def read_samtools_depth_file(depth_file):
    chrom = ""
    x = []
    y = []
    
    with open(depth_file,"r") as depth:
        for line in depth:
            split_line = line.split("\t")
            if chrom == "":
                chrom = split_line[0]
            elif chrom != split_line[0]:
                sys.exit("Error: Samtools depth file has multiple reference contigs.... exiting...\n")
            
            x.append(int(split_line[1]))
            y.append(int(split_line[2]))

    return chrom, x, y


def plot_coverage(chrom, all_pos, all_cov, uniq_pos, uniq_cov, title, height, width, normalize_cov, add_hline, output_file ):

    if normalize_cov:
        all_cov = np.array(all_cov)
        uniq_cov = np.array(uniq_cov)

        all_cov = all_cov/normalize_cov
        uniq_cov = uniq_cov/normalize_cov

    fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=300)

    plot = matplotlib.pyplot

    ax = plot.subplot(1,1,1)

    ax.fill_between(all_pos, all_cov, np.zeros(len(all_cov)), color='grey', step="pre", alpha=0.15)

    ax.fill_between(uniq_pos, uniq_cov, np.zeros(len(uniq_cov)), color='darkgrey', step="pre", alpha=0.4)


    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='y', colors='grey')
    ax.tick_params(axis='x', labelsize=5, length=0, pad=5)
    ax.tick_params(axis='y', labelsize=5, length=0, pad=5)

    ax.set_xlim(min(all_pos),max(all_pos))
    ax.set_ylim(0,max(all_cov))
    plot.xticks([min(all_pos), np.median(all_pos), max(all_pos)])
    plot.yticks([0, int(max(all_cov)/2), max(all_cov)])

    ax.set_xlabel('Position on '+ chrom, fontsize=8)
    

    if title:
        ax.set_title(" "*15 + title, fontsize=10, loc='left')


    ## legend

    labels = []
    elements = []

    elements += [mpatches.Patch(color='darkgrey',alpha=.4)]
    labels.append("Coverage (MAPQ > 0)")
    elements += [mpatches.Patch(color='grey',alpha=.15)]
    labels.append("Coverage (MAPQ = 0)")

    if add_hline is not None:
        ax.set_ylabel('Normalized Coverage', fontsize=8)
        ax.axhline(y=add_hline, color='black',alpha=0.5)
        elements += [matplotlib.pyplot.Line2D([0,0],[0,1],color='black', alpha=0.5)]
        labels.append("Avg. Norm. Cov. ("+format(add_hline, '.2f')+")")
    else:
        ax.set_ylabel('Coverage', fontsize=8)

    plot.legend( elements , labels, loc=1, fontsize = 6, frameon=False)


    return plot

if __name__ == '__main__':
    main()