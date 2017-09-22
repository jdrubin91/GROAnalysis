__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import math
import numpy as np

def format_bp(bp):
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

def run(bam1,bam2,bam3,bam4,chip,filedir,figuredir):
    os.system("sort -k1,1 -k2,2n " + chip + " > " + chip + ".sorted.bed")
    os.system("bedtools multicov -bams " + bam1 + " " + bam2 + " " + bam3 + " " + bam4 + " -bed " + chip + ".sorted.bed >" + filedir + "Boxplot_ChIPPeak_GROcoverage.counts.bed")

    boxplot = list()
    with open(filedir + "Boxplot_ChIPPeak_GROcoverage.counts.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            try:
                boxplot.append((math.log(float(line[-2])/float(line[-1]),2)+math.log(float(line[-4])/float(line[-3]),2))/2.0)
            except:
                pass

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('Log2 Fold-Change over ATF3 Peaks rep1')
    ax.set_ylabel('Log2 Fold-Change t45 DMSO/CA')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.axhline(0, color='black', alpha=0.5)
    bp = ax.boxplot(boxplot, patch_artist=True)
    format_bp(bp)
    F.savefig(figuredir + 'FoldChange_t45_rep1.png', dpi=1200)




if __name__ == "__main__":
    filedir = "/Users/joru1876/scratch_backup/GROAnalysis/files/"
    figuredir = "/Users/joru1876/scratch_backup/GROAnalysis/figures/"
    chip = filedir + 'ATF3_ChIP_rep1and2.bed'
    bamfolder1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/'
    bam1 = bamfolder1 + 'J52_trimmed.flip.fastq.bowtie2.sorted.bam'
    bam2 = bamfolder1 + 'J62_trimmed.flip.fastq.bowtie2.sorted.bam'
    bamfolder2 = '/projects/dowellLab/Taatjes/170825_NB501447_0152_fastq/Demux/cat/trimmed/flipped/bowtie2/sortedbam/'
    bam3 = bamfolder2 + 'J5D451_GTCCGC_S3_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    bam4 = bamfolder2 + 'J6C451_GTGAAA_S4_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    run(bam1,bam2,bam3,bam4,chip,filedir,figuredir)
