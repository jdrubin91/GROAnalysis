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

def parse_chipdir(chipdir):
    d = dict()
    file1 = chipdir + 'metadata.tsv'
    names = list()
    with open(file1) as F:
        header = F.readline().strip('\n').split('\t')
        for line in F:
                line = line.strip('\n').split('\t')
                genome = line[-8]
                filename = line[0]
                fileformat = line[2]
                filetype = line[1]
                name = line[16].split('-')[0]
                rep = line[29]
                techrep = line[30]
                if fileformat == 'peaks' and genome == 'hg19' and 'POLR' not in name and 'CTCF' not in name and 'bed' in filetype:
                    names.append(name)
                    d[filename] = [name,rep,techrep,genome,filetype,fileformat]
    names = list(set(names))
    print names
    print d
    templist = list()
    filelist = list()
    for name in names:
        for filename in d:
            if d[filename][0] == name:
                templist.append(chipdir + filename + '.bed')
        os.system("bedtools intersect -a " + templist[0] + " -b " + " ".join(templist[1:]) + " > " + chipdir + name + ".all_intersect.bed")
        filelist.append(chipdir + name + ".all_intersect.bed")

    return filelist,names



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

    return boxplot

def plot(boxplot,names,figuredir):
    means = [np.mean(x) for x in boxplot]
    indices = [i[0] for i in sorted(enumerate(means), key=lambda x:x[1], reverse=True)]
    boxplot = [boxplot[i] for i in indices]
    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('Log2 Fold-Change over HCT116 Encode ChIP Peaks')
    ax.set_ylabel('Log2 Fold-Change t45 DMSO/CA')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.axhline(0, color='black', alpha=0.5)
    bp = ax.boxplot(boxplot, positions=np.arange(len(boxplot)),patch_artist=True)
    plt.xticks(np.arange(len(boxplot)),names)
    format_bp(bp)
    F.savefig(figuredir + 'FoldChange_HCT116_ChIP_t45.png', dpi=1200)




if __name__ == "__main__":
    filedir = "/Users/joru1876/scratch_backup/GROAnalysis/files/"
    figuredir = "/Users/joru1876/scratch_backup/GROAnalysis/figures/"
    # chip = filedir + 'ATF3_ChIP_rep1and2.bed'
    # chip = filedir + 'JUND_ChIP.bed'
    # chip = filedir + 'SRF_ChIP_rep1and2.bed'
    bamfolder1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/'
    bam1 = bamfolder1 + 'J52_trimmed.flip.fastq.bowtie2.sorted.bam'
    bam2 = bamfolder1 + 'J62_trimmed.flip.fastq.bowtie2.sorted.bam'
    bamfolder2 = '/projects/dowellLab/Taatjes/170825_NB501447_0152_fastq/Demux/cat/trimmed/flipped/bowtie2/sortedbam/'
    bam3 = bamfolder2 + 'J5D451_GTCCGC_S3_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    bam4 = bamfolder2 + 'J6C451_GTGAAA_S4_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    chipdir = '/Users/joru1876/scratch_backup/HCT116_ChIP/'
    filelist,names = parse_chipdir(chipdir)
    boxplot = list()
    for chip in filelist:
        boxplot.append(run(bam1,bam2,bam3,bam4,chip,filedir,figuredir))
    plot(boxplot,names,figuredir)
