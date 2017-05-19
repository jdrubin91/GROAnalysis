__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import sys
import pybedtools
from pybedtools import BedTool
import math

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def get_histone_bed(histones,genes):
    names = list()
    with open(histones) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            names.append(line[-5])

    bed = list()
    with open(genes) as F:
        for line in F:
            line = line.strip().split()
            geneName = line[3].split(';')[0]
            if geneName in names:
                chrom,start,stop = line[:3]
                bed.append([chrom,start,stop])

    return bed

def run(bg1,bg2,genes,histones,figuredir):
    a = BedTool(bg1)
    b = BedTool(bg2)

    bed = BedTool(get_histone_bed(histones,genes))

    m = bed.map(a,c=4,o="sum")
    n = bed.map(b,c=4,o="sum")

    x = list()
    for item in m:
        try:
            x.append(math.log(item[-1]))
        except:
            x.append(0)

    print x

    y = list()
    for item in n:
        try:
            y.append(math.log(item[-1]))
        except:
            y.append(0)

    print len(x),len(y)

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('Histone Gene Transcription')
    ax.set_ylabel('A2780cis Log10(FPKM)')
    ax.set_xlabel('A2780 Log10(FPKM)')
    ax.scatter(x,y)
    plt.savefig(figuredir + 'Histone_genes_fpkm.png',dpi=1200)



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'
    scriptdir = parent_dir(homedir) + '/scripts/'
    outdir = parent_dir(homedir) + '/MEME/'
    genes = filedir + 'refGene.sorted.bed'
    histones = filedir + 'histone_names.txt'

    bg1 = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/A2N_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    bg2 = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/ACN_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'

    run(bg1,bg2,genes,histones,figuredir)