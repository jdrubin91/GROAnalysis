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
import numpy as np

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def convert_meme_to_bed(memefile):
    bed = list()
    boolean = False
    with open(memefile) as F:
        for line in F:
            if 'Motif 1 sites sorted by position p-value' in line:
                boolean = True
            if boolean:
                if 'chr' in line:
                    line = line.strip().split()
                    chrom = line[0].split(':')[0]
                    start = line[0].split(':')[1].split('-')[0]
                    stop = line[0].split(':')[1].split('-')[1]
                    bed.append([chrom,start,stop])
            if 'Motif 1 block diagrams' in line:
                break

    return bed

def convert_deseq_to_bed(deseqfile):
    bed = list()
    with open(deseqfile) as F:
        for line in F:
            line.strip().split()
            if 'id' not in line[0] and 'NA' not in line[0]:
                print line
                item = line[1].split(';')[-1]
                print item
                chrom = item.split(':')[0]
                start = item.split(':')[1].split('-')[0]
                stop = item.split(':')[1].split('-')[1].split('_')[0]
                bed.append([chrom,start,stop])

    return bed


def run(memefile,deseqfile,bg1,bg2,figuredir):
    meme = BedTool(memefile)
    deseq = BedTool(deseqfile)
    meme = deseq.intersect(meme,wa=True)
    a = BedTool(bg1)
    b = BedTool(bg2)


    acm = meme.map(a,c=4,o="sum")
    acd = (deseq-meme).map(a,c=4,o="sum")

    bcm = meme.map(b,c=4,o="sum")
    bcd = (deseq-meme).map(b,c=4,o="sum")

    mf = [float(m)/float(n) for m,n in zip(acm,bcm) if m != '.' and n != '.']
    df = [float(m)/float(n) for m,n in zip(acd,bcd) if m != '.' and n != '.']


    F = plt.figure()
    ax = F.add_subplot(121)
    ax.set_title('AT-rich genes')
    ax.set_ylabel('Count')
    ax.set_xlabel('Fold Change')
    ax.hist(mf,bins = 100)

    ax2 = F.add_subplot(122)
    ax2.set_title('non AT-rich genes')
    ax2.set_ylabel('Count')
    ax2.set_xlabel('Fold Change')
    ax2.hist(df,bins=100)
    plt.savefig(figuredir + 'AT_rich_gene_expression_changes.png',dpi=1200)




if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'
    memedir = parent_dir(homedir) + '/MEME/'
    memefile = memedir + '/Genes/meme.txt'

    deseqfile = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/DE-Seq/A2N_ACN.genes.bed.count.bed.A2NACNnascent.resSig_pvalue.txt'

    bg1 = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/A2N_trimmed.flip.fastq.bowtie2.sorted.reflected.BedGraph'
    bg2 = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/ACN_trimmed.flip.fastq.bowtie2.sorted.reflected.BedGraph'

    memefile = convert_meme_to_bed(memefile)
    deseqfile = convert_deseq_to_bed(deseqfile)
    run(memefile,deseqfile,bg1,bg2,figuredir)
