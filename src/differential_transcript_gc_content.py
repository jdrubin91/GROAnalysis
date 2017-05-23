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

def calculate_gc_content(sequence):
    total = float(len(sequence))
    gc = 0.0
    at = 0.0
    for char in sequence:
        if char == 'G' or char == 'C':
            gc += 1.0
        else:
            at += 1.0

    return gc/total

def get_tss(bedtool):
    newbedtool = list()
    for interval in bedtool:
        chrom,start,stop = interval[:3]
        strand = interval[-1]
        if strand == '+':
            newbedtool.append([chrom,str(int(start)-75),str(int(start)+75)])
        else:
            newbedtool.append([chrom,str(int(stop)-75),str(int(stop)+75)])

    return BedTool(newbedtool)


def split_deseq_file(deseqfile):
    up = list()
    down = list()
    with open(deseqfile) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            if 'NA' not in line[0]:
                interval = line[1].split(';')[-1]
                chrom = interval.split(':')[0]
                start = interval.split(':')[1].split('-')[0]
                stop = interval.split(':')[1].split('-')[0].split('_')[0]
                strand = interval.split('_')[1]
                if float(line[-3]) < 0:
                    down.append([chrom,start,stop,strand])
                else:
                    up.append([chrom,start,stop,strand])

    return up,down

def bulk_gc_content(fastafile):
    gc_content = list()
    for line in open(fastafile.seqfn):
        if '>' not in line:
            gc_content.append(calculate_gc_content(line))

    return gc_content


def split_joey_deseq_file(joeydeseqfile):
    return "something"

def run(hg19fasta,genes,deseqfile,figuredir):
    up,down = split_deseq_file(deseqfile)
    u = get_tss(BedTool(up)).sequence(fi=hg19fasta)
    d = get_tss(BedTool(down)).sequence(fi=hg19fasta)
    g = get_tss(BedTool(genes)).sequence(fi=hg19fasta)

    F = plt.figure()
    ax = F.add_subplot(111)
    # ax.hist(bulk_gc_content(u),alpha=0.5,color='red',bins=100)
    # ax.hist(bulk_gc_content(d),alpha=0.5,color='green',bins=100)
    # plt.savefig(figuredir + 'differential_transcription_gc_content.png',dpi=1200)
    ax.hist(bulk_gc_content(g),bins=100)





if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'
    hg19fasta = '/scratch/Users/joru1876/hg19_reference_files/hg19_all.fa'

    genes = filedir + 'refGene.sorted.bed'
    deseqfile = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/DE-Seq/A2N_ACN.genes.bed.count.bed.A2NACNnascent.resSig_pvalue.txt'

    run(hg19fasta,genes,deseqfile,figuredir)
