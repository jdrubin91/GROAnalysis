__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
from pybedtools import BedTool
import math

def column_add(file1,file2,name):
    outfile = open(name,'w')
    with open(file1) as F:
        with open(file2) as F2:
            for line in F:
                line2 = F2.readline()
                line = line.strip() + '\t' + line2.strip().split()[-1] + '\n'

def run(bedgraph1,bedgraph2,bed,figdir,filedir):
    b1name = bedgraph1.split('/')[-1]
    b2name = bedgraph2.split('/')[-1]
    a = BedTool(bedgraph1)
    b = BedTool(bedgraph2)
    s = BedTool(bed).sort()
    s.map(a,c="4",o="sum").saveas(filedir + b1name + '_gene_counts.bed')
    s.map(b,c="4",o="sum").saveas(filedir + b2name + '_gene_counts.bed')
    column_add(filedir + b1name + '_gene_counts.bed',filedir + b2name + '_gene_counts.bed', filedir + 'both_gene_counts.bed')
    d = dict()
    with open(filedir + 'both_gene_counts.bed') as F:
        for line in F:
            line = line.strip().split()
            try:
                d[line[3]] = math.log(float(line[-2])/float(line[-1]),10)
            except:
                d[line[3]] = 0
    F = plt.figure()
    ax = F.add_subplot(111)
    print d.values()
    ax.hist(d.values(),bins=30)
    plt.savefig(figdir + 'SE_analysis.png')


if __name__ == "__main__":
    gene_list = ['']
    bed = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/refGene.sorted.bed'
    figdir = '/scratch/Users/joru1876/GROAnalysis/figures/'
    filedir = '/scratch/Users/joru1876/GROAnalysis/files/'
    bedgraph1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J52_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    bedgraph2 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J62_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    run (bedgraph1,bedgraph2,bed,figdir,filedir)