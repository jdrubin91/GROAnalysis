__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import pybedtools
from pybedtools import BedTool
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import math

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def make_chromsize_dict(chromsizes):
    d = dict()
    with open(chromsizes) as F:
        line = line.strip().split()
        if len(line[0]) < 6 and line[0] != 'chrM':
            d[line[0]] = [int(line[1]),[]]

    return d

def run(A2N,ACN,chromsizes,figuredir):
    a = BedTool(A2N)
    b = BedTool(ACN)

    print "intersecting..."
    counts1 = (a+b).map(a,c='4',o='sum',null="0")
    counts2 = (a+b).map(b,c='4',o='sum',null="0")

    print counts1[:10]
    print counts2[:10]

    print "done\ncalculating fold change..."
    newbed = list()
    for i in range(len(counts1)):
        int1 = counts1[i]
        int2 = counts2[i]
        if int1[-1] != "0":
            val1 = float(int1[-1])
        else:
            val1 = 0.001
        if int2[-1] != "0":
            val2 = float(int2[-1])
        else:
            val2 = 0.001
        try:
            newbed.append(int1[:-1].append(math.log(val2/val1)))
        except:
            newbed.append(int1[:-1].append(0))

    print "done\ncreating and saving bedtool..."
    newbedtool = BedTool(newbed).saveas(parent_dir(A2N) + '/ACN_A2N_fold_change.BedGraph')

    print "done"



    # d = make_chromsize_dict(chromsizes)
    # print len(d)

    # for i in range(len(counts1)):
    #     chrom = counts1[i][0]
    #     if chrom in d:
    #         d[chrom]


    

    F = plt.figure()
    



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'

    A2N = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/A2N_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    ACN = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/ACN_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'

    chromsizes = '/scratch/Users/joru1876/hg19_reference_files/hg19.chrom.sizes.txt'

    run(A2N,ACN,chromsizes,figuredir)