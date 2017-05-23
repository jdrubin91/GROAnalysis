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

def run(deseqfile,cond1,cond2,figuredir):
    x = list()
    y = list()
    sigx = list()
    sigy = list()
    with open(deseqfile) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            if 'NA' not in line[-1]:
                p = float(line[-2])
                x.append(math.log(float(line[2])))
                y.append(float(line[-3]))
                if p < 0.1:
                    sigx.append(math.log(float(line[2])))
                    sigy.append(float(line[-3]))

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.scatter(x,y,edgecolor='')
    ax.scatter(sigx,sigy,c='r',edgecolor='')
    ax.set_title('Gene Transcription')
    ax.set_ylabel('Log2 Fold Change ' + cond2 + '/' + cond1)
    ax.set_xlabel('Log10 Mean Transcription')
    ax.set_xlim([0,17])
    plt.savefig(figuredir + deseqfile.split('/')[-1] + '.png', dpi=1200)


if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'

    cond1 = 'A2N'
    cond2 = 'ACN'

    deseqfile = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/DE-Seq/'+cond1+'_'+cond2+'.genes.bed.count.bed.'+cond1+cond2+'nascent.res.txt'

    run(deseqfile,cond1,cond2,figuredir)

