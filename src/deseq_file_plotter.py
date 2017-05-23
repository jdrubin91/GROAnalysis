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

def get_histone_bed(histones):
    names = list()
    with open(histones) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            names.append(line[-5])

    return names

def get_cell_cycle_names(cell_cycle):
    names = list()
    with open(cell_cycle) as F:
        F.readline()
        F.readline()
        for line in F:
            names.append(line.strip())

    return names

def run(deseqfile,cond1,cond2,figuredir,histone_names,cell_cycle_names,Sphase_names,DNArepair_names):
    x = list()
    y = list()
    sigx = list()
    sigy = list()
    hisx = list()
    hisy = list()
    ccx = list()
    ccy = list()
    with open(deseqfile) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            if 'NA' not in line[-1]:
                gene = line[1].split(';')[0][1:]
                geneName = line[1].split(';')[1]
                p = float(line[-2])
                x.append(math.log(float(line[2])))
                y.append(float(line[-3]))
                if p < 0.01:
                    sigx.append(math.log(float(line[2])))
                    sigy.append(float(line[-3]))
                # if geneName in cell_cycle_names:
                #     ccx.append(math.log(float(line[2])))
                #     ccy.append(float(line[-3]))
                # if gene in histone_names:
                #     hisx.append(math.log(float(line[2])))
                #     hisy.append(float(line[-3]))
                if geneName in Sphase_names:
                    ccx.append(math.log(float(line[2])))
                    ccy.append(float(line[-3]))
                # if geneName in DNArepair_names:
                #     ccx.append(math.log(float(line[2])))
                #     ccy.append(float(line[-3]))


    name1 = 'A2780'
    if cond1[1] == 'C':
        name1 = name1 + 'cis'
    if cond1[2] == 'C':
        name1 = name1 + ' +Cis'
    if cond1[2] == 'D':
        name1 = name1 + ' +Dox'

    name2 = 'A2780'
    if cond2[1] == 'C':
        name2 = name2 + 'cis'
    if cond2[2] == 'C':
        name2 = name2 + ' +Cis'
    if cond2[2] == 'D':
        name2 = name2 + ' +Dox'

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.scatter(x,y,edgecolor='')
    ax.scatter(sigx,sigy,c='r',edgecolor='')
    # ax.scatter(hisx,hisy,c='g')
    ax.scatter(ccx,ccy,c='y')
    ax.set_title('Gene Transcription ' + name2 + ' vs. ' + name1)
    ax.set_ylabel('Log2 Fold Change ' + name2 + '/' + name1)
    ax.set_xlabel('Log10 Mean Transcription')
    ax.set_xlim([min(x),max(x)])
    # plt.savefig(figuredir + deseqfile.split('/')[-1] + '_histones.png', dpi=1200)
    # plt.savefig(figuredir + deseqfile.split('/')[-1] + '_cell_cycle.png', dpi=1200)
    plt.savefig(figuredir + deseqfile.split('/')[-1] + '_S_phase.png', dpi=1200)
    # plt.savefig(figuredir + deseqfile.split('/')[-1] + '.png', dpi=1200)
    # plt.savefig(figuredir + deseqfile.split('/')[-1] + '_DNA_repair.png', dpi=1200)



if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'
    genes = filedir + 'refGene.sorted.bed'
    histones = filedir + 'histone_names.txt'
    cell_cycle = filedir + 'cell_cycle_genes.txt'
    S_phase = filedir + 'Sphase_genes.txt'
    DNA_repair = filedir + 'DNA_repair.txt'

    cond1 = 'A2N'
    cond2 = 'ACN'

    deseqfile = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/DE-Seq/'+cond1+'_'+cond2+'.genes.bed.count.bed.'+cond1+cond2+'nascent.res.txt'

    histone_names = get_histone_bed(histones)
    cell_cycle_names = get_cell_cycle_names(cell_cycle)
    Sphase_names = get_cell_cycle_names(S_phase)
    DNArepair_names = get_cell_cycle_names(DNA_repair)

    run(deseqfile,cond1,cond2,figuredir,histone_names,cell_cycle_names,Sphase_names,DNArepair_names)

