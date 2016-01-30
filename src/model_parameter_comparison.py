__author = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np

file1 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_DMSO-1_K_models_MLE.tsv'
file2 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_CA-1_K_models_MLE.tsv'
file3 = '/scratch/Users/joru1876/GROAnalysis/files/Master.bed'
savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'

def run(DMSOdiv,CAdiv,Master):
    d1 = dict()
    d2 = dict()
    with open(DMSOdiv) as F1:
        for line in F1:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d1[gene] = [float(fwd), float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()
                    for param in p:
                        if '~' in param:
                            d1[gene].append(float(param.split(',')[1]))
                        elif ',' in param:
                            for moreparam in param.split(','):
                                d1[gene].append(float(moreparam))
                        else:
                            d1[gene].append(float(param))
                            
                            
def run2(file1,file2,file3):
    index = 6
    d1 = dict()
    d2 = dict()
    with open(file1) as F1:
        for line in F1:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d1[gene] = [float(fwd), float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()[index].split(',')[0]
                    d1[gene].append(float(p))
                    
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d2[gene] = [float(fwd),float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()[index].split(',')[0]
                    d2[gene].append(float(p))
    cut = 200
    d3 = dict()
    with open(file3) as F3:
        F3.readline()
        for line in F3:
            line = line.strip().split()
            gene = line[0]
            DMSOgenes, DMSOTSS = line[6:8]
            CAgenes,CATSS = line[9:11]
            DMSOgenes = float(DMSOgenes)
            DMSOTSS = float(DMSOTSS)
            CAgenes = float(CAgenes)
            CATSS = float(CATSS)
            if CAgenes-CATSS != 0 and DMSOgenes-DMSOTSS != 0 and CAgenes > cut and DMSOgenes > cut:
                TR = (CATSS/(CAgenes-CATSS))-(DMSOTSS/(DMSOgenes-DMSOTSS))
                d3[gene] = TR
    x = list()
    y = list()
    for gene in d1:
        if gene in d2:
            for key in d3:
                if gene in key:
                    if d1[gene][0] > cut or d1[gene][1] > cut and d2[gene][0] > cut or d2[gene][1] > cut:
                        x.append(d2[gene][2]-d1[gene][2])
                        y.append(d3[key])
    F = plt.figure()
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    plt.scatter(x,y,c=z,edgecolor="",s=14)
    plt.ylim([-1,1])
    plt.xlim([-1,1])
    plt.savefig(savedir + '/model_parameter_comparison.png')
    
if __name__ == "__main__":
    run2(file1,file2,file3)