__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import math

file1 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_DMSO-1_K_models_MLE.tsv'
file2 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_CA-1_K_models_MLE.tsv'
savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'
index = 6
cut = 50
cut1 = 50000
cut2 = -50000

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def run(file1,file2):
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
                    p = line.strip().split()
                    for param in p:
                        if '~' in param:
                            d1[gene].append(float(param.split(',')[1]))
                        elif ',' in param:
                            for moreparam in param.split(','):
                                d1[gene].append(float(moreparam))
                        else:
                            d1[gene].append(float(param))
                    
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d2[gene] = [float(fwd),float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()
                    for param in p:
                        if '~' in param:
                            d2[gene].append(float(param.split(',')[1]))
                        elif ',' in param:
                            for moreparam in param.split(','):
                                d2[gene].append(float(moreparam))
                        else:
                            d2[gene].append(float(param))
    
    
    for i in range(2,len(p[2:])+2):
        X = list()
        Y = list()
        x = list()
        y = list()
        for key in d1:
            if key in d2:
                if d1[key][0] > cut or d1[key][1] > cut and d2[key][0] > cut or d2[key][1] > cut:
                    if not (math.isinf(d1[key][i]) or math.isinf(d2[key][i]) or math.isnan(d1[key][i]) or math.isnan(d2[key][i])):
                        #if cut1 > d1[key][i] > cut2 and cut1 > d2[key][i] > cut2:
                        x.append(d1[key][i])
                        y.append(d2[key][i])
        #            if d1[key][2] != 0:
        #                if d2[key][2]-d1[key][2] > .25:
        #                    Y.append(key)
        #                else:
        #                    X.append(d2[key][2]-d1[key][2])
        #                
        #print "max: " + str(max(X))
        #print "min: " + str(min(X))
        #print "length: " + str(len(X))
        #print "avg: " + str(sum(X)/len(X))
        #print Y
        
        #plt.hist(X,50)
        print len(x)
        print len(y)
        F = plt.figure()
        plt.scatter(x,y,alpha=0.1)
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        plt.scatter(x,y,c=z,edgecolor="",s=14)
        plt.savefig(savedir + 'tsv_fig' + str(i) + '.png')
    
    return "done"
    
def run2(file1,file2):
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
    
    
    X = list()
    Y = list()
    Z = list()
    for key in d1:
        if key in d2:
            if d1[key][0] > cut or d1[key][1] > cut and d2[key][0] > cut or d2[key][1] > cut:
                if d2[key][2]-d1[key][2] > .25:
                    Y.append((d2[key][2]-d1[key][2],key))
                if d2[key][2]-d1[key][2] < -.25:
                    Z.append((d2[key][2]-d1[key][2],key))
                X.append(d2[key][2]-d1[key][2])
                    
    print "max: " + str(max(X))
    print "min: " + str(min(X))
    print "length: " + str(len(X))
    print "avg: " + str(sum(X)/len(X))
    print "Y: ",Y
    print "Z: ",Z
    print "len(Y): ",len(Y)
    print "len(Z): ",len(Z)
    F = plt.figure()        
    plt.hist(X,50)
    plt.savefig(savedir + 'tsv_fig.png')
    
    return "done"
    

    
if __name__ == "__main__":
    run2(file1,file2)