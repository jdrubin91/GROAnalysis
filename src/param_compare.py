__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import intervals
from scipy.stats import gaussian_kde
import numpy as np

file1 = '/scratch/Shares/dowell/ENCODE/Allen2014_DMSO2_3-19_divergent_classifications.bed'
file2 = '/scratch/Shares/dowell/ENCODE/Allen2014_Nutlin2_3-3_divergent_classifications.bed'
savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'

#Runs interval search over all bed sites in both files, recovers parameters in 
#bidirectional model
def run(file1,file2):

    d1 = dict()
    d2 = dict()
    A = list()
    B = list()
    
    with open(file1) as F1:
        for line in F1:
            if not '#' in line[0]:
                chrom, start, stop, param = line.strip().split()
                param = param.split('|')[1].split(',')
                d1[chrom + ':' + start + '-' + stop] = param
                A.append((int(start),int(stop),'A',chrom))
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                chrom, start, stop, param = line.strip().split()
                param = param.split('|')[1].split(',')
                d2[chrom + ':' + start + '-' + stop] = param
                B.append((int(start),int(stop),'B',chrom))
                
    ST = intervals.comparison((A,B))
    OVERLAPS_0_1 = ST.find_overlaps(0,1)
    print "Overlap Instances: " + str(len(OVERLAPS_0_1))
    for index in range(len(param)):
        x = list()
        y = list()
        for O in OVERLAPS_0_1:
            if not len(O.overlaps.keys()) > 2:
                comparison = list()
                for interval_original in O.overlaps:
                    comparison.append(interval_original.INFO[1])
                if comparison[0] == comparison[1]:
                    for interval_original in O.overlaps:
                        if 'A' in interval_original.INFO:
                            x.append(float(d1[interval_original.INFO[1] + ':' + str(interval_original.start) + '-' + str(interval_original.stop)][index]))
                        elif 'B' in interval_original.INFO:
                            y.append(float(d2[interval_original.INFO[1] + ':' + str(interval_original.start) + '-' + str(interval_original.stop)][index]))
                            
        F = plt.figure()
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        plt.scatter(x,y,c=z,edgecolor="",s=14)
        plt.savefig(savedir + 'figure' + str(index) + '.png')
    
    return x,y
#Creates a scatter plot of x and y values
def plot(x,y):
    F = plt.figure()
    plt.scatter(x,y)
    plt.savefig(savedir)

if __name__ == "__main__":
    x,y = run(file1,file2)
    plot(x,y)