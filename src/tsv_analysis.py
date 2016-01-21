__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

file1 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_DMSO-1_K_models_MLE.tsv'
file2 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_CA-1_K_models_MLE.tsv'
savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'
index = 6

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
    for key in d1:
        if key in d2:
            if d1[key][0] > 10 or d1[key][1] > 10 and d2[key][0] > 10 or d2[key][1] > 10:
                if d2[key][2] != 0:
                    print d1[key][2]
                    print d2[key][2]
                    X.append(d1[key][2]/d2[key][2])
                
    F = plt.figure()
    plt.hist(X,bins=100)
    plt.savefig(savedir + 'tsv_fig.png')
    
    return "done"
    
if __name__ == "__main__":
    run(file1,file2)