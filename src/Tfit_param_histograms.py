__author__ = 'Jonathan Rubin'

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import norm
import numpy as np
import math

#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

#Home directory
homedir = os.path.dirname(os.path.realpath(__file__))

#File directory
filedir = parent_dir(homedir) + '/files'

#Figure directory
figuredir = parent_dir(homedir) + '/figures/'

def run(folder):
    names = ['mu_k', 'sigma_k lambda_k', 'pi_k', 'fp_k', 'w_[p,k]', 'w_[f,k]', 'w_[r,k]', 'b_[f,k]', 'a_[r,k]']
    values = [[] for i in range(len(names))]
    for file1 in os.listdir(folder):
        if 'K_models_MLE.tsv' in file1:
            with open(folder + file1) as F:
                for line in F:
                    if '#' not in line[0]:
                        if '>' in line[0]:
                            i = 0
                        elif i == 2:
                            line = line.strip().split()[1:]
                            w = line[4].split(',')
                            for k in range(len(line)):
                                if k == 5:
                                    for l in range(len(w)):
                                        values[k+l].append(float(w[l]))
                                elif k > 5:
                                    values[k+2].append(float(line[k]))
                                else:
                                    values[k].append(float(line[k]))
                        else:
                            i+=1

        length = len(names)
        subplotmatrix = int(length)
        F = plt.figure()
        F.suptitle(file1, fontsize=14)
        for i in range(length):
            print i,len(values)
            ax = F.add_subplot(subplotmatrix,subplotmatrix,i)
            plt.hist(values[i],bins=100)
            ax.set_title(names[i])

        plt.savefig(figuredir + file1 + '.png')


    

if __name__ == "__main__":
    folder = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/'
    run(folder)