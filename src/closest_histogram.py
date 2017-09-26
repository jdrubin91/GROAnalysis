__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import math

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(file1,file2,filedir,figuredir):
    os.system("sort -k1,1 -k2,2n " + file1 + " > " + file1 + ".sorted.bed")
    os.system("sort -k1,1 -k2,2n " + file2 + " > " + file2 + ".sorted.bed")
    os.system("bedtools closest -d -a " + file1 + ".sorted.bed" + " -b " + file2 + ".sorted.bed" + " > " + filedir + "closestbed_distances.bed")
    print "bedtools closest -d -a " + file1 + " -b " + file2 + " > " + filedir + "closestbed_distances.bed"
    hist = list()
    with open(filedir + "closestbed_distances.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            if float(line[-1]) < 50000:
                hist.append(float(line[-1]))
    F = plt.figure()
    plt.hist(hist,bins=100)
    plt.xlim([0,50000])
    plt.title("Distance between eRNAs DMSO rep1 to rep2 t45 50kb")
    plt.xlabel("Distance (bp)")
    plt.savefig(figuredir + "closest_histogram_DMSOt45_50.png")


if __name__ == "__main__":
    homedir = os.path.dirname(os.path.realpath(__file__))
    figuredir = parent_dir(homedir) + '/figures/'
    filedir =  '../files/'

    folder1 = '../../SerumResponseCA_REP1GROSEQ/Tfit_rep1/'
    folder2 = '../../SerumResponseCA_REP1GROSEQ/Tfit_rep2/'
    filelist = [folder1+'foot_print_testing-7_bidir_predictions.bed',
                folder1+'foot_print_testing-8_bidir_predictions.bed',
                folder1+'foot_print_testing-9_bidir_predictions.bed',
                folder1+'foot_print_testing-10_bidir_predictions.bed',
                folder1+'foot_print_testing-11_bidir_predictions.bed',
                folder1+'foot_print_testing-12_bidir_predictions.bed',
                folder2+'foot_print_testing-7_bidir_predictions.bed',
                folder2+'foot_print_testing-8_bidir_predictions.bed',
                folder2+'foot_print_testing-9_bidir_predictions.bed',
                folder2+'foot_print_testing-10_bidir_predictions.bed',
                folder2+'foot_print_testing-11_bidir_predictions.bed',
                folder2+'foot_print_testing-12_bidir_predictions.bed']

    file1 = filelist[4]
    file2 = filelist[10]
    run(file1,file2,filedir,figuredir)