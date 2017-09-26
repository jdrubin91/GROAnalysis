__author__ = 'Jonathan Rubin'

import os
import sys
import matplotlib
matplotlib.use('Agg')
from pybedtools import BedTool
from matplotlib import pyplot as plt
# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})
import numpy as np

def run(rep1,rep2,rep1bam,rep2bam,figuredir,filedir):
    boxplots = [[],[],[]]

    r1 = BedTool(rep1)
    r2 = BedTool(rep2)
    b1 = BedTool(rep1bam)
    b2 = BedTool(rep2bam)
    (r1-r2).multi_bam_coverage(b1,b2).saveas(filedir + "rep1onlyeRNAs_DMSOt45.counts.bed")
    (r2-r1).multi_bam_coverage(b1,b2).saveas(filedir + "rep2onlyeRNAs_DMSOt45.counts.bed")
    (r1+r2).multi_bam_coverage(b1,b2).saveas(filedir + "rep1and2eRNAs_DMSOt45.counts.bed")

    for file1 in [filedir + "rep1onlyeRNAs_DMSOt45.counts.bed",filedir + "rep2onlyeRNAs_DMSOt45.counts.bed",filedir + "rep1and2eRNAs_DMSOt45.counts.bed"]:
        outfile = open(file1+".normalized.bed",'w')
        with open(file1) as F:
            for line in F:
            line = line.strip('\n').split('\t')
            i=0
            for norm in total_mapped:
                line[-4+i] = str(float(line[-4+i])/norm)
                i += 1
            outfile.write('\t'.join(line) + '\n')

    outfile = open(filedir + "Boxplot_ChIPPeak_GROcoverage.counts.normalized.bed",'w')
    with open(filedir + "Boxplot_ChIPPeak_GROcoverage.counts.bed") as F:
        for line in F:
            line = line.strip('\n').split('\t')
            i=0
            for norm in total_mapped:
                line[-4+i] = str(float(line[-4+i])/norm)
                i += 1


if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #Figure directory
    figuredir = parent_dir(homedir) + '/figures/'
    filedir = parent_dir(homedir) + '/files/'

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

    bamfolder1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/'
    bam1 = bamfolder1 + 'J52_trimmed.flip.fastq.bowtie2.sorted.bam'
    bam2 = bamfolder1 + 'J62_trimmed.flip.fastq.bowtie2.sorted.bam'
    bamfolder2 = '/projects/dowellLab/Taatjes/170825_NB501447_0152_fastq/Demux/cat/trimmed/flipped/bowtie2/sortedbam/'
    bam3 = bamfolder2 + 'J5D451_GTCCGC_S3_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    bam4 = bamfolder2 + 'J6C451_GTGAAA_S4_L007and8_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'

    rep1=filelist[10]
    rep2=filelist[11]
    rep1bam = bam1
    rep2bam = bam3

    run(rep1,rep2,rep1bam,rep2bam,figuredir,filedir)