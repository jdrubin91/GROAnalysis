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

    filelist = [filedir + "rep1onlyeRNAs_DMSOt45.counts.bed",filedir + "rep1and2eRNAs_DMSOt45.counts.bed",filedir + "rep2onlyeRNAs_DMSOt45.counts.bed"]
    total_mapped = [[0,0],[0,0],[0,0]]

    j = 0
    for file1 in filelist:
        with open(file1) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                i = 0
                for val in line[-2:]:
                    total_mapped[j][i] += float(val)
                    i += 1
        j += 1

    j = 0
    for file1 in filelist:
        with open(file1) as F:
            for line in F:
                line = line.strip('\n').split('\t')
                val1 = float(line[-2])/total_mapped[j][0]
                val2 = float(line[-1])/total_mapped[j][1]
                boxplot[j].append(val1-val2)
        j += 1

    names = ["Rep1 Only", "Rep1 and 2", "Rep2 Only"]
    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('FPKM difference between replicate subsets DMSO')
    ax.set_ylabel('Difference in FPKM (rep1-rep2)')
    ax.set_xlabel('Subset')
    bp = ax.boxplot(boxplot, positions=np.arange(len(boxplot)),patch_artist=True)
    bp['boxes'][0].set( facecolor = 'red' )
    bp['boxes'][1].set( facecolor = 'yellow' )
    bp['boxes'][2].set( facecolor = 'green' )
    plt.xticks(np.arange(len(boxplot)),names,rotation=45)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.savefig(figuredir + "eRNA_overlap_fpkm_boxplot_DMSOt45.png")
    


if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #Figure directory
    figuredir = parent_dir(homedir) + '/figures/'
    filedir = parent_dir(homedir) + '/files/'

    folder1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit_run2/'
    folder2 = '/projects/dowellLab/Taatjes/170825_NB501447_0152_fastq/Demux/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit/'
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

    rep1=filelist[4]
    rep2=filelist[10]
    rep1bam = bam1
    rep2bam = bam3

    run(rep1,rep2,rep1bam,rep2bam,figuredir,filedir)