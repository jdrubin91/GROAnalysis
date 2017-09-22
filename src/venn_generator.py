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
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def run(rep1,rep2, figuredir,filedir):
    rep1 = BedTool(rep1)
    rep2 = BedTool(rep2)

    AB = len(rep1+rep2)
    aB = len(rep2) - AB
    Ab = len(rep1) - AB
    # (rep1+rep2).saveas(filedir + 'CA_t45_rep1and2_intersect.bed')
    plt.figure()
    v = venn2(subsets=(Ab,aB,AB), set_labels = ('DMSO','CA'))
    plt.title("Bidirectional Overlap DMSO CA Rep1 and Rep2 t=45")
    plt.savefig(figuredir + 'Bidirectional_t45_DMSO_CA_rep1and2_overlap.png',dpi=1200)


def example():
    plt.figure(figsize=(4,4))
    # v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    v = venn2(subsets = (1,2,3), set_labels = ('A','B'))
    # v.get_patch_by_id('100').set_alpha(1.0)
    # v.get_patch_by_id('100').set_color('white')
    # v.get_label_by_id('100').set_text('Unknown')
    # v.get_label_by_id('A').set_text('Set "A"')
    # c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    # c[0].set_lw(1.0)
    # c[0].set_ls('dotted')
    # plt.title("Sample Venn diagram")
    # plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
                 # ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                 # arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()

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

    rep1=filelist[10]
    rep2=filelist[11]

    rep1=filedir + 'DMSO_t45_rep1and2_intersect.bed'
    rep2=filedir + 'CA_t45_rep1and2_intersect.bed'
    print rep1
    print rep2
    run(rep1,rep2,figuredir,filedir)

