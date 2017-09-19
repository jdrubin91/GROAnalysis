__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
from sklearn import linear_model
import math
import numpy as np
from scipy.stats import gaussian_kde
import numpy.polynomial.polynomial as poly
import pybedtools
import sys

def get_gc(sequence):
    gc = 0.0
    for letter in sequence:
        if letter == 'G' or letter == 'C' or letter == 'g' or letter == 'c':
            gc += 1.0

    return gc/len(sequence)

def show_value(s):
     """
     Convert unicode to str under Python 2;
     all other values pass through unchanged
     """
     if sys.version_info.major == 2:
         if isinstance(s, unicode):
             return str(s)
     return s

def get_factor_norm(countsfile,name_order,header):
    factornorm = [0]*len(name_order)
    with open(countsfile) as F:
        if header:
            F.readline()
        for line in F:
            line = line.strip('\n').split('\t')[6:]
            for i in range(len(line)):
                factornorm[i] += float(line[i])


    factornorm = [factornorm[0]/factornorm[i] for i in range(len(factornorm))]

    return factornorm

def MA_plot(x,y,c,condition1,condition2,figuredir):
    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(x,y,c=c,edgecolor="",s=14,alpha=0.5) 
    ax.set_title('Gene Expression MA Plot')
    ax.set_ylabel('Log2 Fold-Change (' + condition1 + '/' + condition2 + ')')
    ax.set_xlabel('Log10 Average Expression')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.colorbar()

    # plt.show()
    plt.savefig(figuredir + 'GeneExpression_GCcontent_RNASeq_normalized_linear.png',dpi=1200)
    # plt.savefig(figuredir + 'ExonExpression_GCcontent_RNASeq_normalized_poly2.svg')

def gc_scatter(exp,gc_content,condition1,condition2,figuredir):

    # slope,intercept = np.polyfit(gc_content, exp, 1)
    # print 'slope:',slope,' intercept:',intercept

    x_new = np.linspace(0, 1, num=len(gc_content)*10)
    coefs = poly.polyfit(gc_content, exp, 1)
    print coefs
    ffit = poly.polyval(x_new, coefs)
    ffit = poly.Polynomial(coefs)


    F = plt.figure() 
    ax = F.add_subplot(111)
    # print "starting vstack..."
    # xy = np.vstack([gc_content,exp])
    # print "done\nstarting gaussian..."
    # z = gaussian_kde(xy)(xy)
    # print "done\nstarting scatter..."
    plt.scatter(gc_content,exp,edgecolor="",s=14,alpha=0.5) 
    # print "done"
    # ax.plot([0,1],[intercept,slope*1],color = 'r')
    # plt.plot(x_new, ffit(x_new),color='r')
    # plt.plot(x_new, ffit)
    ax.set_title('Expression Fold-Change vs. GC-content')
    ax.set_ylabel('Log2 Fold-Change (' + condition1 + '/' + condition2 + ')')
    ax.set_xlabel('GC-content')
    ax.set_xlim([0, 1])
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.savefig(figuredir + 'Gene_Scatter_FoldChange_GCContent_normalized_linear.png', dpi=1200)
    # plt.savefig(figuredir + 'Scatter_FoldChange_GCContent_poly2.svg')

def run(exonfasta,countsfile,name_order,condition1,condition2,figuredir,header=False):
    fastadict = dict()
    with open(exonfasta) as F:
        for line in F:
            if '>' in line[0]:
                key = line[1:].strip('\n')
            else:
                fastadict[key] = line

    x = list()
    y = list()
    c = list()
    exp = list()
    factornorm = get_factor_norm(countsfile,name_order,header=header)

    with open(countsfile) as F:
        if header:
            F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            index1 = name_order.index(condition1)
            index2 = name_order.index(condition2)
            basemean1 = float(line[6+index1])*factornorm[index1]
            basemean2 = float(line[6+index2])*factornorm[index2]
            avg = (basemean1+basemean2)/2.0
            key = line[0] + ':' + line[1] + '-' + line[2]
            gc = get_gc(fastadict[key])
            # if avg > 10:
            try:
                #Exon Normalizations
                # y.append(math.log(basemean1/(basemean2*(2**((-5.65168341607*gc) + 3.09298861697))),2))
                # y.append(math.log(basemean1/(basemean2*(2**((-4.66315953327*gc) + 2.55580519539))),2))
                # y.append(math.log(basemean1/(basemean2*(2**((-2.55184172504*gc) + 1.36042946199))),2))
                # y.append(math.log(basemean1/(basemean2*(2**((-5.00105778*gc) + 2.71452701))),2))
                # y.append(math.log(basemean1/(basemean2*(2**((10.19644627*(gc**2)) + (-15.5988791*gc) + 5.3616822))),2))
                # y.append(math.log(basemean1/(basemean2*(2**((5.2754932*(gc**3)) + (1.7904469*(gc**2)) + (-11.26449439*gc) + 4.64027801))),2))

                #Gene Normalizations
                y.append(math.log(basemean1/(basemean2*(2**((-5.30326135*gc) + 2.5648209))),2))
                # y.append(math.log(basemean1/basemean2))
                x.append(math.log((basemean1+basemean2)/2.0,10))
                c.append(gc)
            except:
                # pass
                y.append(0)
                x.append(0)
                c.append(gc)
                


    MA_plot(x,y,c,condition1,condition2,figuredir)
    # gc_scatter(y,c,condition1,condition2,figuredir)

#This function does not work right now. It was an early implementation of this method using pybedtools to correct a counts file. It was abandoned because
#pybedtools simply takes too long
def correct_counts_file_pybedtools(name_order, condition1, exonfasta, countsfile, genecountsfile):
    index1 = name_order.index(condition1)
    fastadict = dict()
    with open(exonfasta) as F:
        for line in F:
            if '>' in line[0]:
                key = line[1:].strip('\n')
            else:
                fastadict[key] = line

    a = pybedtools.BedTool(countsfile).cut([0,1,2])
    outfile = open(genecountsfile + '.gc_normalized_NW1.bed','w')
    with open(genecountsfile) as F:
        outfile.write(F.readline())
        counter = 0
        for line in F:
            line=line.strip('\n').split('\t')
            basemean1 = float(line[6+index1])
            genecoord = line[3].split(';')[-1].split(':')[0] + '\t' + line[3].split(';')[-1].split(':')[1].split('_')[0].split('-')[0] + '\t' + line[3].split(';')[-1].split(':')[1].split('_')[0].split('-')[1]
            # print genecoord
            b = pybedtools.BedTool(genecoord,from_string=True)
            c = b.intersect(b=a,stream=True,wb=True)
            print c
            seq = ''
            for interval in c:
                key = show_value(interval.chrom) + ':' + str(show_value(interval.start)) + '-' + str(show_value(interval.stop))
                seq += fastadict[key]
            if len(seq) != 0:
                gc = get_gc(seq)
                basemean1_normalized = basemean1*(2**((10.19644627*(gc**2)) + (-15.5988791*gc) + 5.3616822))
                outfile.write('\t'.join(line[:6+index1]) + '\t' + str(basemean1_normalized) + '\n')

            print counter
            counter += 1


def correct_counts_file(name_order, condition1, exonfasta,gene_exon_countsfile):
    index1 = name_order.index(condition1)
    fastadict = dict()
    with open(exonfasta) as F:
        for line in F:
            if '>' in line[0]:
                key = line[1:].strip('\n')
            else:
                fastadict[key] = line

    outfile = open(gene_exon_countsfile.split('.intersect.')[0]+'.gc_normalized_NW2.bed','w')
    with open(gene_exon_countsfile) as F:
        oldgene = ''
        for line in F:
            line = line.strip('\n').split('\t')
            currentgene = line[3]
            if oldgene == '':
                exons = list()
                exons.append(line[14] + ':' + line[15] + '-' + line[16])
            elif oldgene != currentgene:
                exons = set(list((exons)))
                seq = ''
                for exon in exons:
                    seq += fastadict[exon]
                gc = get_gc(seq)
                basemean1 = float(line[6+index1])
                basemean1_normalized = basemean1*(2**((10.19644627*(gc**2)) + (-15.5988791*gc) + 5.3616822))
                outfile.write('\t'.join(line[:6+index1]) + '\t' + str(int(basemean1_normalized)) + '\n')
                exons = list()
                exons.append(line[14] + ':' + line[15] + '-' + line[16])
            else:
                exons.append(line[14] + ':' + line[15] + '-' + line[16])

            oldgene = currentgene

def correct_exon_counts_file(name_order, condition1, exonfasta,countsfile,header=False):
    index1 = name_order.index(condition1)
    fastadict = dict()
    with open(exonfasta) as F:
        for line in F:
            if '>' in line[0]:
                key = line[1:].strip('\n')
            else:
                fastadict[key] = line

    outfile = open(countsfile.split('.bed.')[0]+'.count.gc_normalized_NW2.bed','w')
    with open(countsfile) as F:
        if header:
            F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            basemean1 = float(line[6+index1])
            key = line[0] + ':' + line[1] + '-' + line[2]
            gc = get_gc(fastadict[key])
            #For Exons
            # basemean1_normalized = basemean1*(2**((10.19644627*(gc**2)) + (-15.5988791*gc) + 5.3616822))

            #For Genes
            basemean1_normalized = basemean1*(2**((-5.30326135*gc) + 2.5648209))
            outfile.write('\t'.join(line[:6+index1]) + '\t' + str(int(basemean1_normalized)) + '\n')


def get_gene_from_exonlist(gene_exon_countsfile, exonlist):
    genelist = list()
    with open(gene_exon_countsfile) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            exon = line[17]
            if exon in exonlist:
                genelist.append(line[3].split(';')[1])

    genelist = set(list(genelist))
    for item in genelist:
        print item



if __name__ == "__main__":
    #The purpose of this script is to compare the expression levels of two replicates when one replicate is thought to have some GC bias.
    #The general thought is that if some GC bias has been introduced into library prep/sequencing, then genes/exons that have higher GC-content will be artificially elevated/lowered
    #The pipeline for determining this and correcting is as follows:
        #1. Generate MA plot where dots are colored by GC-content
        #2. Generate scatter plot of fold change vs. GC-content
        #3. Perform a linear fit of the data - determine equation of line (for JDR 8/3017: y = -5.65168341607x + 3.09298861697)
        #4. Transform affected replicate (for JDR 8/30/17: 
                                            #log2(NW1/NW2) = -5.65168341607(GC-content) + 3.09298861697)
                                            #want: 0 = -5.65168341607(GC-content) + 3.09298861697
                                            #NW1/NW2 = 2^(-5.65168341607(GC-content) + 3.09298861697)
                                            #NW2 is known so we can calculate theoretical NW1 by:
                                            #NW1 = 2^(-5.65168341607(GC-content) + 3.09298861697)*NW2
                                            #Since we want NW1 = NW2 (replicates) then 2^(-5.65168341607(GC-content) + 3.09298861697) is the normalization factor to use
                                            ###For polynomial 2 fit, the normalization factor becomes 2**((14.71036247*(gc**2)) + (-19.9525838*gc) + 6.37484286))
    exonfasta = '../../hg19_reference_files/hg19_exons.fa'
    genefasta = '../../misc/hg19.genes.fa'
    countsfile = '../../misc/hg19_exons.bed.count.bed'
    genecountsfile = '../files/all_samples.genes.bed.count.bed'
    gene_exon_countsfile = '../files/all_samples.genes.count.intersect.exons.count.bed'
    name_order = ['DP1','DP2','DW1','DW2','NP1','NP2','NW1','NW2']
    condition1 = 'NW1'
    condition2 = 'NW2'
    figuredir = '../figures/'
    exonlist = ['uc001sun.5_exon_1_0_chr12_69202988_f', 'uc021yzb.1_exon_3_0_chr6_36653528_f', 'uc003omm.4_exon_1_0_chr6_36651874_f', 'uc002pgf.4_exon_0_0_chr19_47724079_r', 'uc021rai.2_exon_1_0_chr12_69207334_f', 'uc021rai.2_exon_0_0_chr12_69202988_f', 'uc003omn.3_exon_2_0_chr6_36653528_f', 'uc001sui.4_exon_3_0_chr12_69210592_f', 'uc011dtq.2_exon_2_0_chr6_36653528_f', 'uc021rag.2_exon_0_0_chr12_69202988_f', 'uc021rad.2_exon_1_0_chr12_69202988_f', 'uc021rae.1_exon_3_0_chr12_69210592_f', 'uc010ekz.3_exon_0_0_chr19_47724079_r', 'uc001gzq.3_exon_1_0_chr1_203276232_f', 'uc001gwq.4_exon_1_0_chr1_201437469_r', 'uc009zqx.4_exon_2_0_chr12_69207334_f', 'uc001sun.5_exon_2_0_chr12_69207334_f', 'uc031qie.1_exon_2_0_chr12_69233382_f', 'uc009zrh.4_exon_0_0_chr12_69202988_f', 'uc010eky.3_exon_0_0_chr19_47724079_r', 'uc009zrc.4_exon_0_0_chr12_69202988_f', 'uc021rag.2_exon_2_0_chr12_69210592_f', 'uc010xyl.2_exon_0_0_chr19_47724079_r', 'uc009zrb.1_exon_0_0_chr12_69202988_f', 'uc003omn.3_exon_1_0_chr6_36651874_f', 'uc009zqy.1_exon_1_0_chr12_69202988_f', 'uc009zrb.1_exon_1_0_chr12_69207334_f', 'uc021yzb.1_exon_2_0_chr6_36651874_f', 'uc009vmq.3_exon_0_0_chr1_9208346_r', 'uc001gwq.4_exon_0_0_chr1_201434607_r', 'uc009zrd.4_exon_0_0_chr12_69202988_f', 'uc021rag.2_exon_1_0_chr12_69207334_f', 'uc001bps.3_exon_9_0_chr1_28607227_f', 'uc009zra.4_exon_2_0_chr12_69207334_f', 'uc009zrf.4_exon_0_0_chr12_69202988_f', 'uc021raf.2_exon_0_0_chr12_69202988_f', 'uc009zqy.1_exon_2_0_chr12_69207334_f', 'uc009zqx.4_exon_3_0_chr12_69210592_f', 'uc021raf.2_exon_1_0_chr12_69207334_f', 'uc009zrg.4_exon_0_0_chr12_69202988_f', 'uc009zqy.1_exon_3_0_chr12_69210592_f', 'uc001suo.4_exon_0_0_chr12_69202988_f', 'uc031qie.1_exon_1_0_chr12_69207334_f', 'uc009zrc.4_exon_1_0_chr12_69207334_f', 'uc021phf.1_exon_0_0_chr1_201437531_r', 'uc001sui.4_exon_1_0_chr12_69202988_f', 'uc031prm.1_exon_2_0_chr1_201437473_r', 'uc009zrb.1_exon_3_0_chr12_69210592_f', 'uc009zre.4_exon_1_0_chr12_69207334_f', 'uc003omm.4_exon_2_0_chr6_36653528_f', 'uc011dtq.2_exon_1_0_chr6_36651874_f', 'uc021rah.2_exon_0_0_chr12_69202988_f', 'uc021rae.1_exon_1_0_chr12_69202988_f', 'uc021raj.2_exon_0_0_chr12_69202988_f', 'uc021yzc.1_exon_1_0_chr6_36651874_f', 'uc031qie.1_exon_0_0_chr12_69202988_f', 'uc021yzc.1_exon_2_0_chr6_36653528_f', 'uc021rad.2_exon_2_0_chr12_69207334_f', 'uc031prm.1_exon_0_0_chr1_201434607_r', 'uc009zre.4_exon_0_0_chr12_69202988_f', 'uc021raj.2_exon_1_0_chr12_69207334_f', 'uc009zqx.4_exon_1_0_chr12_69202988_f', 'uc021rae.1_exon_2_0_chr12_69207334_f', 'uc021rah.2_exon_1_0_chr12_69210592_f', 'uc001sui.4_exon_2_0_chr12_69207334_f', 'uc009zra.4_exon_1_0_chr12_69202988_f']
    # run(exonfasta,countsfile,name_order,condition1,condition2,figuredir)
    # correct_counts_file_pybedtools(name_order, condition2, exonfasta, countsfile, genecountsfile)
    # correct_counts_file(name_order, condition2, exonfasta,gene_exon_countsfile)
    # correct_exon_counts_file(name_order, condition2, exonfasta,countsfile)
    # get_gene_from_exonlist(gene_exon_countsfile, exonlist)

    # run(genefasta,genecountsfile,name_order,condition1,condition2,figuredir,header=True)
    correct_exon_counts_file(name_order, condition2, genefasta,genecountsfile,header=True)





