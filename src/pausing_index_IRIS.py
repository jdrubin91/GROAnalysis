__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import os

#This script takes 

def split_bed(gene_annotations,split_bed_file,upstream,downstream):
    minlength = (upstream+downstream)*2
    outfile = open(split_bed_file,'w')
    with open(gene_annotations) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop,gene,val,strand = line
            start = int(start)
            stop = int(stop)
            if stop-start > minlength:
                if strand == '+':
                    outfile.write('\t'.join([chrom,str(start-upstream),str(start+downstream),gene,val,strand]) + '\n')
                    outfile.write('\t'.join([chrom,str(start+downstream),str(stop),gene,val,strand]) + '\n')
                else:
                    outfile.write('\t'.join([chrom,str(stop-downstream),str(stop+upstream),gene,val,strand]) + '\n')
                    outfile.write('\t'.join([chrom,str(start),str(stop-downstream),gene,val,strand]) + '\n')


def run(split_bed,bam):
    countsfile = split_bed + ".counts.bed"
    os.system("bedtools multicov -bams " + bam + " -bed " + split_bed + " > " + countsfile)
    pausing_indexes = dict()
    with open(countsfile) as F:
        oldgene = 'none'
        for line in F:
            line = line.strip('\n').split('\t')
            gene = line[3]
            if gene == oldgene:
                denominator = float(line[-1])
                try:
                    pausing_indexes[gene]=numerator/denominator
                except:
                    pass
            else:
                numerator = float(line[-1])
            oldgene = gene

    return pausing_indexes

def plot_vs(pausing_indexes1,pausing_indexes2,figuredir):
    genes = set(pausing_indexes1.keys()) & set(pausing_indexes2.keys())
    x = list()
    y = list()
    for gene in genes:
        if gene.split(';')[1] in ['Apof', 'Pdlim5', 'Plekhg6', 'Lima1', 'Gm12216', 'Nfkbiz', 'Irf1', 'Arid5b', 'Fas', 'Irf9', 'Irf8', 'Clec2d', 'Rhoj', 'Zfp36', '2310001H17Rik', 'Gpx4', 'Helz2', 'Serpina3f', 'Bcl3', 'Nfic', 'Adar', 'Thbs1', 'Stambpl1', 'Car9', 'Olfr56', 'Agfg2', 'Cttnbp2nl', 'Serpina3i', 'Cebpd', 'Plxna3', 'Tubb3', 'Msrb3', 'Serpina3g', 'Stat5a', 'Mgat3', 'Ly6c1', 'Actn1', 'Socs3', 'Ifnar2', 'Anxa3', 'Crem', 'Psmb9', 'Fosl2', 'Bvht', 'Relb', 'Cbx2', 'Stat1', 'Filip1l', 'Tnfrsf1a', 'Arf2', 'Eif2ak2', 'Impa2', 'Eid3', 'Tmem173', 'Ifi47', 'Syn3', 'Tap1', 'Stat3', 'Stat2', 'Ly6a', 'Mb', 'Mir143hg', 'Ap3m1', 'Osmr', 'Tgfb1i1', 'Csrp1', 'Fblim1', 'Il1r1', 'Tagln', 'Fermt2', 'Gm4285', 'Gm16675', 'Parp3', 'Mob3c', 'Mob3a', 'Palm', 'Frmd6', 'Txnrd1', 'Sbno2', 'Litaf', 'Gbp7', 'Pcsk7', 'Mnda', 'Ccrn4l', 'Ttc39b', 'Alpk1', 'Wdr1', 'Actg2', 'Klf6', 'Tldc2', 'Junb', 'Ntn4', 'Agpat4', 'Nlrc5', 'Rnf31', 'Spata13', 'Trerf1', 'Tspan11', 'Azin2', '9330175E14Rik', 'Wdfy1', 'Irgm2', 'Serpine1', 'Irgm1', 'F3', 'Timp3', 'Trim25', 'Trim21', 'Fhl2', 'Palld', 'Fhl1', 'Noc4l', 'Oas1b', 'Dact1', '4930562C15Rik', 'Pls3', 'Myl12b', 'Vcl', 'Apol6', 'Gyltl1b', 'Gm15910', 'Xaf1', 'Gnptab', 'Tgtp1', 'Tgtp2', 'Shisa4', 'Rnf19b', 'Il23a', 'Cdc42ep3', 'Gm12185', 'Cish', 'Btg2', 'A530013C23Rik', '9930111J21Rik2', '9930111J21Rik1', 'Rbm47', 'Prss48', 'Ripk1', 'Taf10', 'Igtp', 'Jak3', 'Bst2']:
            x.append(pausing_indexes1[gene])
            y.append(pausing_indexes2[gene])
    print x,y
    F = plt.figure()
    ax = F.add_subplot(111)
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    ax.scatter(x,y,c=z,edgecolor="",s=14)
    plt.xlabel("30 min IFN + DMSO")
    plt.ylabel("30 min IFN + CA")
    plt.ylim((0,3.5))
    plt.xlim((0,3.5))
    ax.plot([0,50.0],[0,50.0],color='k')
    plt.savefig(figuredir + 'pausing_indexes_IRIS.png')

def plot_boxplots(array_of_pausing_indexes):
    print "not done with this yet"



if __name__ == "__main__":
    upstream = 200
    downstream = 500
    gene_annotations = '/scratch/Users/joru1876/mm10/mm10.refFlat.bed'
    split_bed_file = '/scratch/Users/joru1876/GROAnalysis/files/pausing_index_IRIS_split_bed.bed'
    split_bed(gene_annotations,split_bed_file,upstream,downstream)

    bamdir = '/scratch/Users/joru1876/Taatjes/171026_NB501447_0180_fastq_IRISREP2/Demux/Taatjes-374/trimmed/flipped/bowtie2/sortedbam/'
    bam1 = bamdir + '30_2_S3_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    pausing_indexes1 = run(split_bed_file,bam1)

    bam2 = bamdir + '30_CA_2_S4_R1_001_trimmed.flip.fastq.bowtie2.sorted.bam'
    pausing_indexes2 = run(split_bed_file,bam2)

    figuredir = '/scratch/Users/joru1876/GROAnalysis/figures/'
    plot_vs(pausing_indexes1,pausing_indexes2,figuredir)







