__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
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
    pausing_indexes = list()
    with open(countsfile) as F:
        oldgene = 'none'
        for line in F:
            line = line.strip('\n').split('\t')
            gene = line[3]
            if gene == oldgene:
                denominator = float(line[-1])
                try:
                    pausing_indexes.append((gene,numerator/denominator))
                except:
                    pass
            else:
                numerator = float(line[-1])

    return pausing_indexes

def plot_vs(pausing_indexes1,pausing_indexes2,figuredir):
    x = [l[1] for l in pausing_indexes1]
    y = [l[1] for l in pausing_indexes2]
    print x,y
    F = plt.figure()
    ax = F.add_subplot(111)
    plt.scatter(x,y)
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







