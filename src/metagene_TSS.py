__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np

file1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
file2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
file3 = '/scratch/Users/joru1876/GROAnalysis/files/Master.bed'
outdir = '/scratch/Users/joru1876/GROAnalysis/files'
figout = '/scratch/Users/joru1876/GROAnalysis/figures'
coveragecut = 2000
window = 1000

def run(file1,file2,file3):
    genelist = list()
    with open(file3) as F:
        F.readline()
        for line in F:
            gene,chrom,start,stop,number,strand,coverage = line.strip().split()[0:7]
            if float(coverage) > coveragecut:
                genelist.append(gene)
    
    outfile1 = open(outdir + '/TSS_BP_Intervals.bed','w')
    for gene in genelist:
        start = int(gene.split(';')[2].split('-')[0].split(':')[1])
        chrom = gene.split(';')[2].split('-')[0].split(':')[0]
        for i in range(start-window,start+window):
            outfile1.write(chrom + '\t' + str(i) + '\t' + str(i+1) + '\t' + gene + '\n')
    outfile1.close()
    os.system("sort " + outdir + "/TSS_BP_Intervals.bed -k1,1 -k2,2n > " + outdir + "/TSS_BP_Intervals.sorted.bed")
    os.system("bedtools map -a " + outdir + "/TSS_BP_Intervals.sorted.bed -b " + file1 + " -c 4 -o collapse > " + outdir + "/DMSO_TSS_mapped.bed")
    os.system("bedtools map -a " + outdir + "/TSS_BP_Intervals.sorted.bed -b " + file2 + " -c 4 -o collapse > " + outdir + "/CA_TSS_mapped.bed")
    
    DMSOdict = dict()
    with open(outdir + "/DMSO_TSS_mapped.bed") as F:
        for line in F:
            chrom, start, stop, gene, cov = line.strip().split()
            strand = gene[-1]
            coveragelist = cov.split(',')
            coverage = 0.0
            for item in coveragelist:
                if item is '.':
                    coverage += 0.0
                else:
                    if strand is '-':
                        if float(item) < 0:
                            coverage += -float(item)
                    else:
                        if float(item) > 0:
                            coverage += float(item)
            if gene not in DMSOdict:
                DMSOdict[gene] = np.zeros(window*2)
            TSS = int(gene.split(';')[2].split('-')[0].split(':')[1])
            index = int(start) + window - TSS
            DMSOdict[gene][index] = coverage
            
    CAdict = dict()
    with open(outdir + "/CA_TSS_mapped.bed") as F:
        for line in F:
            chrom, start, stop, gene, cov = line.strip().split()
            strand = gene[-1]
            coveragelist = cov.split(',')
            coverage = 0.0
            for item in coveragelist:
                if item is '.':
                    coverage += 0.0
                else:
                    if strand is '-':
                        if float(item) < 0:
                            coverage += -float(item)
                    else:
                        if float(item) > 0:
                            coverage += float(item)
            if gene not in CAdict:
                CAdict[gene] = np.zeros(window*2)
            TSS = int(gene.split(';')[2].split('-')[0].split(':')[1])
            index = int(start) + window - TSS
            CAdict[gene][index] = coverage
    
    DMSOarray = np.zeros(window*2)
    CAarray = np.zeros(window*2)
    
    for gene in DMSOdict:
        for i in range(len(DMSOdict[gene])):
            DMSOarray[i] += DMSOdict[gene][i]
            
    for gene in CAdict:
        for i in range(len(CAdict[gene])):
            CAarray[i] += CAdict[gene][i]
            
    F = plt.figure()
    plt.plot(DMSOarray)
    plt.plot(CAarray)
    plt.legend(['DMSO', 'CA'], loc='upper left')
    plt.savefig(figout + '/metagene_TSS.png')
    
    
    
    
    return


if __name__ == "__main__":
    run(file1,file2,file3)