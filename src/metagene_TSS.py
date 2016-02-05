__author__ = 'Jonathan Rubin'

import os

file1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
file2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
file3 = '/scratch/Users/joru1876/GROAnalysis/files/Master.bed'
outdir = '/scratch/Users/joru1876/GROAnalysis/files'
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
    
    os.system("bedtools map -a " + outdir + "/TSS_BP_Intervals.bed -b " + file1 + " -c 4 -o sum > " + outdir + "/DMSO_TSS_mapped.bed")
    os.system("bedtools map -a " + outdir + "/TSS_BP_Intervals.bed -b " + file2 + " -c 4 -o sum > " + outdir + "/CA_TSS_mapped.bed")
    
    return


if __name__ == "__main__":
    run(file1,file2,file3)