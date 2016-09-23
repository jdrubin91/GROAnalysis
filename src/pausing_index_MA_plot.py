__author__='Jonathan Rubin'

import sys
from pybedtools import BedTool

def run(genedir,bam1,bam2,figdir,filedir):
    TSS = open(filedir+'TSS.bed','w')
    Body = open(filedir+'Body.bed','w')
    Info = open(filedir+'Gene_info.txt','w')
    with open(genedir) as F:
        for line in F:
            line = line.strip().split()
            chrom,start,stop = line[0:3]
            strand = line[5]
            geneName = line[3].split(';')[1]
            TSS.write(chrom+'\t'+str(int(start)-200)+'\t'+str(int(start)+1000)+'\n')
            Body.write(chrom+'\t'+start+'\t'+stop+'\n')
            Info.write(geneName+'\t'+strand+'\n')

    TSS = BedTool(filedir+'TSS.bed')
    Body = BedTool(filedir+'Body.bed')
    Bam1 = BedTool(bam1)
    Bam2 = BedTool(bam2)

    Bam1.intersect(b=TSS,stream=True).count().saveas(filedir+'1_TSS.count.bed')
    Bam1.intersect(b=Body,stream=True).count().saveas(filedir+'1_Body.count.bed')

    Bam2.intersect(b=TSS,stream=True).count().saveas(filedir+'2_TSS.count.bed')
    Bam2.intersect(b=Body,stream=True).count().saveas(filedir+'2_Body.count.bed')


    

if __name__ == "__main__":
    genedir = '/scratch/Users/joru1876/genome_files/refGene.bed'
    bam1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.bam'
    bam2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.bam'
    figdir = '/scratch/Users/joru1876/GROAnalysis/figures/'
    filedir = '/scratch/Users/joru1876/GROAnalysis/files/'
    run(genedir,bam1,bam2,figdir,filedir)