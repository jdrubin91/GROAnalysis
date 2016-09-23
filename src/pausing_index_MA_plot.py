__author__='Jonathan Rubin'

import sys
import multiprocessing
from pybedtools import BedTool

def intersect(bam,bed,filename):
    print filename
    return pybedtools.BedTool(bam).intersect(b=bed,stream=True).count().saveas(filename)

def run(genedir,bam1,bam2,figdir,filedir):
    TSS = open(filedir+'TSS.bed','w')
    END = open(filedir+'END.bed','w')
    Body = open(filedir+'Body.bed','w')
    Info = open(filedir+'Gene_info.txt','w')
    with open(genedir) as F:
        for line in F:
            line = line.strip().split()
            chrom,start,stop = line[0:3]
            strand = line[5]
            geneName = line[3].split(';')[1]
            TSS.write(chrom+'\t'+str(int(start)-200)+'\t'+str(int(start)+1000)+'\n')
            END.write(chrom+'\t'+str(int(stop)-1000)+'\t'+str(int(stop)+200)+'\n')
            Body.write(chrom+'\t'+start+'\t'+stop+'\n')
            Info.write(geneName+'\t'+strand+'\n')

    TSS.close()
    END.close()
    Body.close()
    Info.close()

    pool = multiprocessing.Pool()
    TSS = filedir+'TSS.bed'
    END = filedir+'END.bed'
    Body = filedir+'Body.bed'
    features = [(bam1,TSS,filedir+'1_TSS.count.bed'),(bam1,END,filedir+'1_END.count.bed'),(bam1,Body,filedir+'1_Body.count.bed'),(bam2,TSS,filedir+'2_TSS.count.bed'),(bam2,END,filedir+'2_END.count.bed'),(bam2,Body,filedir+'2_Body.count.bed')]

    results = pool.map(intersect, features)

    # Bam1.intersect(b=TSS,stream=True).count().saveas(filedir+'1_TSS.count.bed')
    # Bam1.intersect(b=Body,stream=True).count().saveas(filedir+'1_Body.count.bed')

    # Bam2.intersect(b=TSS,stream=True).count().saveas(filedir+'2_TSS.count.bed')
    # Bam2.intersect(b=Body,stream=True).count().saveas(filedir+'2_Body.count.bed')


def plot(TSS1,TSS2,Body1,Body2,genes):
    X = list()
    Y = list()
    with open(TSS1),open(TSS2),open(Body1),open(Body2),open(genes) as F1,F2,F3,F4,F5:
        for line1 in F1:
            TSS1=float(line1.strip().split()[-1])
            TSS2=float(F2.readline().strip().split()[-1])
            Body1=float(F3.readline().strip().split()[-1])
            Body2=float(F4.readline().strip().split()[-1])
            gene,strand=F5.readline().strip().split()





if __name__ == "__main__":
    genedir = '/scratch/Users/joru1876/genome_files/refGene.bed'
    bam1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.bam'
    bam2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.bam'
    figdir = '/scratch/Users/joru1876/GROAnalysis/figures/'
    filedir = '/scratch/Users/joru1876/GROAnalysis/files/'
    run(genedir,bam1,bam2,figdir,filedir)

    TSS1='1_TSS.count.bed'
    TSS2='2_TSS.count.bed'
    Body1='1_Body.count.bed'
    Body2='2_Body.count.bed'
    genes='Gene_info.txt'
    # plot(TSS1,TSS2,Body1,Body2,genes)