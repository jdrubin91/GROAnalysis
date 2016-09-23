__author__='Jonathan Rubin'

import sys
import multiprocessing
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pybedtools import BedTool

def intersect(bam,bed,filename):
    print filename
    return BedTool(bed).sort().map(b=bam,c=4,o="sum").saveas(filename)

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
    features = ((bam1,TSS,filedir+'1_TSS.count.bed'),(bam1,END,filedir+'1_END.count.bed'),(bam1,Body,filedir+'1_Body.count.bed'),(bam2,TSS,filedir+'2_TSS.count.bed'),(bam2,END,filedir+'2_END.count.bed'),(bam2,Body,filedir+'2_Body.count.bed'))


    for item in features:
        intersect(item[0],item[1],item[2])
    # feature1=(bam1,bam1,bam1,bam2,bam2,bam2)
    # feature2=(TSS,END,Body,TSS,END,Body)
    # feature3=(filedir+'1_TSS.count.bed',filedir+'1_END.count.bed',filedir+'1_Body.count.bed',filedir+'2_TSS.count.bed',filedir+'2_END.count.bed',filedir+'2_Body.count.bed')

    # results = pool.map(intersect, feature1, feature2, feature3)

    # Bam1.intersect(b=TSS,stream=True).count().saveas(filedir+'1_TSS.count.bed')
    # Bam1.intersect(b=Body,stream=True).count().saveas(filedir+'1_Body.count.bed')

    # Bam2.intersect(b=TSS,stream=True).count().saveas(filedir+'2_TSS.count.bed')
    # Bam2.intersect(b=Body,stream=True).count().saveas(filedir+'2_Body.count.bed')


def plot(TSS1,TSS2,END1,END2,Body1,Body2,genes,figdir):
    X = list()
    Y = list()
    with open(TSS1) as F1, open(TSS2) as F2, open(Body1) as F3, open(Body2) as F4, open(genes) as F5, open(END1) as F6, open(END2) as F7:
        for line1 in F1:
            TSS1=0.0 if line1.strip().split()[-1] == '.' else float(line1.strip().split()[-1])
            TSS2=0.0 if F2.readline().strip().split()[-1] == '.' else float(F2.readline().strip().split()[-1])
            Body1=0.0 if F3.readline().strip().split()[-1] == '.' else float(F3.readline().strip().split()[-1])
            Body2=0.0 if F4.readline().strip().split()[-1] == '.' else float(F4.readline().strip().split()[-1])
            END1=0.0 if F6.readline().strip().split()[-1] == '.' else float(F6.readline().strip().split()[-1])
            END2=0.0 if F7.readline().strip().split()[-1] == '.' else float(F7.readline().strip().split()[-1])
            gene,strand=F5.readline().strip().split()
            Y.append((abs(Body1)+abs(Body2)/2))
            if strand == '+':
                try:
                    print abs(TSS1/(Body1-TSS1)),abs(TSS2/(Body2-TSS2))
                    X.append(math.log(abs(TSS1/(Body1-TSS1)))-math.log(abs(TSS2/(Body2-TSS2))))
                except ZeroDivisionError,ValueError:
                    X.append(0.0)
            else:
                try:
                    X.append(math.log(abs(END1/(Body1-END1))-math.log(abs(END2/(Body2-END2)))))
                except ZeroDivisionError,ValueError:
                    X.append(0.0)

    F = plt.plot(X,Y)
    F.savefig(figdir+'PI_MA_plot.png')







if __name__ == "__main__":
    genedir = '/scratch/Users/joru1876/genome_files/refGene.bed'
    bam1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
    bam2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
    figdir = '/scratch/Users/joru1876/GROAnalysis/figures/'
    filedir = '/scratch/Users/joru1876/GROAnalysis/files/'
    # run(genedir,bam1,bam2,figdir,filedir)

    TSS1=filedir+'1_TSS.count.bed'
    TSS2=filedir+'2_TSS.count.bed'
    END1=filedir+'1_END.count.bed'
    END2=filedir+'2_END.count.bed'
    Body1=filedir+'1_Body.count.bed'
    Body2=filedir+'2_Body.count.bed'
    genes=filedir+'Gene_info.txt'
    plot(TSS1,TSS2,END1,END2,Body1,Body2,genes,figdir)