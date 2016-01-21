__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import intervals

annotationdir = '/scratch/Users/joru1876/genome_files/refGene.bed'
bed1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
bed2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
TSSinterval = (-30,300)

def getLists(file1,file2,file3):
    genes = list()
    TSS = list()
    with open(file1) as F1:
        for line in F1:
            line = line.strip().split()
            chrom,start,stop = line[0:3]
            strand = line[5]
            geneName = line[3].split(';')[1]
            genes.append((int(start),int(stop),chrom,geneName))
            TSS.append((int(start)-TSSinterval[0],int(start)+TSSinterval[1],chrom,geneName,strand))
        
    return genes,TSS
    
def intervalSearch(bed1,bed2,genes,TSS):
    bed1list = list()
    bed2list = list()
    with open(bed1) as F1:
        for line in F1:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                bed1list.append(int(start),int(stop),float(coverage))
                
    with open(bed2) as F2:
        for line in F2:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                bed2list.append(int(start),int(stop),float(coverage))
                
    
    
    return
    
    
if __name__ == "__main__":
    genes,TSS = getLists(annotationdir)
    intervalSearch(bed1,bed2,genes,TSS)