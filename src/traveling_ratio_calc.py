__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import intervals

annotationdir = '/scratch/Users/joru1876/genome_files/refGene.bed'
bed1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
bed2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
TSSinterval = (-30,300)
threeprimeint = (-100,100)

def getLists(file1,file2,file3):
    genes = list()
    #TSS = list()
    #threeprime = list()
    with open(file1) as F1:
        for line in F1:
            line = line.strip().split()
            chrom,start,stop = line[0:3]
            strand = line[5]
            geneName = line[3].split(';')[1]
            genes.append((int(start)+TSSinterval[0],int(stop)+threeprimeint[1],chrom,geneName,strand))
            #TSS.append((int(start)+TSSinterval[0],int(start)+TSSinterval[1],chrom,geneName))
            #threeprime.append((int(stop)+threeprimeint[0],int(stop)+threeprimeint[1],chrom,geneName))
        
    return genes
    
def intervalSearch(bed1,bed2,genes):
    bed1list = list()
    bed2list = list()
    x = list()
    y = list()
    with open(bed1) as F1:
        for line in F1:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                bed1list.append(int(start),int(stop),chrom,float(coverage))
                
    with open(bed2) as F2:
        for line in F2:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                bed2list.append(int(start),int(stop),chrom,float(coverage))
    
    print "Comparing intervals..."
    ST1 = intervals.comparison((genes,bed1list))
    OVERLAPS_g_1 = ST1.find_overlaps(0,1)
    print "First comparison done, overlap instances: " + str(len(OVERLAPS_g_1))
    ST2 = intervals.comparison((genes,bed2list))
    OVERLAPS_g_2 = ST2.find_overlaps(0,1)
    print "Second comparison done, overlap instances: " + str(len(OVERLAPS_g_2))
    
    #for O in OVERLAPS_g_1:
    #    for interval_original in O.overlaps:
            
        
    
    
    
    return "done"
    
    
def run(bed1,bed2,genes):
    bed1dict = dict()
    bed2dict = dict()
    TRX = list()
    TRY = list()
    ENDX = list()
    ENDY = list()
    
    print "Creating file dicts..."
    with open(bed1) as F1:
        for line in F1:
            chrom, start, stop, coverage = line.strip().split()
            if not chrom in bed1dict:
                bed1dict[chrom] = list()
            bed1dict[chrom].append((int(start),int(stop),float(coverage)))
    with open(bed2) as F2:
        for line in F2:
            chrom, start, stop, coverage = line.strip().split()
            if not chrom in bed2dict:
                bed2dict[chrom] = list()
            bed2dict[chrom].append((int(start),int(stop),float(coverage)))
    
    print "Finished creating file dicts"
    
    for gene in genes:
        start1,stop1,chrom,name,strand = gene
        intersect1 = list()
        intersect2 = list()
        for item in bed1dict[chrom]:
            start2, stop2, coverage = item
            if start2 <= start1 and stop2 <= stop1:
                if '+' in strand:
                    if coverage > 0:
                        intersect1.append(start2,stop2,coverage)
                else:
                    if coverage < 0:
                        intersect1.append(start2,stop2,coverage)
                
        for item in bed2dict[chrom]:
            start2, stop2, coverage = item
            if start2 <= start1 and stop2 <= stop1:
                if '+' in strand:
                    if coverage > 0:
                        intersect2.append(start2,stop2,coverage)
                else:
                    if coverage < 0:
                        intersect2.append(start2,stop2,coverage)
    
    
if __name__ == "__main__":
    genes = getLists(annotationdir)
    intervalSearch(bed1,bed2,genes)