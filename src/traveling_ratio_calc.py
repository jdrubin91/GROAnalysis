__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import intervals
import time

annotationdir = '/scratch/Users/joru1876/genome_files/refGene.bed'
bed1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
bed2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'
TSSinterval = (-30,300)
threeprimeint = (-100,100)

def getLists(file1):
    TSS = list()
    TSSgene = list()
    END = list()
    ENDgene = list()
    with open(file1) as F1:
        for line in F1:
            line = line.strip().split()
            chrom,start,stop = line[0:3]
            strand = line[5]
            geneName = line[3].split(';')[1]
            TSS.append((int(start)+TSSinterval[0],int(start)+TSSinterval[1],chrom,geneName,strand))
            TSSgene.append((int(start)+TSSinterval[1],int(stop),chrom,geneName,strand))
            END.append((int(stop)+threeprimeint[0],int(stop)+threeprimeint[1],chrom,geneName,strand))
            ENDgene.append((int(start),int(stop)+threeprimeint[0],chrom,geneName,strand))
        
    return TSS,TSSgene,END,ENDgene
    
def intervalSearch(bed1,bed2,TSS,TSSgene,END,ENDgene):
    bed1list = list()
    bed2list = list()
    X1 = list()
    X2 = list()
    Y1 = list()
    Y2 = list()
    with open(bed1) as F1:
        for line in F1:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                bed1list.append((int(start),int(stop),chrom,float(coverage)))
                
    with open(bed2) as F2:
        for line in F2:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                bed2list.append((int(start),int(stop),chrom,float(coverage)))
                
    #print bed1list[0:10]
    #print bed2list[0:10]
    #print TSS[0:10]
    #print TSSgene[0:10]
    #print END[0:10]
    #print ENDgene[0:10]
###############################################################################
    print "Performing First Interval Search..."
    start = time.time()
    ST1 = intervals.comparison((TSS,bed1list))
    OVERLAPS_TSS = ST1.find_overlaps(0,1)
    TSScov = dict()
    print "Finished First Interval Search in: ", time.time() - start
    start = time.time()
    for O in OVERLAPS_TSS:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        TSScov[gene] = cov
    
    print "Finished first interval analysis in: ", time.time() - start
    ST2 = intervals.comparison((TSSgene,bed1list))
    OVERLAPS_TSSgene = ST2.find_overlaps(0,1)
    TSSgenecov = dict()
    for O in OVERLAPS_TSSgene:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        TSSgenecov[gene] = cov
    
    ST3 = intervals.comparison((END,bed1list))
    OVERLAPS_END = ST3.find_overlaps(0,1)
    ENDcov = dict()
    for O in OVERLAPS_END:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        ENDcov[gene] = cov
        
    ST4 = intervals.comparison((ENDgene,bed1list))
    OVERLAPS_ENDgene = ST4.find_overlaps(0,1)
    ENDgenecov = dict()
    for O in OVERLAPS_ENDgene:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        ENDgenecov[gene] = cov
###############################################################################
    for gene in TSScov:
        if gene in TSSgenecov:
            X1.append(TSScov[gene]/TSSgenecov[gene])
        if gene in ENDcov and gene in ENDgenecov:
            Y1. append(ENDcov[gene]/ENDgenecov[gene])
###############################################################################
    ST1 = intervals.comparison((TSS,bed2list))
    OVERLAPS_TSS = ST1.find_overlaps(0,1)
    TSScov = dict()
    print "Finished First Interval Search"
    for O in OVERLAPS_TSS:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        TSScov[gene] = cov
    

    ST2 = intervals.comparison((TSSgene,bed2list))
    OVERLAPS_TSSgene = ST2.find_overlaps(0,1)
    TSSgenecov = dict()
    for O in OVERLAPS_TSSgene:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        TSSgenecov[gene] = cov
    
    ST3 = intervals.comparison((END,bed2list))
    OVERLAPS_END = ST3.find_overlaps(0,1)
    ENDcov = dict()
    for O in OVERLAPS_END:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        ENDcov[gene] = cov
        
    ST4 = intervals.comparison((ENDgene,bed2list))
    OVERLAPS_ENDgene = ST4.find_overlaps(0,1)
    ENDgenecov = dict()
    for O in OVERLAPS_ENDgene:
        cov = 0
        chromlist = list()
        for interval_original in O.overlaps:
            if len(interval_original.INFO) > 2:
                chrom,gene,strand = interval_original.INFO
        for interval_original in O.overlaps:
            if len(interval_original.INFO) < 2:
                if chrom in interval_original.INFO:
                    chromlist.append(interval_original)
        for interval_original in chromlist:
            if '+' in strand:
                if interval_original.INFO[0] > 0:
                    cov += interval_original.INFO[0]
            else:
                if interval_original.INFO[0] < 0:
                    cov += -interval_original.INFO[0]
        ENDgenecov[gene] = cov
###############################################################################
    for gene in TSScov:
        if gene in TSSgenecov:
            X2.append(TSScov[gene]/TSSgenecov[gene])
        if gene in ENDcov and gene in ENDgenecov:
            Y2. append(ENDcov[gene]/ENDgenecov[gene])
###############################################################################
    

    bins = 100
    F = plt.figure()
    plt.hist(X1, bins, alpha=0.5, label='DMSO')
    plt.hist(X2, bins, alpha=0.5, label='CA')
    plt.savefig(savedir + 'Traveling Ratio')
    F2 = plt.figure()
    plt.hist(Y1, bins, alpha=0.5, label='DMSO')
    plt.hist(Y2, bins, alpha=0.5, label='CA')
    plt.savefig(savedir + "3' End Ratio")
    
    
    
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
    TSS,TSSgene,END,ENDgene = getLists(annotationdir)
    intervalSearch(bed1,bed2,TSS,TSSgene,END,ENDgene)