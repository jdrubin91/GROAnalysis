__author__ = "Jonathan Rubin"

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import intervals

annotationdir = '/scratch/Users/joru1876/genome_files/refGene.bed'
bed1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
bed2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'
TSSinterval = (-30,300)
threeprimeint = (-100,100)

def getLists(file1):
    TSS = dict()
    TSSgene = dict()
    END = dict()
    ENDgene = dict()
    with open(file1) as F1:
        for line in F1:
            line = line.strip().split()
            chrom,start,stop = line[0:3]
            strand = line[5]
            geneName = line[3].split(';')[1]
            TSS[geneName] = [int(start)+TSSinterval[0],int(start)+TSSinterval[1],chrom,strand]
            TSSgene[geneName] = [int(start)+TSSinterval[1],int(stop),chrom,strand]
            END[geneName] = [int(stop)+threeprimeint[0],int(stop)+threeprimeint[1],chrom,strand]
            ENDgene[geneName] = [int(start),int(stop)+threeprimeint[0],chrom,strand]
        
    return TSS,TSSgene,END,ENDgene
    
def intervalSearch(bed1,bed2,TSS,TSSgene,END,ENDgene):
    bed1list = dict()
    bed2list = dict()
    X1 = list()
    X2 = list()
    Y1 = list()
    Y2 = list()
    with open(bed1) as F1:
        for line in F1:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                if not chrom in bed1list:
                    bed1list[chrom] = list()
                bed1list[chrom].append((int(start),int(stop),float(coverage)))
                
    with open(bed2) as F2:
        for line in F2:
            if '#' not in line[0]:
                chrom, start, stop, coverage = line.strip().split()
                if not chrom in bed2list:
                    bed2list[chrom] = list()
                bed2list[chrom].append((int(start),int(stop),float(coverage)))
                
                          
    print "Performing First Interval Searches..."          
    for gene in TSS:
        cov = 0
        chrom,strand = TSS[gene][2:4]
        ST = intervals.comparison((TSS[gene],bed1list[chrom]))
        OVERLAPS_TSS = ST.find_overlaps(0,1)
        for O in OVERLAPS_TSS:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        TSS[gene].append(cov)
        
    for gene in TSSgene:
        cov = 0
        chrom,strand = TSSgene[gene][2:4]
        ST = intervals.comparison((TSSgene[gene],bed1list[chrom]))
        OVERLAPS_TSSgene = ST.find_overlaps(0,1)
        for O in OVERLAPS_TSSgene:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        TSSgene[gene].append(cov)
        
    for gene in END:
        cov = 0
        chrom,strand = END[gene][2:4]
        ST = intervals.comparison((END[gene],bed1list[chrom]))
        OVERLAPS_END = ST.find_overlaps(0,1)
        for O in OVERLAPS_END:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        END[gene].append(cov)
        
    for gene in ENDgene:
        cov = 0
        chrom,strand = ENDgene[gene][2:4]
        ST = intervals.comparison((ENDgene[gene],bed1list[chrom]))
        OVERLAPS_ENDgene = ST.find_overlaps(0,1)
        for O in OVERLAPS_ENDgene:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        ENDgene[gene].append(cov)
    
    for gene in TSS:
        if gene in TSSgene:
            TSScov = TSS[gene][4]
            TSSgenecov = TSSgene[gene][4]
            X1.append(TSScov/TSSgenecov)
    for gene in END:
        if gene in ENDgene:
            ENDcov = END[gene][4]
            ENDgenecov = ENDgene[gene][4]
            Y1.append(ENDcov/ENDgenecov)
    
    
    print "Performing Second Interval Searches..."          
    for gene in TSS:
        cov = 0
        chrom,strand = TSS[gene][2:4]
        ST = intervals.comparison((TSS[gene],bed2list[chrom]))
        OVERLAPS_TSS = ST.find_overlaps(0,1)
        for O in OVERLAPS_TSS:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        TSS[gene].append(cov)
        
    for gene in TSSgene:
        cov = 0
        chrom,strand = TSSgene[gene][2:4]
        ST = intervals.comparison((TSSgene[gene],bed2list[chrom]))
        OVERLAPS_TSSgene = ST.find_overlaps(0,1)
        for O in OVERLAPS_TSSgene:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        TSSgene[gene].append(cov)
        
    for gene in END:
        cov = 0
        chrom,strand = END[gene][2:4]
        ST = intervals.comparison((END[gene],bed2list[chrom]))
        OVERLAPS_END = ST.find_overlaps(0,1)
        for O in OVERLAPS_END:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        END[gene].append(cov)
        
    for gene in ENDgene:
        cov = 0
        chrom,strand = ENDgene[gene][2:4]
        ST = intervals.comparison((ENDgene[gene],bed2list[chrom]))
        OVERLAPS_ENDgene = ST.find_overlaps(0,1)
        for O in OVERLAPS_ENDgene:
            for interval_original in O.overlaps:
                if not chrom in interval_original.INFO[0]:
                    if '+' in strand:
                        if interval_original.INFO[0] > 0:
                            cov += interval_original.INFO[0]
                    else:
                        if interval_original.INFO[0] < 0:
                            cov += -interval_original.INFO[0]
        ENDgene[gene].append(cov)
        
    for gene in TSS:
        if gene in TSSgene:
            TSScov = TSS[gene][5]
            TSSgenecov = TSSgene[gene][5]
            X2.append(TSScov/TSSgenecov)
    for gene in END:
        if gene in ENDgene:
            ENDcov = END[gene][5]
            ENDgenecov = ENDgene[gene][5]
            Y2.append(ENDcov/ENDgenecov)
    
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