__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from operator import itemgetter

def run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir,figuredir):    
    d = dict()
    with open(DMSOgenes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '1'
            d[gene] = [chrom,start,stop,number,strand,coverage]

    with open(DMSOTSS) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '1'
            d[gene].append(coverage)
            
    with open(DMSOEND) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '1'
            d[gene].append(coverage)
            
    with open(CAgenes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '1'
            d[gene].append(coverage)
            
    with open(CATSS) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '1'
            d[gene].append(coverage)
            
    with open(CAEND) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '1'
            d[gene].append(coverage)
            
    coveragecutoff = 200
    TRlist = list()
    TRgenes = list()
    cutoff1 = 0.25
    ENDlist = list()
    ENDgenes = list()
    cutoff2 = 0.25
    
    outfile = open(filedir + '/Master.bed','w')
    outfile.write('Gene\tChrom\tStart\tStop\tNumber\tStrand\tDMSO gene body\tDMSO TSS\tDMSO END\tCA gene body\tCA TSS\tCA END\n')
    i = 0
    for gene in d:
        i+=1
        outfile.write(gene + '\t')
        for item in d[gene]:
            outfile.write(item + '\t')
        outfile.write('\n')
        DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND = d[gene][5:]
        DMSOgenes = float(DMSOgenes)
        DMSOTSS = float(DMSOTSS)
        DMSOEND = float(DMSOEND)
        CAgenes = float(CAgenes)
        CATSS = float(CATSS)
        CAEND = float(CAEND)
        if CAgenes-CATSS != 0 and DMSOgenes-DMSOTSS != 0 and CAgenes-CAEND != 0 and DMSOgenes-DMSOEND != 0 and DMSOgenes > coveragecutoff and CAgenes > coveragecutoff:
            TR = (CATSS/(CAgenes-CATSS))-(DMSOTSS/(DMSOgenes-DMSOTSS))
            if TR > cutoff1:
                TRgenes.append((gene,TR))
            TRlist.append(TR)
            ER = (CAEND/(CAgenes-CAEND))-(DMSOEND/(DMSOgenes-DMSOEND))
            if ER > cutoff2:
                ENDgenes.append((gene,ER))
            ENDlist.append(ER)
    print "i: ",i
    F1 = plt.figure()
    plt.hist(TRlist,50)
    plt.title("Travelers Ratio")
    plt.savefig(figuredir + '/TravelersRatio.png')
    F2 = plt.figure()
    plt.hist(ENDlist,50)
    plt.title("End Ratio")
    plt.savefig(figuredir + '/EndRatio.png')
    outfile2 = open(filedir + '/GeneList.txt','w')
    for item in sorted(TRgenes, key=itemgetter(1),reverse=True):
        outfile2.write(item[0] + '\t' + str(item[1]) + '\n')
    for item in sorted(ENDgenes, key=itemgetter(1),reverse=True):
        outfile2.write(item[0] + '\t' + str(item[1]) + '\n')
    
    
    
    