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
    DMSOTRgenes = list()
    cutoff1 = 0.25
    ENDlist = list()
    ENDgenes = list()
    DMSOENDgenes = list()
    cutoff2 = 0.25
    cutoff3 = 0.25
    
    outfile = open(filedir + '/Master.bed','w')
    outfile.write('Gene\tChrom\tStart\tStop\tNumber\tStrand\tDMSO gene body\tDMSO TSS\tDMSO END\tCA gene body\tCA TSS\tCA END\n')
    for gene in d:
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
            if TR < -cutoff1:
                DMSOTRgenes.append((gene,TR))
            if not TR > cutoff3 and not TR < -cutoff3:
                TRlist.append(TR)
            ER = (CAEND/(CAgenes-CAEND))-(DMSOEND/(DMSOgenes-DMSOEND))
            if ER > cutoff2:
                ENDgenes.append((gene,ER))
            if ER < -cutoff2:
                DMSOENDgenes.append((gene,ER))
            if not ER > cutoff3 and not ER < -cutoff3:
                ENDlist.append(ER)
    F1 = plt.figure()
    TRlist.sort(reverse=True)
    #plt.hist(TRlist[int(len(TRlist)*.2):int(len(TRlist)*.8)],50)
    plt.hist(TRlist,50)
    plt.title("Travelers Ratio")
    plt.savefig(figuredir + '/TravelersRatio.png')
    F2 = plt.figure()
    ENDlist.sort(reverse=True)
    #plt.hist(ENDlist[int(len(ENDlist)*.2):int(len(ENDlist)*.8)],50)
    plt.hist(ENDlist,50)
    plt.title("End Ratio")
    plt.savefig(figuredir + '/EndRatio.png')
    outfile2 = open(filedir + '/GeneList.txt','w')
    outfile2.write("High CA TR = " + str(len(TRgenes)) + "\nHigh DMSO TR = " + str(len(DMSOTRgenes)) + "\nHigh CA ER = " + str(len(ENDgenes)) + "\nHigh DMSO ER = " + str(len(DMSOENDgenes)) + "\n")
    outfile2.write("High CA TR\n")
    for item in sorted(TRgenes, key=itemgetter(1),reverse=True):
        outfile2.write(item[0] + '\t' + str(item[1]) + '\n')
    outfile2.write("High DMSO TR\n")
    for item in sorted(DMSOTRgenes, key=itemgetter(1),reverse=True):
        outfile2.write(item[0] + '\t' + str(item[1]) + '\n')
    outfile2.write("High CA ER\n")
    for item in sorted(ENDgenes, key=itemgetter(1),reverse=True):
        outfile2.write(item[0] + '\t' + str(item[1]) + '\n')
    outfile2.write("High DMSO ER\n")
    for item in sorted(DMSOENDgenes, key=itemgetter(1),reverse=True):
        outfile2.write(item[0] + '\t' + str(item[1]) + '\n')
    
    
    
    