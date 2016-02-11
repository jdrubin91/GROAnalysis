__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np

file1 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
file2 = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'
file3 = '/scratch/Users/joru1876/GROAnalysis/files/Master.bed'
outdir = '/scratch/Users/joru1876/GROAnalysis/files'
figout = '/scratch/Users/joru1876/GROAnalysis/figures'
coveragecut = 2500
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
    outfile1.close()
    os.system("sort " + outdir + "/TSS_BP_Intervals.bed -k1,1 -k2,2n > " + outdir + "/TSS_BP_Intervals.sorted.bed")
    os.system("bedtools map -a " + outdir + "/TSS_BP_Intervals.sorted.bed -b " + file1 + " -c 4 -o collapse > " + outdir + "/DMSO_TSS_mapped.bed")
    os.system("bedtools map -a " + outdir + "/TSS_BP_Intervals.sorted.bed -b " + file2 + " -c 4 -o collapse > " + outdir + "/CA_TSS_mapped.bed")
    
    DMSOdict = dict()
    DMSOantidict = dict()
    with open(outdir + "/DMSO_TSS_mapped.bed") as F:
        for line in F:
            chrom, start, stop, gene, cov = line.strip().split()
            strand = gene[-1]
            coveragelist = cov.split(',')
            coverage = 0.0
            antisense = 0.0
            for item in coveragelist:
                if item is '.':
                    coverage += 0.0
                else:
                    if strand is '-':
                        if float(item) < 0:
                            coverage += -float(item)
                        if float(item) > 0:
                            antisense += -float(item)
                    else:
                        if float(item) > 0:
                            coverage += float(item)
                        if float(item) < 0:
                            antisense += float(item)
            if gene not in DMSOdict:
                DMSOdict[gene] = np.zeros(window*2)
                DMSOantidict[gene] = np.zeros(window*2)
            TSS = int(gene.split(';')[2].split('-')[0].split(':')[1])
            index = int(start) + window - TSS
            DMSOdict[gene][index] = coverage
            DMSOantidict[gene][index] = antisense
            
    CAdict = dict()
    CAantidict = dict()
    with open(outdir + "/CA_TSS_mapped.bed") as F:
        for line in F:
            chrom, start, stop, gene, cov = line.strip().split()
            strand = gene[-1]
            coveragelist = cov.split(',')
            coverage = 0.0
            antisense = 0.0
            for item in coveragelist:
                if item is '.':
                    coverage += 0.0
                else:
                    if strand is '-':
                        if float(item) < 0:
                            coverage += -float(item)
                        if float(item) > 0:
                            antisense += -float(item)
                    else:
                        if float(item) > 0:
                            coverage += float(item)
                        if float(item) < 0:
                            antisense += float(item)
            if gene not in CAdict:
                CAdict[gene] = np.zeros(window*2)
                CAantidict[gene] = np.zeros(window*2)
            TSS = int(gene.split(';')[2].split('-')[0].split(':')[1])
            index = int(start) + window - TSS
            CAdict[gene][index] = coverage
            CAantidict[gene][index] = antisense
    
    DMSOarray = np.zeros(window*2)
    DMSOantiarray = np.zeros(window*2)
    CAarray = np.zeros(window*2)
    CAantiarray = np.zeros(window*2)
    
    
    for gene in DMSOdict:
        maximum = np.amax(DMSOdict[gene])
        maximumanti = np.amin(DMSOantidict[gene])
        for i in range(len(DMSOdict[gene])):
            if maximum != 0:
                DMSOarray[i] += DMSOdict[gene][i]/maximum
            if maximumanti != 0:
                DMSOantiarray[i] += -DMSOantidict[gene][i]/maximumanti
            
    for gene in CAdict:
        maximum = np.amax(CAdict[gene])
        maximumanti = np.amin(CAantidict[gene])
        for i in range(len(CAdict[gene])):
            if maximum != 0:
                CAarray[i] += CAdict[gene][i]/maximum
            if maximumanti != 0:
                CAantiarray[i] += -CAantidict[gene][i]/maximumanti
    
    DMSOmax = np.amax(DMSOarray)
    DMSOantimax = np.amin(DMSOantiarray)
    CAmax = np.amax(CAarray)
    CAantimax = np.amin(CAantiarray)
    
    for i in range(window*2):
        DMSOarray[i] = DMSOarray[i]/DMSOmax
        DMSOantiarray[i] = -DMSOantiarray[i]/DMSOantimax
        CAarray[i] = CAarray[i]/CAmax
        CAantiarray[i] = -CAantiarray[i]/CAantimax
        
    
    CAtime = 0
    DMSOtime = 0
    CAtot = 0
    DMSOtot = 0
    for val in CAarray-DMSOarray:
        if val > 0:
            CAtime += 1
            CAtot += val
        elif val < 0:
            DMSOtime += 1
            DMSOtot += -val
        
    print "Time in higher CA: ", CAtime
    print "Time in higher DMSO: ", DMSOtime
    print "CA integration: ", CAtot
    print "DMSO integration: ", DMSOtot
    
    CAtimea = 0
    DMSOtimea = 0
    CAtota = 0
    DMSOtota = 0
    for val in CAantiarray-DMSOantiarray:
        if val > 0:
            CAtimea += 1
            CAtota += val
        elif val < 0:
            DMSOtimea += 1
            DMSOtota += -val
            
    print "Time in higher CA antisense: ", CAtimea
    print "Time in higher DMSO antisense: ", DMSOtimea
    print "CA integration antisense: ", CAtota
    print "DMSO integration antisense: ", DMSOtota
                
    F = plt.figure()
    x1 = np.arange(-window,window,1)
    plt.plot(x1,DMSOarray,color='b')
    plt.plot(x1,CAarray,color='g')
    plt.plot(x1,DMSOantiarray,color='r')
    plt.plot(x1,CAantiarray,color='y')
    plt.xlabel('TSS')
    plt.title('TSS Metagene Analysis')
    plt.axvline(x=0.,color='k',ls='dashed')
    plt.legend(['DMSO', 'CA', 'DMSO antisense','CA antisense'], fontsize=8, loc='lower right')
    plt.savefig(figout + '/metagene_TSS.png')
    
    F1 = plt.figure()
    x2 = np.arange(-window,window,1)
    ax1 = F1.add_subplot(2,1,1)
    ax1.plot(x2,CAarray-DMSOarray)
    ax1.set_title('CA-DMSO')
    ax1.text(450,0.15,"Time in higher CA: " + str(CAtime) + "\nTime in higher DMSO: " + str(DMSOtime) + "\nCA integration: " + str(CAtot)[:5] + "\nDMSO integration: " + str(DMSOtot)[:5],fontsize=8)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax2 = F1.add_subplot(2,1,2)
    ax2.plot(x2,CAantiarray-DMSOantiarray,color='r')
    ax2.set_title('CA-DMSO antisense')
    ax2.text(250,0.3,"Time in higher CA antisense: " + str(CAtimea) + "\nTime in higher DMSO antisense: " + str(DMSOtimea) + "\nCA integration antisense: " + str(CAtota)[:5] + "\nDMSO integration antisense: " + str(DMSOtota)[:5],fontsize=8)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    #plt.legend(['CA-DMSO','DMSO-CA antisense'], loc='upper left',fontsize=8)
    plt.xlabel('TSS')
    ax1.axvline(x=0.,color='k',ls='dashed')
    ax1.axhline(y=0.,color='k',ls='solid')
    ax2.axvline(x=0.,color='k',ls='dashed')
    ax2.axhline(y=0.,color='k',ls='solid')
    #plt.title("Coverage Ratio")
    plt.savefig(figout + '/metagene_TSS_ratio.png')
    
    
    
    
    return


if __name__ == "__main__":
    run(file1,file2,file3)