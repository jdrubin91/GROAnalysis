__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
from operator import itemgetter
from scipy.stats import gaussian_kde

import numpy as np

def run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir,figuredir,genelist):    
    d = dict()
    with open(DMSOgenes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '0'
            d[gene] = [chrom,start,stop,number,strand,coverage]

    with open(DMSOTSS) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '0'
            d[gene].append(coverage)
            
    with open(DMSOEND) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '0'
            d[gene].append(coverage)
            
    with open(CAgenes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '0'
            d[gene].append(coverage)
            
    with open(CATSS) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '0'
            d[gene].append(coverage)
            
    with open(CAEND) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            if coverage == '.':
                coverage = '0'
            d[gene].append(coverage)
            
    coveragecutoff = 200
    TRlist = list()
    TRgenes = list()
    DMSOTRgenes = list()
    cutoff1 = 0.01
    ENDlist = list()
    ENDgenes = list()
    DMSOENDgenes = list()
    TRx = list()
    TRy = list()
    ERx = list()
    ERy = list()
    cutoff2 = 0.01
    cutoff3 = 0.25
    i = 0
    
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
        for name in genelist:
            if name == gene.split(';')[1]:
                print name,gene
                if CAgenes-CATSS > CATSS and DMSOgenes-DMSOTSS > DMSOTSS and CAgenes-CAEND > CAEND and DMSOgenes-DMSOEND > DMSOEND and DMSOgenes > coveragecutoff and CAgenes > coveragecutoff:
                    i += 1
                    TRy.append(CATSS/(CAgenes-CATSS))
                    TRx.append(DMSOTSS/(DMSOgenes-DMSOTSS))
                    ERy.append(CAEND/(CAgenes-CAEND))
                    ERx.append(DMSOEND/(DMSOgenes-DMSOEND))
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
    print "Genes: ",i
    
    
    
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
    F3 = plt.figure()
    ax1 = F3.add_subplot(111)
    xy = np.vstack([TRx,TRy])
    z = gaussian_kde(xy)(xy)
    ax1.scatter(TRx,TRy,c=z,edgecolor="",s=14)
    ax1.set_title('Travelers Ratio')
    ax1.set_ylabel('CA')
    ax1.set_xlabel('DMSO')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1])
    ax1.plot([0,1],[0,1],color='k')
    plt.savefig(figuredir + '/TravelersRatio_reflected.png')
    F5 = plt.figure()
    ax2 = F5.add_subplot(111)
    xy = np.vstack([ERx,ERy])
    z = gaussian_kde(xy)(xy)
    ax2.scatter(ERx,ERy,c=z,edgecolor="",s=14)
    ax2.set_title('End Ratio')
    ax2.set_ylabel('CA')
    ax2.set_xlabel('DMSO')
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1])
    ax2.plot([0,1],[0,1],color = 'k')
    plt.savefig(figuredir + '/EndRatio_reflected.png')
    F4 = plt.figure()
    ax1 = F4.add_subplot(111)
    bp1 = ax1.boxplot([TRlist,ENDlist],patch_artist=True)
    ax1.set_xticklabels(['Travelers Ratio','End Ratio'])
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ## change outline color, fill color and linewidth of the boxes
    for box in bp1['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )
    ## change color and linewidth of the whiskers
    for whisker in bp1['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    ## change color and linewidth of the caps
    for cap in bp1['caps']:
        cap.set(color='#7570b3', linewidth=2)
    ## change color and linewidth of the medians
    for median in bp1['medians']:
        median.set(color='#b2df8a', linewidth=2)
    ## change the style of fliers and their fill
    for flier in bp1['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    plt.savefig(figuredir + '/Boxplot.png')
        
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
    
    
    
    