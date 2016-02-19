__author__ = 'Jonathan Rubin'

import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
from operator import itemgetter
from scipy.stats import gaussian_kde
from scipy import stats
import numpy as np
import matplotlib as mpl

def run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir,figuredir):    
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
            
    coveragecutoff = 100
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
    
    distance = list()
    for i in range(len(TRx)):
        x = TRx[i]
        y = TRy[i]
        xy = ((x+y)/2,(x+y)/2)
        d = np.sqrt((x-xy[0])**2+(y-xy[1])**2)
        print x,y,xy,d
        distance.append(d)
        
    grubbs = True
    alpha = 0.05
    i=0
    N = len(distance)
    distancelist = np.zeros(N)
    s = np.sqrt(np.var(distance))
    TRx2 = list()
    TRy2 = list()
    print N,len(TRx),len(TRy)
    for i in range(N):
        if distance[i] > 3*s:
            print i
            TRx2.append(TRx[i])
            TRy2.append(TRy[i])
            
    distance2 = list()
    for i in range(len(ERx)):
        x = ERx[i]
        y = ERy[i]
        xy = ((x+y)/2,(x+y)/2)
        d = np.sqrt((x-xy[0])**2+(y-xy[1])**2)
        print x,y,xy,d
        distance2.append(d)
        
    grubbs = True
    alpha = 0.05
    i=0
    N = len(distance2)
    distancelist = np.zeros(N)
    s = np.sqrt(np.var(distance2))
    ERx2 = list()
    ERy2 = list()
    print N,len(ERx),len(ERy)
    for i in range(N):
        if distance2[i] > 3*s:
            print i
            ERx2.append(ERx[i])
            ERy2.append(ERy[i])
            
    #while grubbs == True:
    #    i+=1
    #    #print i
    #    N = len(distance)
    #    M = np.amax(distance)
    #    t = stats.t.ppf(1-alpha/2, N-1)
    #    s = np.sqrt(np.var(distance))
    #    mean = np.mean(distance)
    #    G = np.absolute(M-mean)/s
    #    grubbs = G > ((N-1)/np.sqrt(N))*np.sqrt((t**2)/(N-2+t**2))
    #    #print N,M,t,s,mean,G,grubbs
    #    #print G,((N-1)/np.sqrt(N))*np.sqrt((t**2)/(N-2+t**2))
    #    if grubbs == True:
    #        index = np.where(distance==M)[0][0]
    #        distancelist[index] = 1
    #        distance = np.delete(distance,distance[index])
    
    
    
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
    ax1 = F3.add_subplot(121)
    cmap, norm = mpl.colors.from_levels_and_colors(levels=levels, colors=colors, extend='max')
    xy = np.vstack([TRx,TRy])
    z = gaussian_kde(xy)(xy)
    ax1.scatter(TRx,TRy,c=z,edgecolor="",s=14)
    ax1.scatter(TRx2,TRy2,c='red',edgecolor="",s=14)
    ax1.set_title('Travelers Ratio')
    ax1.set_ylabel('CA')
    ax1.set_xlabel('DMSO')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1])
    ax1.plot([0,1],[0,1],color='k')
    ax2 = F3.add_subplot(122)
    xy = np.vstack([ERx,ERy])
    z = gaussian_kde(xy)(xy)
    ax2.scatter(ERx,ERy,c=z,edgecolor="",s=14)
    ax2.scatter(ERx2,ERy2,c='red',edgecolor="",s=14)
    ax2.set_title('End Ratio')
    ax2.set_ylabel('CA')
    ax2.set_xlabel('DMSO')
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1])
    ax2.plot([0,1],[0,1],color = 'k')
    plt.savefig(figuredir + '/Scatter_reflected_moregenes.png')
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
    
    
    
    