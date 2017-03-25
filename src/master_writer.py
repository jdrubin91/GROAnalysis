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

def run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir,figuredir): 

    #Populate dictionary d with each key an annotation with coverage values
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
            
    
    #Initiate all required lists to store information
    coveragecutoff = 20
    TRlist = list()
    ENDlist = list()
    TRx = list()
    TRy = list()
    ERx = list()
    ERy = list()
    names = list()
    PIbarplot = list()
    Txnbarplot = list()
    namelist = list()
    pX = list()
    pY = list()
    pNames = list()
    expressionlist = list()
    cdf = list()
    i = 0
    
    #Initiate a Master file to store all pertinent information
    outfile = open(filedir + '/Master.bed','w')
    outfile.write('Gene\tChrom\tStart\tStop\tNumber\tStrand\tDMSO gene body\tDMSO TSS\tDMSO END\tCA gene body\tCA TSS\tCA END\n')
    
    #Iterate through dictionary d
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
        graphcutoff = 100
        name = gene.split(';')[1]
        
        #Populate lists to perform pearsons coefficient test for transcription across all gene bodies that are above a cutoff
        if (DMSOgenes+CAgenes)/2 > 1:
            pX.append(np.log(DMSOgenes))
            pY.append(np.log(CAgenes))
            pNames.append(name)
            
        #Populate lists for loci of interest
        if gene in ['NM_005252;FOS;chr14:75745480-75748937_+','NM_001964;EGR1;chr5:137801180-137805004_+','NM_001136177;EGR2;chr10:64571755-64576126_-','NM_004430;EGR3;chr8:22545173-22550815_-','NM_006981;NR4A3;chr9:102584136-102629173_+']:
            PIbarplot.append(np.log2((CATSS/(CAgenes-CATSS))/(DMSOTSS/(DMSOgenes-DMSOTSS))))
            Txnbarplot.append(np.log2(CAgenes/DMSOgenes))
            namelist.append(name)
        
        #Populate lists for traveler ratio (pausing index) and end ratio
        if CAgenes-CATSS > 0 and DMSOgenes-DMSOTSS > 0 and CAgenes-CAEND > 0 and DMSOgenes-DMSOEND > 0 and DMSOgenes > coveragecutoff and CAgenes > coveragecutoff and CATSS/(CAgenes-CATSS) < graphcutoff and DMSOTSS/(DMSOgenes-DMSOTSS) < graphcutoff and CAEND/(CAgenes-CAEND) < graphcutoff and DMSOEND/(DMSOgenes-DMSOEND) < graphcutoff:
            i += 1
            TRy.append(CATSS/(CAgenes-CATSS))
            TRx.append(DMSOTSS/(DMSOgenes-DMSOTSS))
            ERy.append(CAEND/(CAgenes-CAEND))
            ERx.append(DMSOEND/(DMSOgenes-DMSOEND))
            expressionlist.append((np.log2(DMSOgenes)+np.log2(CAgenes))/2.0)
            TR = (CATSS/(CAgenes-CATSS))-(DMSOTSS/(DMSOgenes-DMSOTSS))
            cdf.append(TR)
            TRlist.append(TR)
            # ENDlist.append(CAEND/(CAgenes-CAEND))-(DMSOEND/(DMSOgenes-DMSOEND))
            names.append(gene.split(';')[1])
            
    print "Genes: ",i
    
    #Perform Pearson coefficient test for all genes
    meanX = np.mean(pX)
    meanY = np.mean(pY)
    num = 0.0
    den1 = 0.0
    den2 = 0.0
    for i in range(len(pX)):
        X = pX[i]
        Y = pY[i]
        num += ((X - meanX)*(Y - meanY))
        den1 += (X - meanX)**2
        den2 += (Y-meanY)**2
    pearsons = num/(np.sqrt(den1)*np.sqrt(den2))
    
    #Determine distance and direction of travel ratio for each gene
    distance = list()
    direction1 = list()
    for i in range(len(TRx)):
        x = TRx[i]
        y = TRy[i]
        if y > x:
            direction1.append(1)
        else:
            direction1.append(0)
        xy = ((x+y)/2,(x+y)/2)
        d = np.sqrt((x-xy[0])**2+(y-xy[1])**2)
        #print x,y,xy,d
        distance.append(d)
    
    
    print "Genes up: ",sum(direction1)
    print "Genes down: ",len(direction1) - sum(direction1)

    #Populate list of significantly higher (>3SD from mean) and lower TR in CA
    N = len(distance)
    s = np.sqrt(np.var(distance))
    TRx2 = list()
    TRy2 = list()
    pX2 = list()
    pY2 = list()
    TRgenesup = list()
    TRgenesdwn = list()
    expressionlist2 = list()
    for i in range(N):
        if distance[i] > 3*s:
            TRx2.append(TRx[i])
            TRy2.append(TRy[i])
            expressionlist2.append(expressionlist[i])
            if direction1[i] == 1:
                TRgenesup.append(names[i])
                index = pNames.index(names[i])
                pX2.append(pX[index])
                pY2.append(pY[index])
            else:
                TRgenesdwn.append(names[i])
    
    
    #Calculate distance and direction for ER
    distance2 = list()
    direction = list()
    for i in range(len(ERx)):
        x = ERx[i]
        y = ERy[i]
        if y > x:
            direction.append(1)
        else:
            direction.append(0)
        xy = ((x+y)/2,(x+y)/2)
        d = np.sqrt((x-xy[0])**2+(y-xy[1])**2)
        distance2.append(d)
        
    #Populate list of significantly different ER genes
    N = len(distance2)
    s = np.sqrt(np.var(distance2))
    ERx2 = list()
    ERy2 = list()
    ERgenesup = list()
    ERgenesdwn = list()
    for i in range(N):
        if distance2[i] > 3*s:
            ERx2.append(ERx[i])
            ERy2.append(ERy[i])
            if direction[i] == 1:
                ERgenesup.append(names[i])
            else:
                ERgenesdwn.append(names[i])
    
    
#Preliminary Code to apply the Grubbs test for outliers to call significantly different TR and ER genes
    #TRx2 = list()
    #TRy2 = list()
    #grubbs = True
    #alpha = 0.05
    #distancelist = distance
    #i=0
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
    #        index2 = np.where(distancelist==M)[0][0]
    #        TRx2.append(TRx[index2])
    #        TRy2.append(TRy[index2])
    #        distance = np.delete(distance,distance[index])
    #        
    #ERx2 = list()
    #ERy2 = list()
    #grubbs = True
    #alpha = 0.05
    #distancelist2 = distance2
    #i=0
    #while grubbs == True:
    #    i+=1
    #    #print i
    #    N = len(distance2)
    #    M = np.amax(distance2)
    #    t = stats.t.ppf(1-alpha/2, N-1)
    #    s = np.sqrt(np.var(distance2))
    #    mean = np.mean(distance2)
    #    G = np.absolute(M-mean)/s
    #    grubbs = G > ((N-1)/np.sqrt(N))*np.sqrt((t**2)/(N-2+t**2))
    #    #print N,M,t,s,mean,G,grubbs
    #    #print G,((N-1)/np.sqrt(N))*np.sqrt((t**2)/(N-2+t**2))
    #    if grubbs == True:
    #        index = np.where(distance2==M)[0][0]
    #        index2 = np.where(distancelist2==M)[0][0]
    #        ERx2.append(ERx[index2])
    #        ERy2.append(ERy[index2])
    #        distance2 = np.delete(distance2,distance2[index])
    
    
#Code to determine slope and intercept of linear regression lines for empirical data
    slope1,intercept1 = np.polyfit(TRx, TRy, 1)
    slope2,intercept2 = np.polyfit(ERx, ERy, 1)
    
    #Histogram of TRCA - TRDMSO
    F1 = plt.figure()
    TRlist.sort(reverse=True)
    #plt.hist(TRlist[int(len(TRlist)*.2):int(len(TRlist)*.8)],50)
    plt.hist(TRlist,50)
    plt.title("Travelers Ratio")
    plt.savefig(figuredir + '/TravelersRatio.png')
    
    # #Histogram of ERCA - ERDMSO
    # F2 = plt.figure()
    # ENDlist.sort(reverse=True)
    # #plt.hist(ENDlist[int(len(ENDlist)*.2):int(len(ENDlist)*.8)],50)
    # plt.hist(ENDlist,50)
    # plt.title("End Ratio")
    # plt.savefig(figuredir + '/EndRatio.png')
    
    #Scatter plot of TR and ER with significant genes colored red
    F3 = plt.figure()
    ax1 = F3.add_subplot(121)
    xy = np.vstack([TRx,TRy])
    z = gaussian_kde(xy)(xy)
    ax1.scatter(TRx,TRy,c=z,edgecolor="",s=14)
    ax1.scatter(TRx2,TRy2,c='red',edgecolor="",s=14)
    ax1.set_title('Pausing Index')
    ax1.set_ylabel('CA')
    ax1.set_xlabel('DMSO')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    #ax1.plot([0,1/slope1],[intercept1,1],color = 'r')
    ax1.set_xlim([0, 20])
    ax1.set_ylim([0, 20])
    ax1.plot([0,50.0],[0,50.0],color='k')
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
    ax2.set_xlim([0, 20])
    ax2.set_ylim([0, 20])
    #ax2.plot([0,1/slope2],[intercept2,1],color = 'r')
    ax2.plot([0,50],[0,50],color = 'k')
    plt.savefig(figuredir + '/Scatter_reflected_moregenes.png')
    
    # #Boxplot of TRCA-TRDMSO and ERCA-ERDMSO
    # F4 = plt.figure()
    # ax1 = F4.add_subplot(111)
    # bp1 = ax1.boxplot([TRlist,ENDlist],patch_artist=True)
    # ax1.set_xticklabels(['Travel Ratio','End Ratio'])
    # ax1.get_xaxis().tick_bottom()
    # ax1.get_yaxis().tick_left()
    # ## change outline color, fill color and linewidth of the boxes
    # for box in bp1['boxes']:
    #     # change outline color
    #     box.set( color='#7570b3', linewidth=2)
    #     # change fill color
    #     box.set( facecolor = '#1b9e77' )
    # ## change color and linewidth of the whiskers
    # for whisker in bp1['whiskers']:
    #     whisker.set(color='#7570b3', linewidth=2)
    # ## change color and linewidth of the caps
    # for cap in bp1['caps']:
    #     cap.set(color='#7570b3', linewidth=2)
    # ## change color and linewidth of the medians
    # for median in bp1['medians']:
    #     median.set(color='#b2df8a', linewidth=2)
    # ## change the style of fliers and their fill
    # for flier in bp1['fliers']:
    #     flier.set(marker='o', color='#e7298a', alpha=0.5)
    # plt.savefig(figuredir + '/Boxplot.png')
    
    #Generate file that has significantly different genes
    outfile2 = open(filedir + '/GeneList.txt','w')
    outfile2.write("High CA TR = " + str(len(TRgenesup)) + "\nHigh DMSO TR = " + str(len(TRgenesdwn)) + "\nHigh CA ER = " + str(len(ERgenesup)) + "\nHigh DMSO ER = " + str(len(ERgenesdwn)) + "\n")
    outfile2.write("High CA TR\n")
    for item in sorted(TRgenesup, key=itemgetter(1),reverse=True):
        outfile2.write(item + '\t' + str(item[1]) + '\n')
    outfile2.write("High DMSO TR\n")
    for item in sorted(TRgenesdwn, key=itemgetter(1),reverse=True):
        outfile2.write(item + '\t' + str(item[1]) + '\n')
    outfile2.write("High CA ER\n")
    for item in sorted(ERgenesup, key=itemgetter(1),reverse=True):
        outfile2.write(item + '\t' + str(item[1]) + '\n')
    outfile2.write("High DMSO ER\n")
    for item in sorted(ERgenesdwn, key=itemgetter(1),reverse=True):
        outfile2.write(item + '\t' + str(item[1]) + '\n')
    
    #Generate pearson plot for transcription of all genes
    F5 = plt.figure()
    ax = F5.add_subplot(111)
    xy = np.vstack([pX,pY])
    z = gaussian_kde(xy)(xy)
    ax.scatter(pX,pY,c=z,edgecolor="",s=14)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_title('Gene Transcription')
    ax.set_xlabel('Log DMSO (FPKM)')
    ax.set_ylabel('Log CA (FPKM)')
    ax.set_xlim([0, 20])
    ax.set_ylim([0, 20])
    ax.plot([0,50000.0],[0,50000.0],color='k')
    ax.text(14,12, "Pearson = " + str(pearsons)[0:5])
    plt.savefig(figuredir + '/Pearson.png')
    
    #Generate pearson's coefficient for TR
    meanX = np.mean(TRx)
    meanY = np.mean(TRy)
    num = 0.0
    den1 = 0.0
    den2 = 0.0
    for i in range(len(TRx)):
        X = TRx[i]
        Y = TRy[i]
        num += ((X - meanX)*(Y - meanY))
        den1 += (X - meanX)**2
        den2 += (Y-meanY)**2
    pearsons = num/(np.sqrt(den1)*np.sqrt(den2))
    
    #Plot TR of all genes, include pearson coefficient and cdf plot
    F6 = plt.figure()
    ax1 = F6.add_subplot(111)
    xy = np.vstack([TRx,TRy])
    z = gaussian_kde(xy)(xy)
    ax1.scatter(TRx,TRy,c=z,edgecolor="",s=expressionlist)
    ax1.scatter(TRx2,TRy2,c='red',edgecolor="",s=expressionlist2)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_title('Pausing Index')
    ax1.set_ylabel('CA')
    ax1.set_xlabel('DMSO')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    #ax1.plot([0,1/slope1],[intercept1,1],color = 'r')
    # ax1.set_xlim([0, 20])
    # ax1.set_ylim([0, 20])
    ax1.plot([0,50.0],[0,50.0],color='k')
    ax1.text(8,18, "Pearson = " + str(pearsons)[0:5])
    # ax2 = F6.add_subplot(122)
    # ax2.plot(np.sort(cdf),np.linspace(0,1,len(cdf)))
    # ax2.plot(stats.norm.cdf(np.linspace(min(cdf),max(cdf)),0,np.var(cdf)),np.linspace(0,1,len(cdf)))
    plt.savefig(figuredir + '/PausingIndex.png')
    
    
    #Generate loci-specific plot
    order = ['FOS','EGR1','EGR2','EGR3','NR4A3']
    PIbarplotsorted = list()
    Txnbarplotsorted = list()
    for item in order:
        PIbarplotsorted.append(PIbarplot[namelist.index(item)])
        Txnbarplotsorted.append(Txnbarplot[namelist.index(item)])
    N = len(order)
    ind = np.arange(N)
    width = 0.3
    fig,ax1 = plt.subplots()
    PI = ax1.bar(ind,PIbarplotsorted,width, color='b')
    Txn = ax1.bar(ind+width,Txnbarplotsorted,width,color='r')
    ax1.plot(ax1.get_xlim(),[0,0],color = 'k')
    ax1.legend((PI,Txn),('Pausing Index','Transcription'))
    ax1.set_ylabel('Log$_2$ Fold Change (CA/DMSO)')
    ax1.set_xticks(ind + width)
    ax1.set_xticklabels(order)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    plt.savefig(figuredir + '/BarPlot.png')
    
    #Generate plot with signficantly different genes
    meanX = np.mean(pX2)
    meanY = np.mean(pY2)
    num = 0.0
    den1 = 0.0
    den2 = 0.0
    for i in range(len(pX2)):
        X = pX2[i]
        Y = pY2[i]
        num += ((X - meanX)*(Y - meanY))
        den1 += (X - meanX)**2
        den2 += (Y-meanY)**2
    pearsons = num/(np.sqrt(den1)*np.sqrt(den2))
    F7 = plt.figure()
    ax = F7.add_subplot(111)
    xy = np.vstack([pX2,pY2])
    z = gaussian_kde(xy)(xy)
    ax.scatter(pX2,pY2,c=z,edgecolor="",s=14)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlim([3, 7])
    ax.set_ylim([3, 7])
    ax.set_title('Gene Transcription')
    ax.set_xlabel('Log DMSO (FPKM)')
    ax.set_ylabel('Log CA (FPKM)')
    ax.plot([0,50000.0],[0,50000.0],color='k')
    ax.text(5.5,4, "Pearson = " + str(pearsons)[0:5])
    plt.savefig(figuredir + '/Transcription_UPGenes.png')
    
    
    