__author__ = 'Jonathan Rubin'

import os
import math
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def plot_MA(x,y,sig1,sig2,sig3,sig4,name,savedir,siglist1,siglist2,name1,name2,genelist):
    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(x,y,edgecolor="",s=14) 
    plt.scatter(sig1,sig2,c='r',edgecolor="",s=14)
    plt.scatter(sig3,sig4,c='g',edgecolor="",s=14)
    ax.set_title('MD Scores ' + name1 + '-' + name2)
    ax.set_ylabel('MD Score Difference (' + name1 + '-' + name2 + ')')
    ax.set_xlabel('Mean Overlap Events (log10)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylim([-0.25,0.25])
    for i in range(len(siglist1)):
        if siglist1[i] == 'VDR':
            ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (4.6,0.06),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        if siglist1[i] == 'RARG':
            ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (3.7,0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        if siglist1[i] == 'RARB':
            ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (3.5,0.07),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        if siglist1[i] == 'RARA':
            ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (4.4,0.09),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    # for i in range(len(genelist)):
    #     if genelist[i] =='VDR' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (4.5,0.05),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if genelist[i] =='RARG' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (3.7,0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if genelist[i] =='RARB' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (3.5,0.07),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if genelist[i] =='RARA' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (4.4,0.09),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))

    print siglist1
    print siglist2
    for item in siglist1:
        print item
    for item in siglist2:
        print item
    # plt.show()
    plt.savefig(savedir + 'MA_plot_JDR_' + name1 + '-' + name2 + '_annotated.png',dpi=1200)



        # for i in range(len(siglist2)):
        #     if siglist2[i] == 'JUN':
        #         ax.annotate(siglist2[i],xy=(sig3[i],sig4[i]),xytext = (2.75,-0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        #     if siglist2[i] == 'JUND':
        #         ax.annotate(siglist2[i],xy=(sig3[i],sig4[i]),xytext = (3,-0.12),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        #     if siglist2[i] == 'FOSL1':
        #         ax.annotate(siglist2[i],xy=(sig3[i],sig4[i]),xytext = (3.5,-0.13),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        #     if siglist2[i] == 'FOSL2':
        #         ax.annotate(siglist2[i],xy=(sig3[i],sig4[i]),xytext = (3.6,-0.09),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        #     if siglist2[i] == 'MAFG':
        #         ax.annotate(siglist2[i],xy=(sig3[i],sig4[i]),xytext = (4,-0.06),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
        #     if siglist2[i] == 'BATF':
    #     #         ax.annotate(siglist2[i],xy=(sig3[i],sig4[i]),xytext = (2.5,-0.08),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    # if name == 'BOTH':
    #     for item in siglist1:
    #         print item
    #     for item in siglist2:
    #         print item
    # plt.show()
    # plt.savefig(savedir + name + 'MA_plot_t45_annotated.png',dpi=1200)

def run(MDS1,MDS2,savedir):
    name1 = MDS1.split('/')[-1].split('_')[0]
    name2 = MDS2.split('/')[-1].split('_')[0]
    d = dict()
    with open(MDS1) as F:
        for line in F:
            if '#Binned' in line:
                break
            if '#' not in line[0]:
                line = line.strip().split()
                d[line[0]] = [line[1].split(','),line[2].split(',')]

    with open(MDS2) as F:
        for line in F:
            if '#Binned' in line:
                break
            if '#' not in line[0]:
                line = line.strip().split()
                d[line[0]].append(line[1].split(','))
                d[line[0]].append(line[2].split(','))

    namelist = ['NON','TSS','BOTH']
    
    for i in range(len(namelist)):
        name = namelist[i]
        X = list()
        Y = list()
        X2 = list()
        Y2 = list()
        X3 = list()
        Y3 = list()
        diff = list()
        siglist = list()
        siglist2 = list()
        genelist = list()
        ps = list()
        zs = list()
        for key in d:
            mdj=d[key][1][i]
            mdk=d[key][3][i]
            diff.append(float(mdj)-float(mdk))
        mean = sum(diff)/len(diff)
        for key in d:
            mdj=float(d[key][1][i])
            mdk=float(d[key][3][i])
            Nj=float(d[key][0][i])
            Nk=float(d[key][2][i])
            # if key.split('.')[0].split('_')[1] == 'SRF':
            #     X2.append(mdj-mdk-mean)
            #     Y2.append(math.log((Nj+Nk)/2.0,10))
            if (Nj+Nk)/2.0 > 10:
                p=((mdj*Nj)+(mdk*Nk))/(Nj+Nk)
                SE=(p*(1-p))*((1/Nj)+(1/Nk))
                Y.append(mdj-mdk-mean)
                X.append(math.log((Nj+Nk)/2.0,10))
                genelist.append(key.split('.')[0].split('_')[1])
                try:
                    z = (mdj-mdk-mean)/math.sqrt(SE)
                except ZeroDivisionError:
                    z=0
                zs.append(z)
                cdf=norm.cdf(z,0,1)
                p=min(cdf,1-cdf)*2
                ps.append(p)
                if p < 0.01:
                    if mdj-mdk-mean > 0:
                        X2.append(math.log((Nj+Nk)/2.0,10))
                        Y2.append(mdj-mdk-mean)
                        siglist.append(key.split('.')[0].split('_')[1])
                        print 'up', key.split('.')[0].split('_')[1], math.log((Nj+Nk)/2.0,10), mdj-mdk-mean
                    else:
                        X3.append(math.log((Nj+Nk)/2.0,10))
                        Y3.append(mdj-mdk-mean)
                        siglist2.append(key.split('.')[0].split('_')[1])
                        print 'down', key.split('.')[0].split('_')[1], math.log((Nj+Nk)/2.0,10), mdj-mdk-mean
        if name == 'NON':
            plot_MA(X,Y,X2,Y2,X3,Y3,name,savedir,siglist,siglist2,name1,name2,genelist)

    # print genelist

if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    MDS1 = parent_dir(homedir) + '/MDS_files/A2D_MDS.tsv'
    MDS2 = parent_dir(homedir) + '/MDS_files/ACD_MDS.tsv'
    savedir = parent_dir(homedir) + '/figures/'
    run(MDS2,MDS1,savedir)

    



