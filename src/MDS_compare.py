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

def plot_MA(x,y,sig1,sig2,name,savedir,siglist):
    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(x,y,edgecolor="",s=14) 
    plt.scatter(sig1,sig2,c='r',edgecolor="",s=14)
    ax.set_title(name + ' MD Scores MA Plot')
    ax.set_ylabel('MD Score Difference (CA-DMSO)')
    ax.set_xlabel('Mean Overlap Events (log10)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    for l in range(len(siglist)):
        ax.annotate(siglist[l], xy=(sig1[l],sig2[l]), xytext=(textx, texty), ha='right',
                    textcoords='offset points', 
                    arrowprops=dict(arrowstyle='->', shrinkA=0))
    plt.show()
    plt.savefig(savedir + name + 'MA_plot.png')

def run(MDS1,MDS2,savedir):
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
        diff = list()
        siglist = list()
        genelist = list()
        alpha = 0.1
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
            if (Nj+Nk)/2.0 > 10:
                p=((mdj*Nj)+(mdk*Nk))/(Nj+Nk)
                SE=(p*(1-p))/((1/Nj)+(1/Nk))
                Y.append(mdj-mdk-mean)
                X.append(math.log((Nj+Nk)/2.0,10))
                genelist.append(key)
                # print SE
                if SE == 0:
                    print SE
                    cdf = norm.cdf(0,0,1)
                else:
                    print (mdj-mdk-mean)/math.sqrt(SE)
                    cdf=norm.cdf((mdj-mdk-mean)/math.sqrt(SE),0,1)
                p=min(cdf,1-cdf)
                if p < 0.05:
                    X2.append(math.log((Nj+Nk)/2.0,10))
                    Y2.append(mdj-mdk-mean)
                    siglist.append(key.split('.')[0].split('_')[0])
        plot_MA(X,Y,X2,Y2,name,savedir,siglist)


if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files'
    MDS1 = parent_dir(homedir) + '/MDS_files/J12_MDS.tsv'
    MDS2 = parent_dir(homedir) + '/MDS_files/J32_MDS.tsv'
    savedir = parent_dir(homedir) + '/figures'
    run(MDS1,MDS2,savedir)

    



