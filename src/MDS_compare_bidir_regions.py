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

def plot_MA(x,y,sig1,sig2,sig3,sig4,savedir,siglist1,siglist2,name1,name2,genelist):
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
    # for i in range(len(siglist1)):
    #     if siglist1[i] == 'VDR':
    #         ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (5.5,0.05),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if siglist1[i] == 'RARG':
    #         ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (4.5,0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if siglist1[i] == 'RARB':
    #         ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (3.7,0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if siglist1[i] == 'RARA':s
    #         ax.annotate(siglist1[i],xy=(sig1[i],sig2[i]),xytext = (5.2,0.09),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    # for i in range(len(genelist)):
    #     if genelist[i] =='VDR' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (5.5,0.05),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if genelist[i] =='RARG' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (4.5,0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if genelist[i] =='RARB' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (3.7,0.1),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))
    #     if genelist[i] =='RARA' and x[i] > 3:
    #         ax.annotate(genelist[i],xy=(x[i],y[i]),xytext = (5.2,0.09),arrowprops=dict(facecolor='black', shrink=0.1,width = 1,headwidth=5))

    # print siglist1
    # print siglist2
    # for item in siglist1:
    #     print item
    # for item in siglist2:
    #     print item
    # plt.show()
    plt.savefig(savedir + 'MA_plot_' + name1 + '-' + name2 + '.png',dpi=1200)

def plot_MA_nosig(x,y,name1,name2,savedir):
    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(x,y,edgecolor="",s=14)
    ax.set_title('MA plot ' + name1 + '/' + name2)
    ax.set_ylabel('Log2 Fold Change (' + name1 + '/' + name2 + ')')
    ax.set_xlabel('Mean Expression (log10)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    # plt.show()
    plt.savefig(savedir + 'MA_plot_' + name1 + '-' + name2 + '_bidir_predictions.png',dpi=1200)

#Takes in a list of ints and returns h/H and H
def compute_MDS(histlist):
    windowsize = 150
    H = sum(histlist)
    middle = len(histlist)/2
    h = sum(histlist[middle-windowsize:middle]) + sum(histlist[middle:middle+windowsize])
    if H != 0:
        return float(h)/float(H), H
    else:
        return 0,0


def run(count_file,savedir):
    name1 = count_file.split('/')[-1].split('_')[0]
    name2 = count_file.split('/')[-1].split('_')[1]
    totals = [0,0]
    with open(count_file) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            totals[0] += float(line[-2])
            totals[1] += float(line[-1])

    norm = totals[0]/totals[1]
    x = list()
    y = list()
    with open(count_file) as F:
        F.readline()
        for line in F:
            line = line.strip().split()
            try:
                val1 = math.log((float(line[-2]) + float(line[-1])*norm)/2,2)
                val2 = math.log(float(line[-2])/(float(line[-1])*norm),10)
                x.append(val1)
                y.append(val2)
            except:
                print "Some Math error"

    plot_MA_nosig(x,y,name1,name2,savedir)

    

    # X = list()
    # Y = list()
    # X2 = list()
    # Y2 = list()
    # X3 = list()
    # Y3 = list()
    # diff = list()
    # siglist = list()
    # siglist2 = list()
    # genelist = list()
    # ps = list()
    # zs = list()
    # for key in d:
    #     mdj=d[key][0]
    #     mdk=d[key][2]
    #     diff.append(float(mdj)-float(mdk))
    # mean = sum(diff)/len(diff)
    # for key in d:
    #     mdj=float(d[key][0])
    #     mdk=float(d[key][2])
    #     Nj=float(d[key][1])
    #     Nk=float(d[key][3])
    #     if (Nj+Nk)/2.0 > 10:
    #         p=((mdj*Nj)+(mdk*Nk))/(Nj+Nk)
    #         SE=(p*(1-p))*((1/Nj)+(1/Nk))
    #         Y.append(mdj-mdk-mean)
    #         X.append(math.log((Nj+Nk)/2.0,10))
    #         genelist.append(key.split('.')[0].split('_')[1])
    #         try:
    #             z = (mdj-mdk-mean)/math.sqrt(SE)
    #         except ZeroDivisionError:
    #             z=0
    #         zs.append(z)
    #         cdf=norm.cdf(z,0,1)
    #         p=min(cdf,1-cdf)*2
    #         ps.append(p)
    #         if p < 0.01:
    #             if mdj-mdk-mean > 0:
    #                 X2.append(math.log((Nj+Nk)/2.0,10))
    #                 Y2.append(mdj-mdk-mean)
    #                 siglist.append(key.split('.')[0].split('_')[1])
    #                 print 'up', key.split('.')[0].split('_')[1], math.log((Nj+Nk)/2.0,10), mdj-mdk-mean
    #             else:
    #                 X3.append(math.log((Nj+Nk)/2.0,10))
    #                 Y3.append(mdj-mdk-mean)
    #                 siglist2.append(key.split('.')[0].split('_')[1])
    #                 print 'down', key.split('.')[0].split('_')[1], math.log((Nj+Nk)/2.0,10), mdj-mdk-mean
    # plot_MA(X,Y,X2,Y2,X3,Y3,savedir,siglist,siglist2,name1,name2,genelist)


if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    savedir = parent_dir(homedir) + '/figures/'

    count_file = filedir + 'ACN_ACD_bidir.bed.count.bed'
    run(count_file,savedir)