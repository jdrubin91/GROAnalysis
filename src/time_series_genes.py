__author__ = 'Jonathan Rubin'

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import operator

def run(folder,savedir):
    pvalcut = 0.01
    timepoints = ['15','45']
    i = 0
    d = dict()
    for file1 in folder:
        time = timepoints[i]
        i+=1
        with open(file1) as F:
            F.readline()
            for line in F:
                line = line.strip().split()
                try:
                    pval = float(line[-2])
                except ValueError:
                    pval = 1
                try:
                    log2change = float(line[-3])
                except ValueError:
                    log2change = 0
                gene = line[1]
                if gene not in d:
                    d[gene] = list()
                d[gene].append((time,log2change,pval))
    genes = list()
    x = list()
    y = list()
    g15up = list()
    x15up = list()
    y15up = list()
    g15do = list()
    x15do = list()
    y15do = list()
    g45up = list()
    x45up = list()
    y45up = list()
    g45do = list()
    x45do = list()
    y45do = list()
    pvals = list()
    for key in d:
        if len(d[key]) > 1:
            gene = key.split(';')[1]
            genes.append(gene)
            x.append(0)
            y.append(0)
            for i in range(len(d[key])):
                tup = d[key][i]
                if i==0:
                    if tup[1] > 0:
                        if tup[2] < pvalcut:
                            if gene not in g15up:
                                g15up.append(gene)
                                x15up.append(tup[0])
                                x15up.append(d[key][i+1][0])
                                y15up.append(tup[1])
                                y15up.append(d[key][i+1][1])
                    if tup[1] < 0:
                        if tup[2] < pvalcut:
                            if gene not in g15do:
                                g15do.append(gene)
                                x15do.append(tup[0])
                                x15do.append(d[key][i+1][0])
                                y15do.append(tup[1])
                                y15do.append(d[key][i+1][1])
                if i==1:
                    if tup[1] > 0:
                        if tup[2] < pvalcut:
                            if gene not in g15up and gene not in g15do:
                                g45up.append(gene)
                                x45up.append(tup[0])
                                x45up.append(d[key][i-1][0])
                                y45up.append(tup[1])
                                y45up.append(d[key][i-1][1])
                    if tup[1] < 0:
                        if tup[2] < pvalcut:
                            if gene not in g15up and gene not in g15do:
                                g45do.append(gene)
                                x45do.append(tup[0])
                                x45do.append(d[key][i-1][0])
                                y45do.append(tup[1])
                                y45do.append(d[key][i-1][1])


    for item in g45up:
        print item
    # F = plt.figure()
    # ax = F.add_subplot(111)
    # ax.set_title('Serum Induction Timecourse')
    # ax.set_ylabel('Log$_2$ Fold Change Relative to Starvation')
    # ax.set_xlabel('Time after Serum Induction (min)')
    # ax.get_xaxis().tick_bottom()
    # ax.get_yaxis().tick_left()
    # plt.scatter(x,y,c='k',edgecolor="",alpha=0.1)
    # plt.scatter(x15up,y15up,c='g',edgecolor="",alpha=0.1)
    # plt.scatter(x15do,y15do,c='r',edgecolor="",alpha=0.1)
    # plt.scatter(x45up,y45up,c='b',edgecolor="",alpha=0.1)
    # plt.scatter(x45do,y45do,c='m',edgecolor="",alpha=0.1)
    # for i in range(0,len(x15up),2):
    #     line1, = ax.plot([0,x15up[i]],[0,y15up[i]],color='g',alpha=0.1)
    #     ax.plot([x15up[i],x15up[i+1]],[y15up[i],y15up[i+1]],color='g',alpha=0.1)

    # for i in range(0,len(x15do),2):
    #     line2, = ax.plot([0,x15do[i]],[0,y15do[i]],color='r',alpha=0.1)
    #     ax.plot([x15do[i],x15do[i+1]],[y15do[i],y15do[i+1]],color='r',alpha=0.1)

    # for i in range(0,len(x45up),2):
    #     line3, = ax.plot([0,x45up[i+1]],[0,y45up[i+1]],color='b',alpha=0.1)
    #     ax.plot([x45up[i+1],x45up[i]],[y45up[i+1],y45up[i]],color='b',alpha=0.1)

    # for i in range(0,len(x45do),2):
    #     line4, = ax.plot([0,x45do[i+1]],[0,y45do[i+1]],color='m',alpha=0.1)
    #     ax.plot([x45do[i+1],x45do[i]],[y45do[i+1],y45do[i]],color='m',alpha=0.1)
    # green_patch = mpatches.Patch(color='green', label='Genes up at 15min')
    # red_patch = mpatches.Patch(color='red', label='Genes down at 15min')
    # blue_patch = mpatches.Patch(color='blue', label='Genes up at 45min')
    # magenta_patch = mpatches.Patch(color='magenta', label='Genes down at 45min')
    # ax.legend([green_patch,red_patch,blue_patch,magenta_patch],['Genes up at 15min','Genes down at 15min','Genes up at 45min','Genes down at 45min'],loc=3,fontsize=10)
    # F.savefig(savedir + 'serum_timecourse.png', dpi=1200)





if __name__ == "__main__":
    folder = ['/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/refGene.bed.count.bed.starvinducednascent.res.txt','/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/starvinduced45.genes.bed.count.bed.starvinduced45nascent.res.txt']
    savedir = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/'
    run(folder,savedir)
