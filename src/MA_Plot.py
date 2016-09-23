__author__ = 'Jonathan Rubin'

# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import norm
import numpy as np
import math

def plot(file1,file2,file3):
    d1 = dict()
    with open(file1) as F1:
        i=0
        for line in F1:
            i+=1
            if '#' in line[0] and i > 20:
                break
            elif '#' not in line[0]:
                line = line.strip().split()
                d1[line[0]] = [float(line[1].split(',')[0]),float(line[2].split(',')[0])]
                    
    with open(file2) as F2:
        i=0
        for line in F2:
            i+=1
            if '#' in line[0] and i > 20:
                break
            elif '#' not in line[0]:
                line = line.strip().split()
                d1[line[0]].append(float(line[1].split(',')[0]))
                d1[line[0]].append(float(line[2].split(',')[0]))   

    silaclist=list()
    with open(file3) as F3:
        header=F3.readline().strip().split('\t')
        geneindex=header.index('Gene names')
        for line in F3:
            silaclist.append(line.strip().split('\t')[geneindex].split(';')[0])
    masterTFlist=['ID1','SREBF1','ONECUT3','ELF3','SMAD3','UBTF','TEAD3','ZIC5','ZBTB7B','MNT','JUNB','NR2F1','MYC','FOSB','HOXB8','ID3','NFIC','TGIF1','RARA','SP1','ZNF213','GLIS2','BHLHE40','TCF7L2','HES4','RARG','SOX12']

    X = list()
    Y = list()
    X2 = list()
    Y2 = list()
    diff = list()
    siglist = list()
    genelist = list()
    alpha = 0.1
    for key in d1:
        mdj=d1[key][1]
        mdk=d1[key][3]
        diff.append(mdj-mdk)
    mean = sum(diff)/len(diff)
    for key in d1:
        mdj=d1[key][1]
        mdk=d1[key][3]
        Nj=d1[key][0]
        Nk=d1[key][2]
        if (Nj+Nk)/2.0 > 10:
            p=((mdj*Nj)+(mdk*Nk))/(Nj+Nk)
            SE=(p*(1-p))/((1/Nj)+(1/Nk))
            Y.append(mdj-mdk-mean)
            X.append(math.log((Nj+Nk)/2.0,10))
            genelist.append(key)
            cdf=norm.cdf((mdj-mdk-mean)/math.sqrt(SE),0,1)
            p=min(cdf,1-cdf)
            if key.split('.')[0].split('_')[0] in masterTFlist:
                X2.append(math.log((Nj+Nk)/2.0,10))
                Y2.append(mdj-mdk-mean)
                siglist.append(key.split('.')[0].split('_')[0])

    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(X,Y,c='b',edgecolor="",s=14) 
    plt.scatter(X2,Y2,c='r',edgecolor="",s=14)
    ax.set_title('MD Scores MA Plot')
    ax.set_ylabel('MD Score Difference (CA-DMSO)')
    ax.set_xlabel('Mean Overlap Events (log10)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    for l in range(len(X2)):
        if 3 < X2[l] < 5:
            textx = 30
            texty = 30
            if Y2[l] < 0:
                texty=-30
            if siglist[l] == 'MNT':
                textx=-10
            if siglist[l] == 'RARG':
                textx=-7
                texty=70
            if siglist[l] == 'RARA':
                textx=20
                texty=60
            if siglist[l] == 'FOSB':
                texty=-50
                textx=0
            if siglist[l] == 'JUNB':
                texty=-40
            if siglist[l] == 'SP1':
                texty=-30
            ax.annotate(siglist[l], xy=(X2[l],Y2[l]), xytext=(textx, texty), ha='right',
                textcoords='offset points', 
                arrowprops=dict(arrowstyle='->', shrinkA=0))
        if 2 < X2[l] < 2.5 and Y2[l] < -0.06:
            textx = 30
            texty = -40
            ax.annotate(siglist[l], xy=(X2[l],Y2[l]), xytext=(textx, texty), ha='right',
                textcoords='offset points', 
                arrowprops=dict(arrowstyle='->', shrinkA=0))
    #     if 2.5 < X2[l] < 3 and (0.06 < Y2[l] or Y2[l] < -0.08):
    #         textx = 40
    #         texty = 50
    #         if siglist[l] == 'TFAP4':
    #             texty=30
    #             textx=60
    #         ax.annotate(siglist[l], xy=(X2[l],Y2[l]), xytext=(textx, texty), ha='right',
    #             textcoords='offset points', 
    #             arrowprops=dict(arrowstyle='->', shrinkA=0))
    #     if 3.5 < X2[l] < 4 and (0.06 < Y2[l] or Y2[l] < -0.02):
    #         textx = 40
    #         texty = -20
    #         ax.annotate(siglist[l], xy=(X2[l],Y2[l]), xytext=(textx, texty), ha='right',
    #             textcoords='offset points', 
                # arrowprops=dict(arrowstyle='->', shrinkA=0))
    plt.show()
    # plt.savefig(savedir + 'tsv_fig.png')
    
    return "done"

if __name__ == "__main__":
    file1='/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/motif_hits_2/Jonathan_CA_enrichment_stats.tsv'
    file2='/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/motif_hits_2/Jonathan_DMSO_enrichment_stats.tsv'
    file3='/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/HCT116_Serinduction_phospho_rep1.sorted.txt'
    plot(file1,file2,file3)
