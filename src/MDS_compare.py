__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def plot_MA(x,y,sig1,sig2,name,savedir):
    F = plt.figure() 
    ax = F.add_subplot(111)
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    plt.scatter(x,y,c=z,edgecolor="",s=14) 
    plt.scatter(sig1,sig2,c='r',edgecolor="",s=14)
    ax.set_title(name + ' MD Scores MA Plot')
    ax.set_ylabel('MD Score Difference (CA-DMSO)')
    ax.set_xlabel('Mean Overlap Events (log10)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.savefig(savedir + name + 'MA_plot.png')

def run(MDS1,MDS2,savedir):
    d = dict()
    with open(MDS1) as F:
        for line in F:
            if '#' not in line[0]:
                line = line.strip().split()
                d[line[0]] = [line[1].split(','),line[2].split(',')]

    with open(MDS2) as F:
        for line in F:
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
            diff.append(mdj-mdk)
        mean = sum(diff)/len(diff)
        for key in d:
            mdj=d[key][1][i]
            mdk=d[key][3][i]
            Nj=d[key][0][i]
            Nk=d[key][2][i]
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
        plot_MA(X,Y,X2,Y2,name,savedir)


if __name__ == "__main__":
    MDS1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J12_MDS.tsv'
    MDS2 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J32_MDS.tsv'
    savedir = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/'
    run(MDS1,MDS2,savedir)

    



