__author__ = "Jonathan Rubin"

# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import norm
import numpy as np
import math

#file1 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_DMSO-1_K_models_MLE.tsv'
#file2 = '/scratch/Shares/dowell/ENCODE/Rubin2016_genes_CA-1_K_models_MLE.tsv'
#savedir = '/scratch/Users/joru1876/GROAnalysis/figures/'
#genes = '/scratch/Users/joru1876/genome_files/refGene.bed'
# file1 = 'C:/cygwin64/home/Jonathan/Rubin2016_genes_DMSO-1_K_models_MLE.tsv'
# file2 = 'C:/cygwin64/home/Jonathan/Rubin2016_genes_CA-1_K_models_MLE.tsv'
# savedir = 'C:/cygwin64/home/Jonathan/GROAnalysis/figures/'
# genes = 'C:/cygwin64/home/Jonathan/refGene.bed'
file1 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/Rubin2016_genes_DMSO-1_K_models_MLE.tsv'
file2 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/Rubin2016_genes_CA-1_K_models_MLE.tsv'
savedir = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/GROAnalysis/figures/'
genes = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/refGene.bed'

index = 6
cut = 16
cut1 = 50000
cut2 = -50000

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def gene_dict(genes):
    genedict = dict()
    with open(genes) as F1:
        for line in F1:
            key,gene,position = line.strip().split()[3].split(';')[0:3]
            genedict[key] = (gene,position)
            
    return genedict

def run(file1,file2):
    d1 = dict()
    d2 = dict()
    with open(file1) as F1:
        for line in F1:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d1[gene] = [float(fwd), float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()
                    for param in p:
                        if '~' in param:
                            d1[gene].append(float(param.split(',')[1]))
                        elif ',' in param:
                            for moreparam in param.split(','):
                                d1[gene].append(float(moreparam))
                        else:
                            d1[gene].append(float(param))
                    
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d2[gene] = [float(fwd),float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()
                    for param in p:
                        if '~' in param:
                            d2[gene].append(float(param.split(',')[1]))
                        elif ',' in param:
                            for moreparam in param.split(','):
                                d2[gene].append(float(moreparam))
                        else:
                            d2[gene].append(float(param))
    
    
    for i in range(2,len(p[2:])+2):
        X = list()
        Y = list()
        x = list()
        y = list()
        for key in d1:
            if key in d2:
                if d1[key][0] > cut or d1[key][1] > cut and d2[key][0] > cut or d2[key][1] > cut:
                    if not (math.isinf(d1[key][i]) or math.isinf(d2[key][i]) or math.isnan(d1[key][i]) or math.isnan(d2[key][i])):
                        #if cut1 > d1[key][i] > cut2 and cut1 > d2[key][i] > cut2:
                        x.append(d1[key][i])
                        y.append(d2[key][i])
        #            if d1[key][2] != 0:
        #                if d2[key][2]-d1[key][2] > .25:
        #                    Y.append(key)
        #                else:
        #                    X.append(d2[key][2]-d1[key][2])
        #                
        #print "max: " + str(max(X))
        #print "min: " + str(min(X))
        #print "length: " + str(len(X))
        #print "avg: " + str(sum(X)/len(X))
        #print Y
        
        #plt.hist(X,50)
        print len(x)
        print len(y)
        F = plt.figure()
        plt.scatter(x,y,alpha=0.1)
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        plt.scatter(x,y,c=z,edgecolor="",s=14)
        plt.savefig(savedir + 'tsv_fig' + str(i) + '.png')
    
    return "done"
    
def run2(file1,file2):
    d1 = dict()
    d2 = dict()
    with open(file1) as F1:
        for line in F1:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d1[gene] = [float(fwd), float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()[index].split(',')[0]
                    d1[gene].append(float(p))
                    
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    fwd, rev = line[2].split(',')
                    d2[gene] = [float(fwd),float(rev)]
                elif '1' in line[1]:
                    p = line.strip().split()[index].split(',')[0]
                    d2[gene].append(float(p))
    
    
    X = list()
    Y = list()
    Z = list()
    for key in d1:
        if key in d2:
            if d1[key][0] > cut or d1[key][1] > cut and d2[key][0] > cut or d2[key][1] > cut:
                if d2[key][2]-d1[key][2] > .25:
                    #Y.append((key,d2[key][2]-d1[key][2]))
                    Y.append(key)
                if d2[key][2]-d1[key][2] < -.25:
                    #Z.append((key,d2[key][2]-d1[key][2]))
                    Z.append(key)
                X.append(d2[key][2]-d1[key][2])
    
    genedict = gene_dict(genes)
    Y1 = list()
    for item in Y:
        if item in genedict:
            Y1.append(genedict[item])
    Z1 = list()
    for item in Z:
        if item in genedict:
            Z1.append(genedict[item])
    
    
    print "max: " + str(max(X))
    print "min: " + str(min(X))
    print "length: " + str(len(X))
    print "avg: " + str(sum(X)/len(X))
    #print "Y: ",Y1
    print ','.join(Y1)
    print "================================================================="
    print ','.join(Z1)
    #for item in Z1:
    #    print item
    #print "Z: ",Z1
    #print "Y: ",sorted(Y, key=lambda x: x[1])
    #print "Z: ",sorted(Z, key=lambda x: x[1])
    print "len(Y): ",len(Y1)
    print "len(Z): ",len(Z1)
    F = plt.figure()        
    plt.hist(X,50)
    plt.savefig(savedir + 'tsv_fig.png')
    
    return "done"
    

def run3(file1,file2):
    d1 = dict()
    d2 = dict()
    with open(file1) as F1:
        for line in F1:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    position = line[1]
                    fwd, rev = line[2].split(',')
                    d1[gene] = [float(fwd), float(rev), position]
                elif '1' in line[1]:
                    p = line.strip().split()[index].split(',')[0]
                    d1[gene].append(float(p))
                    
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                if '>' in line[0]:
                    line = line.strip().split('|')
                    gene = line[0][1:]
                    position = line[1]
                    fwd, rev = line[2].split(',')
                    d2[gene] = [float(fwd),float(rev), position]
                elif '1' in line[1]:
                    p = line.strip().split()[index].split(',')[0]
                    d2[gene].append(float(p))
    
    
    X = list()
    Y = list()
    Z = list()
    genelist = list()
    for key in d1:
        if key in d2:
            if (d1[key][0] > cut or d1[key][1] > cut) and (d2[key][0] > cut or d2[key][1] > cut) and d2[key][3] > 0 and d1[key][3] > 0:
                #if d2[key][2]-d1[key][2] > .25:
                    #Y.append((key,d2[key][2]-d1[key][2]))
                    #Y.append(key)
                if d2[key][3]-d1[key][3] < -.25:
                    #Z.append((key,d2[key][2]-d1[key][2]))
                    Z.append(key)
                Y.append(math.log(d2[key][3]/d1[key][3],2))
                X.append(math.log((sum((d2[key][0],d2[key][1]))+sum((d1[key][0],d1[key][1])))/2.0,2))
                genelist.append(key)
                
    
    M = max(X)
    m = min(X)
    Z = list()
    size = ((M-m)/5.0)
    for i in range(int(int(M-m)/0.1)):
        window = (i*0.1+m,i*0.1+size+m)
        Z.append([])
        for j in range(len(X)):
            if X[j] > window[0] and X[j] <= window[1]:
                Z[i].append(Y[j])
    X1 = list()
    for item in Z:
        var = list()
        for val in item:
            var.append(val**2)
        if len(var) > 2:
            X1.append(sum(var)/(len(var)-1))
        else:
            X1.append(0.0)
    
    X2 = list()
    Y2 = list()
    genedict = gene_dict(genes)
    siglist = list()
    bedlist = list()
    for i in range(len(X)):
        if int((X[i] - m)*10) < len(X1):
            k = int((X[i] - m)*10)
            var = X1[k]
        else:
            k = len(X1)-1
            var = X1[k]
        cdf = norm.cdf(Y[i],0,math.sqrt(var))
        p = min(cdf,1-cdf)*len(Z[k])
        if p < 0.1:
            X2.append(X[i])
            Y2.append(Y[i])
            if genelist[i] in genedict:
                siglist.append(genedict[genelist[i]][0])
                bedlist.append(genedict[genelist[i]][1])
            else:
                siglist.append(genelist[i])
                bedlist.append(d1[genelist[i]][2])
            
    outfile = open('/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/GROAnalysis/files/tsv_analysis.bed','w')
    window = (200)
    for item in bedlist:
        item = item.split(':')
        chrom = item[0]
        start = str(int(item[1].split('-')[0])-window)
        stop = str(int(item[1].split('-')[0])+window)
        outfile.write(chrom + '\t' + start + '\t' + stop + '\n')
    F = plt.figure() 
    ax = F.add_subplot(111)
    xy = np.vstack([X,Y])
    z = gaussian_kde(xy)(xy)
    plt.scatter(X,Y,c=z,edgecolor="",alpha=0.5,s=14) 
    plt.scatter(X2,Y2,c='r',edgecolor="",s=14)
    ax.set_xlim([4,20])
    ax.set_title('Travel Ratio MA Plot')
    ax.set_ylabel('log$_2$(TR$_{CA}$/TR$_{DMSO}$)')
    ax.set_xlabel('log$_2$(Expression)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    for l in range(len(X2)):
        textx = 40
        texty = 15
        if X2[l] < 16 and Y2[l] < 0:
            textx = 50
            texty = -15
        if siglist[l] == 'NR4A2':
            texty = 15
        if siglist[l] == 'PMF1-BGLAP' or siglist[l] == 'C2orf16':
            texty = -25
        ax.annotate(siglist[l], xy=(X2[l],Y2[l]), xytext=(textx, texty), ha='right',
                textcoords='offset points', 
                arrowprops=dict(arrowstyle='->', shrinkA=0))
    #ax.set_xscale('log', basex=2)
    #ax.set_yscale('log', basey=2)
    #ax2 = F.add_subplot(212)
    #ax2.plot(X1)
    plt.show()
    plt.savefig(savedir + 'tsv_fig.png')
    
    return "done"

    
if __name__ == "__main__":
    run3(file1,file2)