__author__ = 'Jonathan Rubin'

import operator

filedir = 'C:/cygwin64/home/Jonathan/GROAnalysis/files/'

def run(filedir):
    posgenes = list()
    neggenes = list()
    with open(filedir + 'Master.bed') as F:
        for line in F:
            if line[0] is 'N':
                gene,chrom,start,stop,number,strand,DMSOgene,exp2,exp3,CAgene,exp4,exp5 = line.strip().split()
                if strand is '+': 
                    posgenes.append((gene,(float(DMSOgene)+float(CAgene))/2))
                else:
                    neggenes.append((gene,(float(DMSOgene)+float(CAgene))/2))
    posgenes.sort(key=operator.itemgetter(1))
    neggenes.sort(key=operator.itemgetter(1))
    
    d = dict()
    outfile1 = open(filedir + 'train.pos.txt','w')
    for gene in posgenes[int(len(posgenes)*.85)-15:int(len(posgenes)*0.85)]:
        chrom, [start,stop] = gene[0].split(';')[2].split(':')[0],gene[0].split(';')[2].split(':')[1].split('_')[0].split('-')
        d[chrom + '\t' + start + '\t' + stop + '\t1' + '\n'] = list()
    for key in d:
        outfile1.write(key)
    
    d = dict()
    for gene in posgenes[0:15]:
        chrom, [start,stop] = gene[0].split(';')[2].split(':')[0],gene[0].split(';')[2].split(':')[1].split('_')[0].split('-')
        d[chrom + '\t' + start + '\t' + stop + '\t0' + '\n'] = list()
    for key in d:
        outfile1.write(key)
    
    outfile2 = open(filedir + 'train.neg.txt','w')
    d = dict()
    for gene in neggenes[int(len(neggenes)*.85)-15:int(len(neggenes)*0.85)]:
        chrom, [start,stop] = gene[0].split(';')[2].split(':')[0],gene[0].split(';')[2].split(':')[1].split('_')[0].split('-')
        d[chrom + '\t' + start + '\t' + stop + '\t1' + '\n'] = list()
    for key in d:
        outfile2.write(key)
    
    d = dict()
    for gene in neggenes[0:15]:
        chrom, [start,stop] = gene[0].split(';')[2].split(':')[0],gene[0].split(';')[2].split(':')[1].split('_')[0].split('-')
        d[chrom + '\t' + start + '\t' + stop + '\t0' + '\n'] = list()
    for key in d:
        outfile2.write(key)

if __name__ == "__main__":
    run(filedir)