__author__ = 'Jonathan Rubin'

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import operator
import math

def run(file1):
    outfile = open('/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/15min_long_genes.bed','w')
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
                start = int(gene.split(':')[-1].split('-')[0])
                stop = int(gene.split(':')[-1].split('-')[1].split('_')[0])
                chrom = gene.split(';')[2].split(':')[0]
                strand = gene.split('_')[1]
                length = math.fabs(stop-start)
                if pval < 0.01 and length > 100000:
                    if start > stop:
                        outfile.write('\t'.join([chrom,str(start),str(stop)]) +'\n')
                    else:
                        outfile.write('\t'.join([chrom,str(stop),str(start)]) +'\n')
                    print gene
                # if gene not in d:
                #     d[gene] = list()
                # d[gene].append((log2change,pval))


if __name__ == "__main__":
    file1 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/refGene.bed.count.bed.starvinducednascent.res.txt'
    run(file1)