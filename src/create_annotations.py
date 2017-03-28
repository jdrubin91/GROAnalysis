__author__ = 'Jonathan Rubin'

import os

def run(genes,filedir):
    TSS = (0,1000)
    END = (0,1000)
    list1 = list()
    list2 = list()
    with open(genes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand = line.strip().split()
            if strand is '+':
                list1.append((chrom,int(start) + TSS[0],int(start) + TSS[1],gene,number,strand))
                list2.append((chrom,int(stop) + END[0],int(stop) + END[1],gene,number,strand))
            else:
                list2.append((chrom,int(start) + TSS[0],int(start) + TSS[1],gene,number,strand))
                list1.append((chrom,int(stop) + END[0],int(stop) + END[1],gene,number,strand))
    outfile = open(filedir + '/TSS.bed', 'w')
    for item in list1:
        outfile.write(item[0] + '\t' + str(item[1]) + '\t' + str(item[2]) + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')
    outfile2 = open(filedir + '/END.bed', 'w')
    for item in list2:
        outfile2.write(item[0] + '\t' + str(item[1]) + '\t' + str(item[2]) + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')
    outfile.close()
    outfile2.close()
    os.system("sort " + filedir + "/TSS.bed -k1,1 -k2,2n > " + filedir + "/TSS.sorted.bed")
    os.system("sort " + filedir + "/END.bed -k1,1 -k2,2n > " + filedir + "/END.sorted.bed")