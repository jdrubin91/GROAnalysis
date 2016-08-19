__author__ = 'Jonathan Rubin'

import os

#Create intersect temp files
def run(DMSO,CA,genes,TSS,END,filedir):
    os.system("bedtools map -a " + genes + " -b " + DMSO + " -c 4 -o sum > " + filedir + "/DMSO.genes.bed")
    os.system("bedtools map -a " + TSS + " -b " + DMSO + " -c 4 -o sum > " + filedir + "/DMSO.TSS.bed")
    os.system("bedtools map -a " + END + " -b " + DMSO + " -c 4 -o sum > " + filedir + "/DMSO.END.bed")
    os.system("bedtools map -a " + genes + " -b " + CA + " -c 4 -o sum > " + filedir + "/CA.genes.bed")
    os.system("bedtools map -a " + TSS + " -b " + CA + " -c 4 -o sum > " + filedir + "/CA.TSS.bed")
    os.system("bedtools map -a " + END + " -b " + CA + " -c 4 -o sum > " + filedir + "/CA.END.bed")
    