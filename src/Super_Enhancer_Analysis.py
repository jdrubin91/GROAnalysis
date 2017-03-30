__author__ = 'Jonathan Rubin'

from pybedtools import BedTool

def column_add(file1,file2,name):
    outfile = open(name,'w')
    with open(file1) as F:
        with open(file2) as F2:
            for line in F:
                line2 = F2.readline()
                line = line.strip() + '\t' + line2.strip().split()[-1] + '\n'
                outfile.write(line)

def run(bedgraph1,bedgraph2,SEs,figdir,filedir):
    b1name = bedgraph1.split('/')[-1]
    b2name = bedgraph2.split('/')[-1]
    b1 = BedTool(bedgraph1)
    b2 = BedTool(bedgraph2)
    s = BedTool(SEs).sort()
    s.map(b1,o="sum").saveas(filedir + b1name + '_SE.bed')
    s.map(b2,o="sum").saveas(filedir + b2name + '_SE.bed')
    column_add(filedir + b1name + '_SE.bed',filedir + b2name + '_SE.bed', filedir + 'SE_Counts.bed')


if __name__ == "__main__":
    SEs = '/scratch/Users/joru1876/HCT-116_Super_Enhancers.bed'
    figdir = '/scratch/Users/joru1876/GROAnalysis/figures/'
    filedir = '/scratch/Users/joru1876/GROAnalysis/files/'
    bedgraph1 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J52_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    bedgraph2 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J62_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    run (bedgraph1,bedgraph2,SEs,figdir,filedir)