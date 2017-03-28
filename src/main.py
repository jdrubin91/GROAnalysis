__author__ = 'Jonathan Rubin'

import sys
import os
import closest_gene
import create_annotations
import reflect_coverage
import bedtools_create_intersects
import master_writer
import pausing_index_boxplots

# #Specify DMSO treated bedgraph directory
# DMSO = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'

# #Specify CA treated bedgraph directory
# CA = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'

# #Specify gene annotations
# genes = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/refGene.sorted.bed'

# #ChIP directory to get closest genes - currently set at SRF ChIP dataset
# chipdir = '/scratch/Shares/dowell/ENCODE/old/HCT116/SRF/peak_files/ENCFF001UEM.bed'

#Specify list of genes to examine
#genelist = ['CYR61','NR4A3','FOS','ATF3','EGR2','FOSB','JUN','NR4A1','DUSP1','NR4A2','DUSP2','EGR3','BTG2','WEE1','THBS1','ZFP36','SNF1LK','EGR1','JUNB','BHLHB2','AXUD1','PTG52','IER2','DUSP5','PLK2','GEM','GDF15','KLF6','SNORD102']

#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

#Home directory
homedir = os.path.dirname(os.path.realpath(__file__))

#File directory
filedir = parent_dir(homedir) + '/files'

#Figure directory
figuredir = parent_dir(homedir) + '/figures'


def run(DMSO,CA):
    # #Specify DMSO treated bedgraph directory
    # DMSO = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J32_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'

    # #Specify CA treated bedgraph directory
    # CA = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J42_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    
    #Specify gene annotations
    genes = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/refGene.sorted.bed'

    #ChIP directory to get closest genes - currently set at SRF ChIP dataset
    chipdir = '/scratch/Shares/dowell/ENCODE/old/HCT116/SRF/peak_files/ENCFF001UEM.bed'

    #Specify list of genes to examine
    #genelist = ['CYR61','NR4A3','FOS','ATF3','EGR2','FOSB','JUN','NR4A1','DUSP1','NR4A2','DUSP2','EGR3','BTG2','WEE1','THBS1','ZFP36','SNF1LK','EGR1','JUNB','BHLHB2','AXUD1','PTG52','IER2','DUSP5','PLK2','GEM','GDF15','KLF6','SNORD102']

    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files'

    #Figure directory
    figuredir = parent_dir(homedir) + '/figures'
    
    # print "Calculating closest genes to ChIP file..."
    # genelist = closest_gene.run(chipdir,genes,filedir)
    # genes = genelist


    print "Creating annotation files..."
    create_annotations.run(genes,filedir)
    TSS = filedir + '/TSS.sorted.bed'
    END = filedir + '/END.sorted.bed'
    print "done\nReflecting coverage values..."
    reflect_coverage.run(DMSO,CA,filedir)
    print "done\nCreating intersect files..."
    DMSOreflect = filedir + '/DMSO.bedgraph'
    CAreflect = filedir + '/CA.bedgraph'
    bedtools_create_intersects.run(DMSOreflect,CAreflect,genes,TSS,END,filedir)
    os.system("rm " + filedir + "/DMSO.bedgraph")
    os.system("rm " + filedir + "/CA.bedgraph")
    print "done\nGenerating files..."
    DMSOgenes = filedir + '/DMSO.genes.bed'
    DMSOTSS = filedir + '/DMSO.TSS.bed'
    DMSOEND = filedir + '/DMSO.END.bed'
    CAgenes = filedir + '/CA.genes.bed'
    CATSS = filedir + '/CA.TSS.bed'
    CAEND = filedir + '/CA.END.bed'
    d1 = master_writer.run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir,figuredir)
    return d1
    print "done"

if __name__ == "__main__":
    DMSO = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J12_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    CA = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J22_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    d1 = run(DMSO,CA)

    DMSO = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J32_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    CA = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J42_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    d2 = run(DMSO,CA)

    DMSO = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J52_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    CA = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/J62_trimmed.flip.fastq.bowtie2.sorted.BedGraph.mp.BedGraph'
    d3 = run(DMSO,CA)



    pausing_index_boxplots.run(d1,d2,d3,figuredir)

    outfile = open(filedir + '/all_timepoints.txt','w')
    outfile.write(d1)
    outfile.write(d2)
    outfile.write(d3)
    
    