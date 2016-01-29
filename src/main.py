__author__ = 'Jonathan Rubin'

import sys
import os
import create_annotations
import reflect_coverage
import bedtools_create_intersects
import master_writer

#Specify DMSO treated bedgraph directory
DMSO = sys.argv[1]

#Specify CA treated bedgraph directory
CA = sys.argv[2]

#Specify gene annotations
genes = sys.argv[3]

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

#Directories to reference files


def run():
    print "Creating annotation files..."
    create_annotations.run(genes,filedir)
    TSS = filedir + '/TSS.bed'
    END = filedir + '/END.bed'
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
    master_writer.run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir,figuredir)
    print "done"
    
    