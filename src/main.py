__author__ = 'Jonathan Rubin'

import sys
import os
import reflect_coverage
import bedtools_create_intersects

#Specify DMSO treated bedgraph directory
DMSO = sys.argv[1]

#Specify CA treated bedgraph directory
CA = sys.argv[2]

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
genes = filedir + '/refGene.sorted.bed'
TSS = filedir + '/refTSS.sorted.bed'
END = filedir + '/ref3END.sorted.bed'

def run():
    print "Reflecting coverage values..."
    reflect_coverage.run(DMSO,CA,filedir)
    print "done\nCreating intersect files..."
    DMSOreflect = filedir + '/DMSO.bedgraph'
    CAreflect = filedir + '/CA.bedgraph'
    bedtools_create_intersects.run(DMSOreflect,CAreflect,genes,TSS,END,filedir)
    print "done"
    