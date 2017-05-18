__author__ = 'Jonathan Rubin'

import os
import sys
import pybedtools
from pybedtools import BedTool

def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

def convert_deseqgenes_to_tssbed(file1):
    bedfile = list()
    with open(file1) as F:
        for line in F:
            if "id" not in line[0] and "NA" not in line[0]:
                line = line[1].split(';')[-1]
                print line
                chrom = line.split(':')[0]
                if line.split('_')[-1] == '+':
                    start = str(int(line.split(':')[1].split('-')[0])-200)
                    stop = str(int(line.split(':')[1].split('-')[0])+200)
                else:
                    start = str(int(line.split(':')[1].split('-')[1])-200)
                    stop = str(int(line.split(':')[1].split('-')[1])+200)

                bedfile.append([chrom,start,stop])

    return bedfile


def run_MEME(fastafile,outdir,scriptdir):
    os.system("qsub -v fastafile=" + fastafile + ",output=" + outdir + " " + scriptdir + "MEME_run.sh")

def run(bed,hg19fasta,outdir,filedir,scriptdir):
    a = BedTool(bed)
    a_fasta = bed + '.fasta'
    hg19 = BedTool(hg19fasta)


    a.sequence(fi=hg19fasta).save_seqs(a_fasta)
    run_MEME(a_fasta,outdir,scriptdir)

if __name__ == "__main__":
    #Home directory
    homedir = os.path.dirname(os.path.realpath(__file__))

    #File directory
    filedir = parent_dir(homedir) + '/files/'
    figuredir = parent_dir(homedir) + '/figures/'
    scriptdir = parent_dir(homedir) + '/scripts/'
    outdir = parent_dir(homedir) + '/MEME/'
    hg19fasta = '/scratch/Users/joru1876/hg19_reference_files/hg19_all.fa'

    bed = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/DE-Seq/A2N_ACN.genes.bed.count.bed.A2NACNnascent.resSig_pvalue.txt'

    bedfile = convert_deseqgenes_to_tssbed(bed)

    run(bedfile,hg19fasta,outdir,filedir,scriptdir)
