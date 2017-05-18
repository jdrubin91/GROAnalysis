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
            line = line.strip().split()
            if "id" not in line[0] and "NA" not in line[0]:
                line = line[1].split(';')[-1]
                chrom = line.split(':')[0]
                if line.split('_')[-1] == '+':
                    start = str(int(line.split(':')[1].split('_')[0].split('-')[0])-200)
                    stop = str(int(line.split(':')[1].split('_')[0].split('-')[0])+200)
                else:
                    start = str(int(line.split(':')[1].split('_')[0].split('-')[1])-200)
                    stop = str(int(line.split(':')[1].split('_')[0].split('-')[1])+200)

                bedfile.append([chrom,start,stop])

    BedTool(bedfile).saveas(file1 + '.tss.bed')

    return file1 + '.tss.bed'

def convert_joeydeseq_to_bed(file1):
    bedfile = list()
    with open(file1) as F:
        for line in F:
            line = line.strip().split(',')
            if 'baseMean' not in line[0] and 'NA' not in line[-1]:
                if float(line[-1]) < 0.01:
                    chrom = line[0].split(':')[0]
                    start = line[0].split(':')[1].split('-')[0]
                    stop = line[0].split(':')[1].split('-')[1]
                    bedfile.append([chrom,start,stop])

    print bedfile[0:10]
    BedTool(bedfile).saveas(file1 + '.bed')

    return file1 + '.bed'


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

    # bed = '/projects/dowellLab/Taatjes/170413_K00262_0087_AHJLW5BBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/DE-Seq/A2N_ACN.genes.bed.count.bed.A2NACNnascent.resSig_pvalue.txt'
    # bedfile = convert_deseqgenes_to_tssbed(bed)

    joeyfile = '/scratch/Users/joru1876/ButcherDrugRes/DEseq/out/Enhancers.csv'
    bedfile = convert_joeydeseq_to_bed(joeyfile)

    run(bedfile,hg19fasta,outdir,filedir,scriptdir)
