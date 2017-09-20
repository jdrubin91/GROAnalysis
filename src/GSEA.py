__author__ = 'Jonathan Rubin'

import os
import sys
import math

def add_header(header,file1):
    lines = list()
    with open(file1) as F:
        for line in F:
            lines.append(line)

    outfile = open(file1,'w')
    outfile.write(header)
    for line in lines:
        outfile.write(line)

def concatenate_files(file1,file2):
    all_lines = list()
    with open(file1) as F:
        with open(file2) as F2:
            for line in F:
                line2 = F2.readline()
                line = line.strip('\n').split('\t')
                line2 = line2.strip('\n').split('\t')
                line.append(line2[-1])
                all_lines.append(line)

    outfile = open(file1,'w')
    for line in all_lines:
        outfile.write('\t'.join(line) + '\n')

def create_bidir_interval_file(filelist,filedir,condition1bam,condition2bam):
    os.system("cat " + " ".join(filelist) + " > " + filedir + "all_preliminary_bidir.bed")
    os.system("sort -k1,1 -k2,2n " + filedir + "all_preliminary_bidir.bed > " + filedir + "all_preliminary_bidir.sort.bed")
    os.system("bedtools merge -i " + filedir + "all_preliminary_bidir.sort.bed > " + filedir + "all_preliminary_bidir.merge.bed")
    os.system("sort -k1,1 -k2,2n " + filedir + "all_preliminary_bidir.merge.bed > " + filedir + "all_preliminary_bidir.merge.sort.bed")
    os.system("bedtools multicov -bams " + condition1bam + " " + condition2bam + " -bed " + filedir + "all_preliminary_bidir.merge.sort.bed >" + filedir + "all_preliminary_bidir.merge.sort.count.bed")

def create_intersect_file(interval_file,path_to_PSSMs,filedir):
    header = ['chrom','start','stop','cond1counts','cond2counts']
    i = 0
    for TF in os.listdir(path_to_PSSMs):
        header.append(TF)
        if i == 0:
            os.system("bedtools intersect -c -a " + interval_file + " -b " + path_to_PSSMs + TF + " > " + filedir + "all_preliminary_bidir.merge.sort.count.intersect.bed")
        else:
            os.system("bedtools intersect -c -a " + interval_file + " -b " + path_to_PSSMs + TF + " > " + filedir + "temp.bed")
            concatenate_files(filedir + "all_preliminary_bidir.merge.sort.count.intersect.bed", filedir + "temp.bed")
        i += 1

    add_header('\t'.join(header)+'\n', filedir + "all_preliminary_bidir.merge.sort.count.intersect.bed")

def sort_intersect_file(intersect_file,filedir):
    fcs = list()
    intervals = list()
    with open(intersect_file) as F:
        TFnames = F.readline().strip('\n').split('\t')[5:]
        TFvals = [[] for i in range(len(TFnames))]
        for line in F:
            line = line.strip('\n').split('\t')
            intervals.append('\t'.join(line[:3]))
            try:
                log2fc = math.log(float(line[3])/float(line[4]),2)
                fcs.append(log2fc)
            except:
                fcs.append(0)
            for i in range(len(line[5:])):
                val = line[5+i]
                if val != 0:
                    TFvals[i].append(1)
                else:
                    TFvals[i].append(0)

    indices = [i[0] for i in sorted(enumerate(fcs), key=lambda x:x[1],reverse=True)]

    outfile = open(filedir + "all_preliminary_bidir.merge.sort.count.intersect.sorted.bed",'w')
    outfile.write('chrom\tstart\tstop\tlog2foldchange\t' + '\t'.join(TFnames) + '\n')
    for i in indices:
        outfile.write(intervals[i] + '\t' + str(fcs[i]))
        for vals in TFvals:
            outfile.write('\t' + vals[i])
        outfile.write('\n')






if __name__ == "__main__":
    folder1 = '/projects/dowellLab/Taatjes/170825_NB501447_0152_fastq/Demux/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit/'
    folder2 = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/Tfit_run2/'
    filelist = [folder1+'foot_print_testing-1_prelim_bidir_hits.bed',
                folder1+'foot_print_testing-2_prelim_bidir_hits.bed',
                folder1+'foot_print_testing-3_prelim_bidir_hits.bed',
                folder1+'foot_print_testing-4_prelim_bidir_hits.bed',
                folder1+'foot_print_testing-5_prelim_bidir_hits.bed',
                folder1+'foot_print_testing-6_prelim_bidir_hits.bed',
                folder2+'foot_print_testing-1_prelim_bidir_hits.bed',
                folder2+'foot_print_testing-2_prelim_bidir_hits.bed',
                folder2+'foot_print_testing-3_prelim_bidir_hits.bed',
                folder2+'foot_print_testing-4_prelim_bidir_hits.bed',
                folder2+'foot_print_testing-5_prelim_bidir_hits.bed',
                folder2+'foot_print_testing-6_prelim_bidir_hits.bed']
    filedir = "/Users/joru1876/scratch_backup/GROAnalysis/files/"
    figuredir = "/Users/joru1876/scratch_backup/GROAnalysis/figures/"
    bamfolder = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/'
    condition1bam = bamfolder + 'J52_trimmed.flip.fastq.bowtie2.sorted.bam'
    condition2bam = bamfolder + 'J62_trimmed.flip.fastq.bowtie2.sorted.bam'
    path_to_PSSMs = '/scratch/Shares/dowell/md_score_paper/PSSM_hits_genome_wide/pval-7/'

    # create_bidir_interval_file(filelist,filedir,condition1bam,condition2bam)
    # interval_file = filedir + "all_preliminary_bidir.merge.sort.count.bed"
    # create_intersect_file(interval_file,path_to_PSSMs,filedir)



    intersect_file = filedir + "all_preliminary_bidir.merge.sort.count.intersect.bed"
    sort_intersect_file(intersect_file,filedir)


