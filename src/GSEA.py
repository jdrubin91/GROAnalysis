__author__ = 'Jonathan Rubin'

import os
import sys

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
    os.system("bedtools merge -i " + filedir + "all_preliminary_bidir.bed > all_preliminary_bidir.merge.bed")
    os.system("sort -k1,1 -k2,2n " + filedir + "all_preliminary_bidir.merge.bed > " + filedir + "all_preliminary_bidir.merge.sort.bed")
    os.system("bedtools multicov -bams " + condition1bam + " " + condition2bam + " -bed " + filedir + "all_preliminary_bidir.merge.sort.bed >" + filedir + "all_preliminary_bidir.merge.sort.count.bed")

    return filedir + "all_preliminary_bidir.merge.sort.count.bed"

def create_intersect_file(interval_file,path_to_PSSMs,filedir):
    header = list()
    i = 0
    for TF in os.listdir(path_to_PSSMs):
        header.append(TF)
        if i == 0:
            os.system("bedtools intersect -u -a " + interval_file + " -b " + path_to_PSSMs + TF + " > " + filedir + "all_preliminary_bidir.merge.sort.count.intersect.bed")
        else:
            os.system("bedtools intersect -u -a " + interval_file + " -b " + path_to_PSSMs + TF + " > " + filedir + "temp.bed")
            concatenate_files(filedir + "all_preliminary_bidir.merge.sort.count.intersect.bed", filedir + "temp.bed")
        i += 1



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
    filedir = "../files/"
    bamfolder = '/projects/dowellLab/Taatjes/170207_K00262_0069_AHHMHVBBXX/cat/trimmed/flipped/bowtie2/sortedbam/'
    condition1bam = 'J52_trimmed.flip.fastq.bowtie2.sorted.bam'
    condition2bam = 'J62_trimmed.flip.fastq.bowtie2.sorted.bam'



    interval_file = create_bidir_interval_file(filelist,filedir,condition1bam,condition2bam)
    create_intersect_file(interval_file,path_to_PSSMs,filedir)



