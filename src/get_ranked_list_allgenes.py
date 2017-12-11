import sys
import math

#Use DE-Seq output .res.txt
file1 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/hg19.genes.all_replicates.count.bed.DMSO45CA45nascent.res.txt'

outfile = open(file1+'.ranked_list.rnk','w')

names = list()
fcs = list()
with open(file1) as F:
	F.readline()
	for line in F:
		line = line.strip('\n').split('\t')
		name = line[1].split(';')[1]
		fc = line[5]
		if name not in names:
			names.append(name)
			if fc != 'NA' and fc != 'Inf':
				fcs.append(fc)
			else:
				fcs.append(0)

for i in range(len(names)):
	name = names[i]
	fc = str(fcs[i])
	outfile.write(name + '\t' + fc + '\n')
