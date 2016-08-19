__author__ = 'Jonathan Rubin'

	def run(DMSO,CA,filedir):
	list1 = list()
	with open(DMSO) as F1:
		for line in F1:
			chrom, start, stop, coverage = line.strip().split()
			if float(coverage) < 0:
				coverage = -float(coverage)
			else:
				coverage = float(coverage)
			list1.append((chrom,start,stop,coverage))
	outfile1 = open(filedir + '/DMSO.bedgraph', 'w')
	for item in list1:
		outfile1.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + str(item[3]) + '\n')
	list1 = list()
	with open(CA) as F2:
		for line in F2:
			chrom, start, stop, coverage = line.strip().split()
			if float(coverage) < 0:
				coverage = -float(coverage)
			else:
				coverage = float(coverage)
			list1.append((chrom,start,stop,coverage))
	outfile2 = open(filedir + '/CA.bedgraph', 'w')
	for item in list1:
		outfile2.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + str(item[3]) + '\n')