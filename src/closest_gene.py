__author__ = 'Jonathan Rubin'

from pybedtools import BedTool
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
import reflect_coverage
from operator import itemgetter
from scipy.stats import gaussian_kde
from scipy import stats
import numpy as np

def run(chipdir,refseq,filedir,DMSO,CA):
	TSS = (-200,1000)

	a = BedTool(chipdir)
	b = a.closest(refseq, d=True)
	b.cut([9,10,11,12,13,14,21]).saveas(filedir + '/SRF_closest.bed')
	d = dict()
	with open(filedir + '/SRF_closest.bed') as F:
		for line in F:
			line = line.strip().split()
			chrom,start,stop = line[0:3]
			d[chrom + '\t' + start + '\t' + stop + '\t'] = '\t'.join(line[3:])
	outfile = open(filedir + '/SRF_closest.rmdup.bed','w')
	for key in d:
		if '.' not in key.split():
			outfile.write(key+d[key]+'\n')
	outfile.close()
	# os.system("sort -k1,1 -k2,2n " + filedir + "/SRF_closest.rmdup.bed > " + filedir + "/SRF_closest.rmdup.sorted.bed")
	a = BedTool(filedir + '/SRF_closest.rmdup.bed')
	a.sort().saveas(filedir + '/SRF_closest.rmdup.sorted.bed')
	outfile = open(filedir + '/SRF.TSS.bed','w')
	outfile2 = open(filedir + '/SRF.gene.bed','w')
	with open(filedir + '/SRF_closest.rmdup.sorted.bed') as F:
		for line in F:
			chrom,start,stop,gene,number,strand,distance = line.strip().split()
			if int(stop) - int(start) > 2000 and int(distance) > 10000:
				if strand is '+':
					outfile.write(chrom+'\t'+str(int(start)+TSS[0])+'\t'+str(int(start)+TSS[1])+'\n')
					outfile2.write(chrom+'\t'+str(int(start)+TSS[1])+'\t'+stop+'\n')
				else:
					outfile.write(chrom+'\t'+str(int(stop)-TSS[1])+'\t'+str(int(stop)-TSS[0])+'\n')
					outfile2.write(chrom+'\t'+start+'\t'+str(int(stop)-TSS[1])+'\n')
	outfile.close()
	outfile2.close()
	a = BedTool(filedir + '/SRF.TSS.bed')
	a.sort().saveas(filedir + '/SRF.TSS.bed')
	a = BedTool(filedir + '/SRF.gene.bed')
	a.sort().saveas(filedir + '/SRF.gene.bed')

	TSS = filedir + '/SRF.TSS.bed'
	genes = filedir + '/SRF.gene.bed'

	os.system("bedtools map -a " + genes + " -b " + DMSO + " -c 4 -o sum > " + filedir + "/DMSO.genes.bed")
	os.system("bedtools map -a " + TSS + " -b " + DMSO + " -c 4 -o sum > " + filedir + "/DMSO.TSS.bed")
	os.system("bedtools map -a " + genes + " -b " + CA + " -c 4 -o sum > " + filedir + "/CA.genes.bed")
	os.system("bedtools map -a " + TSS + " -b " + CA + " -c 4 -o sum > " + filedir + "/CA.TSS.bed")

	TRx = list()
	TRy = list()
	expressionlist = list()

	with open(filedir + "/DMSO.genes.bed") as a, open(filedir + "/DMSO.TSS.bed") as b, open(filedir + "/CA.genes.bed") as c, open(filedir + "/CA.TSS.bed") as d:
		for line in a:
			bline = b.readline().strip().split()[-1]
			cline = c.readline().strip().split()[-1]
			dline = d.readline().strip().split()[-1]
			if line.strip().split()[-1] is '.':
				DMSOgene = 0.0
			else:
				DMSOgene = float(line.strip().split()[-1])
			if bline is '.':
				DMSOTSS = 0.0
			else:
				DMSOTSS = float(bline)
			if cline is '.':
				CAgene = 0.0
			else:
				CAgene = float(cline)
			if dline is '.':
				CATSS = 0.0
			else:
				CATSS = float(dline)
			if DMSOgene == 0.0:
				TRx.append(0.0)
			else:
				TRx.append((DMSOTSS/DMSOgene))
			if CAgene == 0.0:
				TRy.append(0.0)
			else:
				TRy.append((CATSS/CAgene))
			expressionlist.append((np.log2(DMSOgene)+np.log2(CAgene))/2.0)

	F6 = plt.figure()
	ax1 = F6.add_subplot(111)
	xy = np.vstack([TRx,TRy])
	z = gaussian_kde(xy)(xy)
	ax1.scatter(TRx,TRy,c=z,edgecolor="")
	# ax1.scatter(TRx2,TRy2,c='red',edgecolor="",s=expressionlist2)
	ax1.set_title('Pausing Index')
	ax1.set_ylabel('CA')
	ax1.set_xlabel('DMSO')
	ax1.get_xaxis().tick_bottom()
	ax1.get_yaxis().tick_left()
	#ax1.plot([0,1/slope1],[intercept1,1],color = 'r')
	ax1.set_xlim([0, 20])
	ax1.set_ylim([0, 20])
	ax1.plot([0,50.0],[0,50.0],color='k')
	# ax1.text(8,18, "Pearson = " + str(pearsons)[0:5])
	# ax2 = F6.add_subplot(122)
	# ax2.plot(np.sort(cdf),np.linspace(0,1,len(cdf)))
	# ax2.plot(stats.norm.cdf(np.linspace(min(cdf),max(cdf)),0,np.var(cdf)),np.linspace(0,1,len(cdf)))
	plt.savefig(figuredir + '/PausingIndex.png')

			


	


if __name__ == "__main__":
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

	chipdir = '/scratch/Shares/dowell/ENCODE/old/HCT116/SRF/peak_files/ENCFF001UEM.bed'
	refseq = '/scratch/Users/joru1876/hg19_reference_files/refFlat_hg19.bed'
	#Specify DMSO treated bedgraph directory
	DMSO = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'

	#Specify CA treated bedgraph directory
	CA = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/genomecoveragebed/fortdf/JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.BedGraph.mp.BedGraph'


	reflect_coverage.run(DMSO,CA,filedir)
	DMSOreflect = filedir + '/DMSO.bedgraph'
	CAreflect = filedir + '/CA.bedgraph'
	run(chipdir,refseq,filedir,DMSOreflect,CAreflect)