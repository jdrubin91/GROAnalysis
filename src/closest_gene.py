__author__ = 'Jonathan Rubin'

from pybedtools import BedTool

def run(chipdir,refseq,filedir):
	a = BedTool(chipdir)
	b = a.closest(refseq)
	b.saveas(filedir + '/SRF_closest.bed')

	return filedir + '/SRF_closest.bed'
			


	


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
	refseq = '/scratch/Users/joru1876/refFlat_hg19.bed'
	run(chipdir,refseq,filedir)