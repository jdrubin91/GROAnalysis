__author__ = 'Jonathan Rubin'

import os
import sys
import pybedtools

def run(folder,genes):
    for file1 in os.listdir(folder):
        if '.mp.BedGraph' in file1 and 'tdf' not in file1:
            a = pybedtools.BedTool(genes).sort()
            b = pybedtools.BedTool(folder + file1).sort()
            a.map(b,0="sum").saveas(folder + file1 + '.gene_counts.bed')

    


if __name__ == "__main__":
    folder = sys.argv[1]
    genes = ''
    run(folder,genes)

