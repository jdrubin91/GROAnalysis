__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
import math

def run(file1,figuredir):
    fcs = list()
    PEPs = list()
    Intensities = list()
    genes = list()
    rankedlist = list()
    with open(file1) as F:
        F.readline()
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            fc = math.log(float(line[1]),2)
            try:
                PEP = -math.log(float(line[9]),10)
            except ValueError:
                PEP = 1
            Intensity = math.log(float(line[15]),10)
            gene = line[23].split(';')[0]
            fcs.append(fc)
            PEPs.append(PEP)
            Intensities.append(Intensity)
            genes.append(gene)
            if math.isnan(float(line[1])):
                realfc = 1
            else:
                realfc = float(line[1])
            rankedlist.append((1/realfc,gene))
            if PEP > 100:
                print gene

    outfile = open('/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/misc/HCT116_Serinduction_phospho_rep1.fc_sorted_reverse_test.rnk','w')
    print len(rankedlist)
    for tup in sorted(rankedlist, key=lambda x: x[0],reverse=True):
        if not len(tup[1]) < 1: 
            outfile.write(tup[1] + '\t' + str(tup[0]) + '\n')

    # F = plt.figure()
    # ax = F.add_subplot(111)
    # plt.scatter(Intensities,fcs,c=PEPs,edgecolor="")
    # plt.show()

if __name__ == "__main__":
    file1 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/misc/HCT116_Serinduction_phospho_rep1.txt'
    figuredir = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/GROAnalysis/figures/'
    run(file1,figuredir)
