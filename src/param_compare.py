__author__ = "Jonathan Rubin"

import intervals

file1 = '/scratch/Shares/dowell/ENCODE/Rubin2015_DMSO-2_divergent_classifications.bed'
file2 = '/scratch/Shares/dowell/ENCODE/Rubin2015_DMSO-1_divergent_classifications.bed'

#Runs interval search over all bed sites in both files, recovers parameters in 
#bidirectional model
def run(file1,file2):
    param = 0
    x = list()
    y = list()
    d1 = dict()
    d2 = dict()
    A = list()
    B = list()
    
    with open(file1) as F1:
        for line in F1:
            if not '#' in line[0]:
                chrom, start, stop, param = line.strip().split()
                param = param.split('|')[1].split(',')
                d1[chrom + ':' + start + '-' + stop] = param
                A.append((int(start),int(stop),'A',chrom))
    with open(file2) as F2:
        for line in F2:
            if not '#' in line[0]:
                chrom, start, stop, param = line.strip().split()
                param = param.split('|')[1].split(',')
                d2[chrom + ':' + start + '-' + stop] = param
                B.append((int(start),int(stop),'B',chrom))
                
    ST = intervals.comparison((A,B))
    OVERLAPS_0_1 = ST.find_overlaps(0,1)
    print "Overlap Instances: " + str(len(OVERLAPS_0_1))
    for O in OVERLAPS_0_1:
        if not len(O.overlaps.keys()) > 2:
            for interval_original in O.overlaps:
   	        if 'A' in interval_original.INFO:
              	    x.append(d1[interval_original.INFO[1] + ':' + interval_original.start + '-' + interval_original.stop][param])
              	elif 'B' in interval_original.INFO:
              	    y.append(d2[interval_original.INFO[1] + ':' + interval_original.start + '-' + interval_original.stop][param])
    
    return x,y
#Creates a scatter plot of x and y values
#def plot(x,y):
#    return

if __name__ == "__main__":
    x,y = run(file1,file2)
    print "x: " + str(len(x)) + "\ny: " + str(len(y))
    #plot(x,y)