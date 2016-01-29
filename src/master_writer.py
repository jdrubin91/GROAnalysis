__author__ = 'Jonathan Rubin'

def run(DMSOgenes,DMSOTSS,DMSOEND,CAgenes,CATSS,CAEND,filedir):    
    d = dict()
    with open(DMSOgenes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            d[gene] = [chrom,start,stop,number,strand,coverage]
    with open(DMSOTSS) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            d[gene].append(coverage)
            
    with open(DMSOEND) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            d[gene].append(coverage)
            
    with open(CAgenes) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            d[gene].append(coverage)
            
    with open(CATSS) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            d[gene].append(coverage)
            
    with open(CAEND) as F1:
        for line in F1:
            chrom,start,stop,gene,number,strand,coverage = line.strip().split()
            d[gene].append(coverage)
    
    outfile = open(filedir + '/Master.bed','w')
    outfile.write('Gene\tChrom\tStart\tStop\tNumber\tStrand\tDMSO gene body\tDMSO TSS\tDMSO END\tCA gene body\tCA TSS\tCA END\n')
    for gene in d:
        outfile.write(gene + '\t')
        for item in d[gene]:
            outfile.write(item + '\t')
        outfile.write('\n')
    
    