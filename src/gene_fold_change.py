__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math

def setBoxColors(bp):    
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

def run(countsfile,savedir,genes):
    with open(countsfile) as F:
        samples = F.readline().strip('\n').split('\t')[6:]
        factornorm = [0]*len(samples)
        for line in F:
            line = line.strip('\n').split('\t')[6:]
            for i in range(len(line)):
                factornorm[i] += float(line[i])

    print factornorm
    # for i in range(len(factornorm)):
    #     factornorm[i] = factornorm[i]/factornorm[0]
    factornorm = [factornorm[0]/factornorm[i] for i in range(len(factornorm))]
    print factornorm

    with open(countsfile) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            gene = line[3].split(';')[1]
            if gene in genes:
                print gene
                vals = [a*b for a,b in zip([float(i) for i in line[6:]],factornorm)]
                fc = [i/((vals[2]+vals[3])/2) for i in vals]
                print fc
                # fc = [np.mean([fc[0],fc[1]]),np.mean([fc[2],fc[3]]),np.mean([fc[4],fc[5]]),np.mean([fc[6],fc[7]])]
                ind = np.arange(len(fc))
                F = plt.figure()
                ax = F.add_subplot(111)
                ax.bar(ind,fc,0.3)
                plt.savefig(savedir + gene + '_NutlinPetide_fc.svg')

def make_boxplots(countsfile,savedir,genelist):
    with open(countsfile) as F:
        samples = F.readline().strip('\n').split('\t')[6:]
        factornorm = [0]*len(samples)
        for line in F:
            line = line.strip('\n').split('\t')[6:]
            for i in range(len(line)):
                factornorm[i] += float(line[i])

    factornorm = [factornorm[0]/factornorm[i] for i in range(len(factornorm))]


    peptide = list()
    water = list()
    with open(countsfile) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            gene = line[3].split(';')[1]
            if gene in genelist:
                try:
                    vals = [a*b for a,b in zip([float(i) for i in line[6:]],factornorm)]
                    fc = [math.log(i/((vals[2]+vals[3])/2),2) for i in vals]
                    water.append(fc[-1])
                    water.append(fc[-2])
                    peptide.append(fc[-3])
                    peptide.append(fc[-4])
                except:
                    pass
                # fc = [np.mean([fc[0],fc[1]]),np.mean([fc[2],fc[3]]),np.mean([fc[4],fc[5]]),np.mean([fc[6],fc[7]])]


    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('p53 Gene Targets Fold Change')
    ax.set_ylabel('Log2 Fold-Change (Nutlin/DMSO)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.axhline(0, color='black',alpha=0.5)
    ax.set_axisbelow(True)
    bp = ax.boxplot(water, positions = [1], patch_artist=True)
    setBoxColors(bp)
    bp2 = ax.boxplot(peptide, positions = [3], patch_artist=True)
    setBoxColors(bp2)
    ax.set_xlim([0, 4])
    plt.xticks([1,3], ['-Peptide','+Peptide'])

    F.savefig(savedir + 'p53_genetargets_fc_Nutlin_peptide.png', dpi=1200)
    F.savefig(savedir + 'p53_genetargets_fc_Nutlin_peptide.svg')

def make_MA_plot(countsfile,savedir,name_order,condition1,condition2):
    with open(countsfile) as F:
        samples = F.readline().strip('\n').split('\t')[6:]
        factornorm = [0]*len(samples)
        for line in F:
            line = line.strip('\n').split('\t')[6:]
            for i in range(len(line)):
                factornorm[i] += float(line[i])

    factornorm = [factornorm[0]/factornorm[i] for i in range(len(factornorm))]


    x = list()
    y = list()
    index1 = name_order.index(condition1)
    index2 = name_order.index(condition2)
    with open(countsfile) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            gene = line[3].split(';')[1]
            vals = [a*b for a,b in zip([float(i) for i in line[6:]],factornorm)]
            val1 = vals[index1]
            val2 = vals[index2]
            try:
                x.append(math.log((val1+val2)/2,10))
            except:
                x.append(0)
            try:
                y.append(math.log(val1/val2,2))
            except:
                y.append(0)

    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(x,y,edgecolor="",s=14,alpha=0.5) 
    ax.set_title('Gene Expression MA Plot ' + condition1 + '-' + condition2)
    ax.set_ylabel('Log2 Fold-Change (' + condition1 + '/' + condition2 + ')')
    ax.set_xlabel('Log10 Average Expression')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # plt.show()
    plt.savefig(savedir + 'Gene_expression_MA_plot_'+condition1+ '_' + condition2+'.png',dpi=1200)
    # plt.savefig(figuredir + 'ExonExpression_GCcontent_RNASeq_normalized_poly2.svg')



if __name__ == "__main__":
    countsfile = '../files/all_samples.genes.bed.count.bed'
    savedir = '../figures/'
    # genes = ['BAX','BTG2','CDKN1A','MDM2']
    # run(countsfile,savedir,genes)

    genelistrep1 = ['PGAP1', 'SYTL1', 'CPM', 'ZNF337', 'PLK3', 'ZFPM1', 'ZNF79', 'FDXR', 'TM7SF3', 'BLOC1S2', 'SERPINB5', 'WDR63', 'HHAT', 'XPC', 'PAPPA', 'REV3L', 'OSBPL3', 'MDM2', 'TRAF3IP2-AS1', 'RNF19B', 'NTPCR', 'CCNG1', 'RPS27L', 'IKBIP', 'CDKN1A', 'C12orf5', 'PHLDA3', 'AEN', 'PANK1', 'PPM1D', 'TRIAP1', 'RRM2B', 'KSR1', 'GDF15', 'CCDC90B', 'PRKAB1', 'FAM212B', 'TP53I3', 'ASTN2', 'TP53INP1', 'TAF3', 'GADD45A', 'ASCC3', 'CEP85L', 'SESN1', 'SESN2', 'TNFRSF10C', 'TNFRSF10B', 'SERHL2', 'TNFRSF10D', 'PSTPIP2', 'POLH', 'DRAM1', 'BTG2', 'CD82', 'FHL2', 'FBXO22', 'DDB2', 'APAF1', 'GPR87', 'E2F7', 'PLCL2', 'PTCHD4', 'CMBL', 'NEAT1', 'SMEK1', 'RPPH1', 'CYFIP2', 'SLC44A5', 'NADSYN1', 'BAX', 'ZMAT3', 'BBC3', 'PRKX', 'ACER2', 'FAS', 'SULF2', 'DGKA', 'ANXA4', 'KIAA0247', 'TANC1', 'RHOBTB2', 'BBS9', 'PLEKHG1', 'LOC100128505', 'GRIN2C', 'PTPRE', 'PTPRD', 'ATF3']
    # make_boxplots(countsfile,savedir,genelistrep1)

    name_order = ['DP1','DP2','DW1','DW2','NP1','NP2','NW1','NW2']
    condition1 = 'NW1'
    condition2 = 'NW2'
    # make_MA_plot(countsfile,savedir,name_order,condition1,condition2)

    genelistrep1and2 = ['CD82', 'SYTL1', 'CYFIP2', 'CDKN1A', 'ZFPM1', 'NADSYN1', 'TNFRSF10C', 'ZNF79', 'BTBD19', 'FDXR', 'MFNG', 'BBC3', 'PHLDA3', 'POLH', 'ACER2', 'PLXNB3', 'ADAMTS7', 'PPM1D', 'BTG2', 'FAS', 'SLC22A20', 'SULF2', 'HFM1', 'SLCO5A1', 'CDC42BPG', 'DOK7', 'DDB2', 'DGKA', 'FAS-AS1', 'PTCHD4', 'FRMPD2', 'COL7A1', 'LOC100128505', 'GRIN2C', 'ZNF506']
    print len(genelistrep1and2)
    normalized_countsfile = '../files/all_samples.genes.count.gc_normalized_NW2.bed'
    # make_boxplots(countsfile,savedir,genelistrep1and2)

    allen2014genelist = ['ABCA1','ABCA12','ABCB9','ACER2','ACTA2','ADAMTS7','ADIG','AEN','ALOX5','AMZ2','AMZ2P1','ANKRA2','ANKUB1','ANXA4','APAF1','APOBEC3C','APOBEC3D','APOBEC3H','ARHGAP30','ASCC3','ASS1','ASTN2','ATP2B2','BAX','BBC3','BLOC1S2','BTG2','C12orf5','C17orf82','C2orf71','C3orf20','C6orf58','CBS','CCBP2','CCDC90B','CCNG1','CD82','CDC42BPG','CDKN1A','CEP85L','CFL1P1','CHAC1','CMBL','COL17A1','COL20A1','COL4A1','COL5A1','COL6A3','CPE','CPN1','CPO','CYFIP2','DAO','DCP1B','DDB2','DDIT4','DGKA','DGKI','DOCK8','DPEP1','DRAM1','DUOX1','EBI3','EFNB1','EPN3','EPS8L2','ESRRB','FAM183A','FAM210B','FAM212B','FAS','FBXO22','FBXW7','FDXR','FLJ43681','FLJ44511','FRMPD2','GBE1','GDF15','GJB5','GPR56','GPR87','GRIN2C','HES2','HHAT','HSPA4L','ICAM1','IKBIP','INPP5D','ISCU','ITGA3','ITGA9','KANK3','KCNC2','KCNN4','KDM4B','KIAA0247','KLHDC7A','KRT18P55','KSR1','LAMA3','LAPTM5','LDLRAD1','LINC00158','LINC00589','LINC00663','LOC100294145','LOC100506343','LOC283761','LOC284080','LOC284385','LOC643401','LRP1','LYNX1','MDM2','MFSD4','MGAM','MICALL1','MIR1204','MIR3189' ,'MIR34A','MIR4679-1','MIR4692','MYBPC3','NADSYN1','NINJ1','NTF3','NTPCR','NTRK2','ORAI3','PANK1','PARD6G','PGAP1','PGPEP1','PHLDA3','PKP1','PLCL2','PLEKHG1','PLK2','PLXNB3','PML','PNLIPRP2','POLH','PPM1D','PRDM1','PRKAB1','PRKX','PROM2','PSTPIP2','PTCHD4','PTP4A1','PTPRE','PTPRU','PVRL4','PVT1','RFNG','RHOD','RINL','RNF157-AS1','RPS27L','RRM2B','SAC3D1','SCN2A','SDC1','SERPINB5','SESN1','SESN2','SF3B14','SLC30A1','SLC44A5','SMAD5-AS1','SULF2','SYTL1','TAF3','TAKR','TCERG1L','TCP11L1','TM7SF3','TMEM63B','TNFRSF10B','TNFRSF10C','TOB1','TP53I11','TP53I3','TP53INP1','TRAF3IP2-AS1','TRAF4','TRIAP1','TSKU','UNC5B','VCAN','WDR63','WRAP73','XPC','ZMAT3','ZNF219','ZNF337','ZNF79']
    print len(allen2014genelist)
    # make_boxplots(normalized_countsfile,savedir,allen2014genelist)



