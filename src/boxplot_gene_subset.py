__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import kstest
import numpy as np

def setBoxColors(bp):    
    ## change outline color, fill color and linewidth of the boxes
    # change outline color
    for box in bp['boxes']:
        box.set( color='#7570b3', linewidth=2)
    # change fill color
    bp['boxes'][0].set( facecolor = 'green' )
    bp['boxes'][1].set( facecolor = 'blue' )
    bp['boxes'][2].set( facecolor = 'magenta' )
    bp['boxes'][3].set( facecolor = 'red' )

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


def run(g15up, g15do, g45up, g45do):
    pvalcut = 0.01
    timepoints = ['0','15','45']
    i = 0
    d = dict()
    for file1 in folder:
        time = timepoints[i]
        i+=1
        with open(file1) as F:
            F.readline()
            for line in F:
                line = line.strip().split()
                gene = line[1]
                try:
                    pval = float(line[-2])
                except ValueError:
                    pval = 1
                try:
                    log2change = float(line[-3])
                except ValueError:
                    log2change = 0
                if gene not in d:
                    d[gene] = list()
                d[gene].append((time,log2change,pval))


    t0 = [[] for i in range(4)]
    for key in d:
        for time,log2change,pval in d[key]:
            if time == '0':
                gene = key.split(';')[1]
                if gene in g15up:
                    t0[0].append(log2change)
                if gene in g45up:
                    t0[1].append(log2change)
                if gene in g45do:
                    t0[2].append(log2change)
                if gene in g45up:
                    t0[3].append(log2change)

    t15 = [[] for i in range(4)]
    for key in d:
        for time,log2change,pval in d[key]:
            if time == '15':
                gene = key.split(';')[1]
                if gene in g15up:
                    t15[0].append(log2change)
                if gene in g45up:
                    t15[1].append(log2change)
                if gene in g45do:
                    t15[2].append(log2change)
                if gene in g45up:
                    t15[3].append(log2change)


    t45 = [[] for i in range(4)]
    for key in d:
        for time,log2change,pval in d[key]:
            if time == '45':
                gene = key.split(';')[1]
                if gene in g15up:
                    t45[0].append(log2change)
                if gene in g45up:
                    t45[1].append(log2change)
                if gene in g45do:
                    t45[2].append(log2change)
                if gene in g45up:
                    t45[3].append(log2change)


    F = plt.figure()
    ax = F.add_subplot(111)
    ax.set_title('CA-dependent Gene Fold-Changes')
    ax.set_ylabel('Log$_2$ Fold Change (CA/DMSO)')
    ax.set_xlabel('Time after Serum Induction (min)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    # plt.axhline(0, color='black')
    bp = ax.boxplot(t0, positions = [1,2,3,4], patch_artist=True,whis=5)
    setBoxColors(bp)
    bp2 = ax.boxplot(t15, positions = [6,7,8,9], patch_artist=True,whis=5)
    setBoxColors(bp2)
    bp3 = ax.boxplot(t45, positions = [11,12,13,14], patch_artist=True,whis=5)
    setBoxColors(bp3)
    ax.set_xlim([0, 15])
    plt.xticks([2.5,7.5,12.5], timepoints)
    plt.axhline(0, color='black',alpha=0.3)
    ax.set_axisbelow(True)

    green_patch = mpatches.Patch(color='green', label='UP at 15min')
    red_patch = mpatches.Patch(color='red', label='DOWN at 15min')
    blue_patch = mpatches.Patch(color='blue', label='UP at 45min')
    magenta_patch = mpatches.Patch(color='magenta', label='DOWN at 45min')
    ax.legend([green_patch,blue_patch,magenta_patch,red_patch],['Genes up at 15min','Genes up at 45min','Genes down at 45min','Genes down at 15min'],bbox_to_anchor=(0, 0),loc=3,fontsize=10)
    F.savefig(savedir + 'CA_Gene_Changes.png', dpi=1200)



if __name__ == "__main__":
    g15up = ['MYH9', 'LIF', 'NFKB2', 'KRT80', 'LAMB3', 'SORBS2', 'ARPC4', 'MIR4530', 'GUK1', 'CBX4', 'TPM4', 'ZFP36', 'WEE1', 'BTG2', 'RSU1', 'INSIG1', 'CGA', 'MAFF', 'LOC100506810', 'YWHAZ', 'DUSP5', 'LIMA1', 'CREM', 'ODC1', 'GPR3', 'KCNK1', 'CASP9', 'ECE1', 'MLF1', 'EZR', 'ZNF548', 'LMNA', 'GDF15', 'AREG', 'TRIM39', 'PDE10A', 'FAM122C', 'SKIL', 'UBC', 'MYPN', 'PER1', 'MIR22HG', 'CYP24A1', 'SLC25A25', 'KLF6', 'PLAUR', 'LONRF2', 'MEG9', 'GLA', 'EGR3', 'MYADM', 'ETV3', 'DCUN1D3', 'MNT', 'MYL9', 'GNAT3', 'UAP1', 'YPEL5', 'SLC2A3', 'FOSB', 'CYR61', 'SNX32', 'ATF3', 'TGFBI', 'LATS2', 'MCL1', 'RIC3', 'DDX5', 'LOC100129427', 'ZNF323', 'CRIP2', 'KIFC1', 'ALPK2', 'TXN2', 'FOSL1', 'VCL', 'GPRC5A', 'MTCH1', 'CAV1', 'LINC00346', 'ANXA3', 'GNAL', 'TRIP13', 'TEX26-AS1', 'HOMER1', 'RPS16P5', 'NR4A1', 'PNPLA8', 'KIAA1683', 'C6orf99', 'MYL6', 'CHMP1B', 'CKS2', 'HSP90B1', 'NLRP2', 'TMEM66', 'PLEKHG2', 'SERTAD3', 'ARL5B', 'TULP2', 'PDLIM7', 'PTGES3L', 'LINC00475', 'CCDC80', 'PDLIM3', 'MYO1E', 'ARHGDIB', 'ELMSAN1', 'IBA57', 'COL1A1', 'ACTB', 'ELL2', 'DCLK1', 'LOC100507066', 'VASP', 'C1orf110', 'PPP1R11', 'SYAP1', 'TUFT1', 'LINC00152', 'ACTG1', 'MTHFD1L', 'FHL2', 'KIF9', 'ITGB1', 'MAP2K3', 'PVRL2', 'RARRES1', 'CYCS', 'CBX8', 'ATP6V0A1', 'GEM', 'GRK6', 'MAPKAPK2', 'C17orf110', 'SYTL3', 'GP6', 'DNAJB2', 'ME3', 'CCDC173', 'TRIB1', 'ARHGAP23', 'DDC', 'MTFP1', 'TNFRSF10B', 'ARF4', 'USP11', 'IFRD1', 'NFKBIZ', 'GPRC5C', 'PRKAR1A', 'TSC22D3', 'JOSD1', 'EREG', 'SGIP1', 'CCNB1IP1', 'DNAJB4', 'ETF1', 'SIK1', 'NEAT1', 'KCTD10', 'TMEM40', 'LPP-AS2', 'PDLIM5', 'ATP1B1', 'CNN2', 'MAP3K8', 'PLK3', 'ANAPC15', 'SERPINE1', 'CYSTM1', 'LDB3', 'MUS81', 'ELOVL5', 'USP36', 'LDLR', 'IDI1', 'ASTL', 'PER2', 'RTN2', 'COQ7', 'CFL1', 'DAPK3', 'MEF2D', 'CSRNP1', 'RAB11B', 'CD55', 'NR4A3', 'GBP2', 'IER2', 'EGR2', 'NR1H2', 'AVPI1', 'DUSP4', 'SLC20A1', 'C14orf80', 'SELK', 'GPR112', 'LOC541471', 'LOC100505782', 'LOC201651', 'SORCS3', 'UCA1', 'REL', 'EMD', 'KCNK3', 'ATP1B3', 'SEPT10', 'ACTN4', 'TTC19', 'LOC100131434', 'DUSP2', 'FLNA', 'JUN', 'LOC284454', 'ARID5A', 'TRAPPC2P1', 'MIR3132', 'IL8', 'DDAH1', 'CCDC9', 'ZYX', 'GDA', 'VPS37B', 'LOC100128675', 'THBS1', 'C12orf44', 'PLOD1', 'MYL12A', 'NAMPT', 'SCARA5', 'ERRFI1', 'BHLHE40', 'LINC00589', 'ARSG', 'ISG20L2', 'DUSP3', 'TPM2', 'FAM100A', 'SPTSSA', 'ANXA1', 'RGCC', 'ANKRD37', 'CAPN13', 'HERPUD1', 'FLOT1', 'PPP1R15A', 'ALG13', 'MAFK', 'RAB11B-AS1', 'MIR23A', 'STRAP', 'MIR554', 'DUSP1', 'LOC644936', 'TSPYL2', 'VMP1', 'LIX1L', 'RELB', 'EGR1', 'KLF4', 'BCL10', 'MFSD12', 'SRF', 'PTP4A1', 'C10orf107', 'LOC731223', 'LOC338758', 'BCAS1', 'DUSP8', 'LOC100216001', 'EDN1', 'CDCP1', 'TBX18', 'IRF2BP2', 'JUNB', 'SGMS2', 'QPCTL', 'FBXO46', 'WDR1', 'GABARAPL1', 'DSTN', 'MAP3K13', 'DNAJB11', 'RCAN1', 'RBMXL1', 'RRP12', 'ISG15', 'NR4A2', 'MAP3K14', 'TCTEX1D4', 'FOSL2', 'FGFBP1', 'MIR21', 'PTGS2', 'FLJ43681', 'SOWAHC', 'ECH1', 'FOS', 'KLF2', 'TNFAIP1', 'LINC00473', 'BTBD19', 'CCNL1', 'OSTCP1', 'ISG20', 'CTGF', 'ZP3']
    g15do = ['ASB4', 'NRL', 'MESP1', 'MGC45922', 'CALHM2', 'KAT2A', 'RAG2', 'THYN1', 'GBX2', 'PSRC1', 'ERBB3', 'A2M-AS1', 'AKIP1', 'ATP6V1F', 'MRGPRF', 'LOC143666', 'GJB2', 'HOXA7', 'LOC389791', 'APOBEC3H', 'MIS18A', 'DNAJC28', 'IGFBP6', 'TRAPPC5', 'FCRLB', 'LAMTOR2', 'MGC3771', 'TMEM19', 'S1PR5', 'LOC257396', 'WDR87', 'MESP2', 'C19orf45', 'HOTTIP', 'MPZ', 'LOC100131347', 'LOC100288748', 'NPDC1', 'HILPDA', 'P2RY1', 'ZNF503-AS2', 'SAP30', 'KLLN', 'LOC440700', 'DPY19L2P3', 'LOC201617', 'TMEM160', 'FAM178A', 'RDM1', 'ATP12A', 'PRSS36', 'PLK1', 'RFXAP', 'EPN3', 'LOC100996455', 'CBY3', 'ZNRF2P2', 'METTL7B', 'RIMBP3B', 'LOC283214', 'PIGV', 'GJB5', 'HLTF-AS1', 'AQP10', 'ZNF599', 'FZD2', 'MDFIC', 'C18orf1', 'CDCA3', 'ENPP4', 'ST3GAL6-AS1', 'IDH1-AS1', 'THRA', 'HOXA9', 'LOC154822', 'SP9', 'RAD21-AS1', 'FLJ33534', 'RHPN1-AS1', 'RIMBP3C', 'SIX1', 'TMEM105', 'PCDHAC1', 'ISCA2', 'KRT222', 'STOML1', 'ZIC2', 'NUDT6', 'LGALS8-AS1', 'HIST1H2BM', 'PRR15', 'HOXA-AS3', 'CCDC115']
    g45up = ['IL24', 'CREM', 'ACTBL2', 'IL24', 'GPR126', 'LINC00312', 'SCHIP1', 'KRT81', 'LAMC1', 'ATP2B1', 'RBKS', 'MACF1', 'C3orf32', 'CREB5', 'SEMA7A', 'MFI2', 'CCL20', 'ACTG2', 'PPAP2B', 'HIVEP2', 'SMS', 'NF2', 'SLC20A2', 'GRK5', 'GPR126', 'ELOVL5', 'AGPAT9', 'MUC2', 'RCAN1', 'MB21D2', 'WWC1', 'KRTAP4-9', 'NT5E', 'NHS', 'PSG4', 'COL12A1', 'RXFP1', 'LINC00319', 'NEDD4L', 'CNN1', 'NF2', 'BACH2', 'SYNJ2', 'SPAG9', 'KRTAP2-3', 'RHOH', 'RXFP1', 'IL1RAP', 'LINC00592', 'FN1', 'SLC9A2', 'NF2', 'LINC00473', 'SCHIP1', 'MAP3K14-AS1', 'PMEPA1', 'PMEPA1', 'BCAR3', 'SCHIP1', 'ZSWIM6', 'SAMD4A', 'TNC', 'SPAG9', 'CREB5', 'MICAL2', 'RXFP1', 'WDR25', 'NEDD4L', 'CLU', 'C19orf71', 'FERMT2', 'CREB5', 'NF2', 'RIMKLB', 'AKAP12', 'SOGA2', 'C10orf55', 'ATP13A3', 'PSG4', 'LOXHD1', 'ACTG2', 'AP3M1', 'TES', 'NF2', 'SH3KBP1', 'MUM1L1', 'CAV3', 'NPAS2', 'TRIM55', 'AGPAT9', 'CPA4', 'GPR126', 'TNFAIP3', 'TNFAIP3', 'FLJ42393', 'AP3S1', 'SLC20A2', 'C3orf32', 'RXFP1', 'RXFP1', 'KCTD20', 'C3orf32', 'RXFP1', 'COBL', 'IQCJ-SCHIP1', 'PSG8', 'KIAA0226', 'TRIM55', 'IL24', 'PSG8', 'LOC100130880', 'CMIP', 'RXFP1', 'NF2', 'MGLL', 'PAPL', 'COTL1', 'IL1RAP', 'ATXN7', 'SMS', 'WWC1', 'PSG5', 'PLS3', 'CAV3', 'WWC1', 'NT5E', 'GPR126', 'CFLAR-AS1', 'SDC4', 'IFLTD1', 'AP3M1', 'TRIM55', 'ITGBL1', 'TRIM55', 'BCAR3', 'PMEPA1', 'PLEC', 'CREB5', 'NEDD4L', 'CELA2B', 'RTN4', 'NF2', 'IL1RAP', 'SAMD4A', 'EXT1', 'NEDD4L', 'RXFP1', 'ARSJ', 'ABHD2', 'SEMA7A', 'ETS1', 'FERMT2', 'DGKD', 'RTN4', 'ATXN7', 'LBH', 'IL1RAP', 'NCEH1', 'LOC100505583', 'PLEC', 'MGLL', 'KLHL30', 'ESYT2', 'SH3KBP1', 'PLAU', 'LOC100652768', 'KIAA0513', 'TAGLN', 'CTNNAL1', 'IL1RAP', 'MYL7', 'SNORD114-28', 'MIR3918', 'MUM1L1', 'IL1RAP', 'PSG9', 'RXFP1', 'LOC283403', 'MN1', 'ABL2', 'STK24', 'C3orf32', 'LCP1', 'PLEC', 'MIR661', 'SGMS2', 'SGMS2', 'ABL2', 'CRIM1', 'NKD2', 'TNFAIP3', 'IL18', 'PLS3', 'PFKP', 'SNORD114-27', 'ALOX5AP', 'CORO1C', 'RTN4', 'PSG8', 'RXFP1', 'C8orf86', 'FBLIM1', 'SEMA7A', 'NEDD4L', 'NF2', 'CCL20', 'NCEH1', 'BACH2', 'HSPB8', 'SCHIP1', 'NCEH1', 'CXCL2', 'LINC00602', 'NCEH1', 'PLEC', 'RXFP1', 'CRYGN', 'RTN4', 'RXFP1', 'PLAU', 'F3', 'LMCD1', 'PMEPA1', 'PDE4D', 'NEDD4L', 'PADI1', 'RTN4', 'SH3KBP1', 'BCAR3', 'PSG5', 'GJC2', 'TECTA', 'CLU', 'PTPN1', 'PMEPA1', 'EMR3', 'IQCJ-SCHIP1', 'DEC1', 'IL18', 'CEBPE', 'CELA2A', 'CPA4', 'KRTAP4-12', 'ITGA6', 'SLC20A2', 'F3', 'ITGA6', 'NKD2', 'CMIP', 'PDLIM2', 'PLEC', 'SPAG9', 'TAGLN', 'FERMT2', 'KRT17', 'NF2', 'NUAK2', 'AGPAT9', 'DKFZp434J0226', 'ATP2B1', 'MGLL', 'CT62', 'NEDD4L', 'GBP1', 'RXFP1', 'WDR25', 'IL24', 'PLS3', 'ABL2', 'MYO16', 'ETS1', 'AKAP12', 'LURAP1L', 'DGKD', 'SPAG9', 'MIR630', 'TRIO', 'COL12A1', 'ABHD2', 'MYO16']
    g45do = ['THRA', 'NFATC4', 'NDRG2', 'TRIB2', 'APOF', 'PDCD4', 'HMMR', 'FGF20', 'NFATC4', 'SPRY1', 'HMMR', 'YPEL3', 'MPEG1', 'HFE', 'HFE', 'HOXA-AS3', 'YPEL2', 'ORAI3', 'PDK2', 'METTL21CP1', 'BTN3A2', 'S1PR5', 'PTCH1', 'C3orf18', 'BTN3A2', 'KLRC3', 'HFE', 'CDNF', 'CXXC4', 'PTTG1', 'UNKL', 'NDRG2', 'LOC100129722', 'KLRC3', 'OAS1', 'THRA', 'NDRG2', 'C9orf173', 'ABTB1', 'PCOLCE-AS1', 'FAM100B', 'YPEL3', 'PCOLCE', 'HOXA10', 'NFATC4', 'HIST1H2BF', 'PTCH1', 'ZNF792', 'HFE', 'BTN3A2', 'SMTNL1', 'KLRC4', 'ASPM', 'TCP11L2', 'CCDC74B', 'CCNG2', 'NKAPP1', 'GRB7', 'NDRG2', 'NFATC4', 'MGAT3', 'BTN3A1', 'MLN', 'C9orf173', 'IL22RA1', 'FRAT1', 'SEPP1', 'ITGB7', 'SEPP1', 'NDRG2', 'ARL6IP5', 'ASPM', 'PDCD4', 'PDK2', 'HMMR', 'BTN3A2', 'BTBD8', 'ICAM4', 'PDCD4', 'ERAP2', 'FLJ37035', 'KLRC2', 'C3orf18', 'AXIN2', 'ICAM4', 'PIF1', 'SEPT5', 'OAS1', 'TOP2A', 'ICAM4', 'CTDSP2', 'LOC100506368', 'TRIB2', 'HFE', 'SPRY1', 'NFATC4', 'KLHL24', 'OAS1', 'CCDC74B', 'CYP39A1', 'CKAP2L', 'HOXA10', 'SEPP1', 'GRB7', 'GRB7', 'NDRG2', 'TMEM244', 'C3orf18', 'C9orf173', 'HFE', 'CCDC152', 'UCP1', 'PHOX2A', 'ERAP2', 'HELB', 'TSHZ1', 'ABTB1', 'BTN3A1', 'C9orf173', 'HFE', 'BMF', 'PTCH1', 'MLN', 'NDRG2', 'LOC100129213', 'PDK2', 'HMMR', 'PDK2', 'ABTB1', 'BTN3A1', 'FLJ30403', 'MLN', 'PTCH1', 'EPHB6', 'HFE', 'BTN3A1', 'C9orf173', 'HFE', 'SPRY1', 'FRAT2', 'HOXA10-HOXA9', 'BTN3A2', 'CRIPT', 'NDRG2', 'C9orf173', 'PTCH1', 'GRB7']
    folder = ['/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/0.genes.bed.count.bed.DMSOCAnascent.res.txt','/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/15.genes.bed.count.bed.DMSOCAnascent.res.txt','/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/45.genes.bed.count.bed.DMSOCAnascent.res.txt']
    savedir = savedir = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/'
    run(g15up, g15do, g45up, g45do)
