__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import os

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def run(file1,protein_sequences,TF1,TF2,figuredir):
    d = dict()
    d2 = dict()
    sequences = dict()
    with open(protein_sequences) as F:
        for line in F:
            line = line.strip().split()
            protein = line[1].split('_')[0]
            sequence = line[-2]
            if protein in sequences:
                print 'True'
            sequences[protein] = sequence


    with open(file1) as F:
        F.readline()
        F.readline()
        for line in F:
            line = line.strip().split()
            key = line[-4]
            for TF in TF1:
                if TF == key.split(';')[0]:
                    print TF, line[3]+line[14], line[1]
                    # if float(line[1]) > 1.1 or float(line[1]) < 0.9:
                    # if TF in d:
                    #     d[TF].append((float(line[1]),list(set(line[14].split(';'))),sequences[TF],line[3]))
                    # else:
                    #     d[TF] = [(float(line[1]),list(set(line[14].split(';'))),sequences[TF],line[3])]

            for TF in TF2:
                if TF == key.split(';')[0]:
                    print TF, line[3]+line[14], line[1]
                    # if float(line[1]) > 1.1 or float(line[1]) < 0.9:
                    # if TF in d2:
                    #     d2[TF].append((float(line[1]),list(set(line[14].split(';'))),sequences[TF],line[3]))
                    # else:
                    #     d2[TF] = [((float(line[1])),list(set(line[14].split(';'))),sequences[TF],line[3])]

    c = mcolors.ColorConverter().to_rgb
    rvb = make_colormap(
        [c('green'), c('blue'), 0.5, c('blue'), c('red'), 1.0, c('red')])
    val = [1 for i in range(len(d)+len(d2))]    # the bar lengths
    gap = 4
    offset = 1.5
    positions = [i*gap for i in range(-len(d2),len(d)+1)]
    positions.remove(0)  # the bar centers on the y axis
    x = list()
    y = list()
    colors = list()

    F = plt.figure()
    ax = F.add_subplot(111)
    ax.barh(positions,val, align='center',height=1,color='white',fc=(1, 1, 1, 0))
    # ax.barh(positions,val, align='center',height=1,color='white',fc=(0, 0, 0, 0.5),left=val)
    labels = list()
    i = len(d)
    fontsize = 8
    for TF in TF1:
        if TF in d:
            labels.append(TF)
            ax.text(0,i*gap-offset,'1',fontsize = 10)
            ax.text(1,i*gap-offset,str(len(d[TF][0][2])),fontsize = 10)
            for val,p,seq,a in d[TF]:
                for pos in p:
                    pos = int(pos)
                    # if seq[pos-1] == 'S' or seq[pos-1] == 'T':
                    pos = float(pos)
                    length = float(len(seq))
                    x.append(pos/length)
                    y.append(i*gap)
                    colors.append(val)
                    if 0.75 > val or val > 1.25:
                        if TF == 'MYC':
                            print a + str(int(pos))
                            if a + str(int(pos)) == 'S62':
                                ax.text(pos/length-0.01,i*gap+1.7,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                            elif a + str(int(pos)) == 'T58':
                                ax.text(pos/length-0.01,i*gap+1,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                            else:
                                ax.text(pos/length-0.01,i*gap+1,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                        elif TF == 'JUNB':
                            if a + str(int(pos)) == 'S259':
                                ax.text(pos/length-0.01,i*gap+1.7,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                            elif a + str(int(pos)) == 'T58':
                                ax.text(pos/length-0.01,i*gap+1,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                            else:
                                ax.text(pos/length-0.01,i*gap+1,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                        else:
                            ax.text(pos/length-0.01,i*gap+1,a + str(int(pos)),fontsize=fontsize)
            i = i-1

    i = -1
    for TF in TF2:
        if TF in d2:
            labels.append(TF)
            ax.text(0,i*gap-offset,'1',fontsize = 10)
            ax.text(1,i*gap-offset,str(len(d2[TF][0][2])),fontsize=10)
            boolean = True
            for val,p,seq,a in d2[TF]:
                for pos in p:
                    pos = int(pos)
                    # if seq[pos-1] == 'S' or seq[pos-1] == 'T':
                    pos = float(pos)
                    length = float(len(seq))
                    x.append(pos/length)
                    y.append(i*gap)
                    colors.append(val)
                    if 0.75 > val or val > 1.25:
                        if TF == 'MKL1':
                            if seq[int(pos)-1] + str(int(pos)) == 'S454':
                                if boolean:
                                    ax.text(pos/length-0.03,i*gap+1.5,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                                boolean = False
                            elif seq[int(pos)-1] + str(int(pos)) == 'S511':
                                ax.text(pos/length-0.03,i*gap-1,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                            else:
                                ax.text(pos/length-0.03,i*gap+1,seq[int(pos)-1] + str(int(pos)),fontsize=fontsize)
                        else:
                            ax.text(pos/length-0.01,i*gap+1,a + str(int(pos)),fontsize=fontsize)
            i = i-1
    plt.scatter(x,y,zorder=1,s=60,edgecolor="",c=colors, cmap=rvb,vmin=0.5,vmax=1.5,alpha=0.7)
    ax.set_title('Detected phosphorylation sites')
    ax.set_xlabel('Relative Amino Acid Position')
    plt.yticks(positions, list(reversed(labels)))
    ax.set_xticks([])
    ax.set_xlim([-0.1,1.1])
    ax.set_ylim([min(positions)-gap,max(positions)+gap])
    plt.colorbar()
    # F.savefig(figuredir + '/Phosphoproteomics_FOSL1_JUND_Interactors.png', dpi=1200)
    # plt.show()

def find_connectors(folder,file2):
    outfile = open(folder + 'network.txt','w')
    d = dict()
    for file1 in os.listdir(folder):
        if 'string_interactions' in file1 and ('SRF' in file1 or 'JUN' in file1):
            key = file1.split('_')[-1].split('.')[0]
            d[key] = list()
            with open(folder+file1) as F:
                for line in F:
                    line = line.strip().split()
                    if key in line[0]:
                        d[key].append(line[1])
                    elif key in line[1]:
                        d[key].append(line[0])

    print d['JUN']
    print d['SRF']
    d2 = dict()
    with open(file2) as F:
        F.readline()
        F.readline()
        for line in F:
            line = line.strip().split()
            phosphopeptide = line[-4].split(';')[0]
            for TF in d:
                for protein in d[TF]:
                    if phosphopeptide == protein and float(line[1]) < 1:
                        if 'CDK8/19' not in d2:
                            d2['CDK8/19'] = list()
                        if protein not in d2:
                            d2[protein] = list()
                        if protein in d2['CDK8/19']:
                            index = d2['CDK8/19'].index(protein)+1
                            if float(d2['CDK8/19'][index]) > float(line[1]):
                                d2['CDK8/19'][index] = line[1]
                        else:
                            d2['CDK8/19'].append(protein)
                            d2['CDK8/19'].append(line[1])
                        d2[protein].append(TF)
    for source in d2:
        if source == 'CDK8/19':
            for i in range(0,len(d2[source]),2):
                outfile.write(source + '\t' + d2[source][i] + '\t' + d2[source][i+1] + '\n')
        elif len(d2[source]) > 0:
            d2[source] = list(set(d2[source]))
            for target in d2[source]:
                outfile.write(source + '\t' + target + '\t' + '1.0' + '\n')






if __name__ == "__main__":
    file1 = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/misc/HCT116_Serinduction_phospho_rep1.txt'
    protein_sequences = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/GROAnalysis/files/uniprot-proteome%3AUP000005640.tab'
    # TF1 = ['RARG', 'HES7', 'PAX5', 'RARB', 'TLX1', 'ID4', 'THA', 'THA', 'PITX2', 'AP2D', 'THB', 'HEY2', 'ZN219', 'IKZF1', 'P73', 'RXRB', 'TBX2', 'HEY1', 'PAX5', 'HES5', 'RARA']
    # TF2 = ['FOSL1','JUN','PITX3','FOSL2','BATF','FOSB','FOS','JUND']
    # TF1 = ['ETV2', 'SP3', 'KLF15', 'RARG', 'PAX5', 'ITF2', 'SNAI2', 'RARB', 'TLX1', 'ID4', 'RFX5', 'LMX1A', 'THA', 'VDR', 'THA', 'CPEB1', 'SRF', 'FUBP1', 'PITX2', 'IRX2', 'AP2D', 'FOXM1', 'FOXJ2', 'THB', 'TFDP1', 'FOXO1', 'FOXG1', 'HMGA1', 'SRY', 'ZN219', 'FOXJ3', 'FOXJ3', 'ZEB1', 'IKZF1', 'HTF4', 'P73', 'FOXL1', 'MAFA', 'SP4', 'RXRB', 'TBX2', 'MAZ', 'TBX5', 'SP2', 'FOXD3', 'KLF16', 'BPTF', 'HEY1', 'PAX5', 'HEY2', 'NKX25', 'MESP1', 'ONEC2', 'TGIF1', 'RARA']
    # TF2 = ['JUN', 'FOSL1', 'FOSL2', 'MAFG', 'BATF', 'NF2L1', 'JUND']
    # TF1 = ['FOS', 'ATF3', 'FOSL2', 'FOSL1', 'SMAD3', 'MAPK8', 'EP300', 'MAPK9', 'MAPK10', 'ATF2']
    # TF2 = ['MKL1', 'MYOCD', 'ELK4', 'FOS', 'NKX2-5', 'KS6A1', 'GATA4', 'HOPX', 'TAGLN', 'KS6A3']
    # TF1 = ['AHR', 'BACH1', 'TF7L2', 'EMX2', 'ALX4', 'BATF3', 'MLX', 'NFKB1', 'MEIS3', 'SOX3', 'E2F1', 'E2F4', 'STA5A', 'AIRE', 'MEF2C', 'TGIF2', 'PO6F1', 'HXD4', 'HLF', 'ATOH1', 'ZN423', 'PBX2', 'SP1', 'MEIS2', 'ESX1', 'MNT', 'SP1', 'DLX5', 'HNF4G', 'TBX15', 'HMX2', 'COT2', 'RXRA', 'GLI1', 'RORG', 'ETV2', 'COT2', 'JUN', 'CEBPD', 'RHXF1', 'NR2F6', 'GLI2', 'ETV5', 'RUNX2', 'PO3F4', 'PO4F2', 'ZN282', 'OVOL1', 'GSX2', 'PDX1', 'DLX3', 'ISX', 'TFEB', 'CREB1', 'BCL6B', 'HNF6', 'PRDM1', 'NDF2', 'MEOX1', 'DLX1', 'TBX1', 'ZN713', 'GLIS2', 'HXD12', 'PITX3', 'GFI1B', 'HXB2', 'KLF3', 'MAF', 'MSX1', 'ZFX', 'LHX9', 'GLIS1', 'NKX61', 'SP3', 'KLF15', 'RARG', 'NFAC2', 'GABPA', 'FOSL1', 'SOX2', 'TBP', 'HES7', 'CXXC1', 'ENOA', 'SPI1', 'NR6A1', 'HNF1A', 'GMEB2', 'SMAD3', 'IRF1', 'BC11A', 'HXB7', 'GSC2', 'GABP1', 'BSH', 'LHX4', 'NR4A2', 'HMBX1', 'SOX13', 'NR2E1', 'NR2C2', 'ALX1', 'STAT4', 'HNF1B', 'PROX1', 'ARX', 'NFIA', 'NFAC4', 'TEAD4', 'BARH1', 'NFIA', 'HMX3', 'HXD3', 'PROP1', 'FIGLA', 'PRRX1', 'RUNX3', 'FOXF2', 'PAX5', 'CR3L2', 'NKX21', 'SHOX', 'MEIS1', 'SOX7', 'FOSL2', 'TYY1', 'P5F1B', 'HSF2', 'THAP1', 'HXD13', 'NR1H4', 'TBX21', 'ZKSC1', 'HSF1', 'VSX1', 'DUXA', 'MTF1', 'ZBT49', 'FOXF1', 'ZBT18', 'MZF1', 'PO2F3', 'ZBT7A', 'ITF2', 'ARI5B', 'SNAI2', 'NRF1', 'RARB', 'CDX2', 'NKX31', 'CREB3', 'MNX1', 'LEF1', 'SPZ1', 'TLX1', 'LHX6', 'CEBPA', 'TLX1', 'ZN232', 'ID4', 'PRD14', 'SOX1', 'NKX28', 'EMX1', 'RFX5', 'FOXO4', 'TAL1', 'LBX2', 'STAT1', 'NR5A2', 'ZN384', 'VENTX', 'MYOG', 'EVX2', 'STAT1', 'TAL1', 'P53', 'GCR', 'EVI1', 'ZBT7B', 'LMX1A', 'GCR', 'IRF5', 'MCR', 'REST', 'OTX2', 'REL', 'HEN1', 'GCM2', 'MIXL1', 'TF2LX', 'NFYB', 'VAX2', 'AP2A', 'FOXH1', 'ZIC3', 'ZIC2', 'NRL', 'MEF2B', 'NFIC', 'E2F8', 'GFI1', 'PO2F1', 'SPIB', 'SNAI1', 'BHE23', 'FLI1', 'LHX3', 'MGAP', 'VDR', 'RORA', 'THA', 'ETS1', 'FOXQ1', 'BCL6', 'BHE41', 'PPARD', 'HXB3', 'VDR', 'THA', 'NR0B1', 'ZN143', 'CPEB1', 'MAFG', 'HIF1A', 'E2F2', 'SOX11', 'MAFG', 'SRF', 'FUBP1', 'BATF', 'SPDEF', 'HXA9', 'NR1H2', 'TF65', 'SOX21', 'ZN589', 'NR1I2', 'NR1I3', 'PITX2', 'IRX2', 'PO3F3', 'TEAD3', 'ETV7', 'NR1I2', 'NR1I3', 'ELF3', 'ATF6A', 'PHX2B', 'PAX1', 'NF2L1', 'PO5F1', 'PPARG', 'AP2B', 'NKX62', 'NFE2', 'AP2D', 'BHE40', 'SUH', 'ZSC16', 'FOXO3', 'PAX7', 'STA5B', 'SOX5', 'DDIT3', 'ETV3', 'NOTO', 'PO2F2', 'THB', 'ZEP2', 'SMAD2', 'E4F1', 'FOXM1', 'FOXJ2', 'MYC', 'RELB', 'THB', 'SOX15', 'HXD9', 'TFDP1', 'RFX3', 'GSX1', 'HXC8', 'ARNT2', 'FOSB', 'TFDP1', 'TEF', 'P63', 'COT1', 'HMGA2', 'BHE22', 'FOXO1', 'NFAC3', 'PLAG1', 'STAT3', 'ONEC3', 'ZSCA4', 'MSX2', 'PLAG1', 'FOXG1', 'SOX10', 'SOX8', 'HXA5', 'CEBPB', 'GLIS3', 'SRBP1', 'E2F3', 'NFAC1', 'DLX2', 'NFKB2', 'PO6F2', 'DLX4', 'CLOCK', 'HIC2', 'NFAC1', 'USF2', 'LMX1B', 'FOS', 'HMGA1', 'OTX1', 'ZBTB6', 'HOMEZ', 'OLIG2', 'E2F5', 'FOXC2', 'KLF8', 'HXA2', 'WT1', 'PLAL1', 'ERF', 'SCRT1', 'SRY', 'JDP2', 'ZIC4', 'CREB5', 'GBX2', 'HESX1', 'ELF1', 'HSFY1', 'NFAT5', 'RREB1', 'NKX22', 'DMBX1', 'TCF7', 'FOXA3', 'FOXA1', 'EGR3', 'MUSC', 'GATA4', 'ZN219', 'SOX17', 'CRX', 'HXD10', 'HBP1', 'GSC', 'PKNX1', 'ESR2', 'CEBPE', 'RFX1', 'FOXP3', 'HXB8', 'MEF2D', 'FOXP2', 'ATF7', 'RFX4', 'ESR2', 'ARI3A', 'MECP2', 'SCRT2', 'FOXJ3', 'ATF3', 'UNC4', 'SMRC1', 'FOXC1', 'ZFHX3', 'ARI3A', 'PRDM4', 'FOXJ3', 'ZIC1', 'NOBOX', 'ZEB1', 'NF2L2', 'ELK1', 'IKZF1', 'USF1', 'PEBB', 'ETV1', 'GLI3', 'MLXPL', 'NR1D1', 'ISL1', 'CUX1', 'RUNX1', 'HNF4A', 'TFE2', 'INSM1', 'HINFP', 'HES1', 'GBX1', 'FOXD1', 'HTF4', 'PAX8', 'VAX1', 'FOXI1', 'BARH2', 'P63', 'HXC11', 'RX', 'TF7L1', 'IRX3', 'KAISO', 'CR3L1', 'PITX1', 'ERR2', 'ZN350', 'KAISO', 'MAFF', 'TFE3', 'NR2E3', 'TEAD1', 'P73', 'PBX3', 'FOXL1', 'TYY2', 'ZN410', 'TWST1', 'MAFA', 'P73', 'MEF2A', 'PAX6', 'E2F7', 'ANDR', 'EGR4', 'VSX2', 'KLF13', 'KLF12', 'SP4', 'HMX1', 'PO4F1', 'PAX3', 'SRBP2', 'ZN740', 'RXRB', 'TBX2', 'COE1', 'NR4A1', 'TFAP4', 'HXC12', 'PURA', 'BATF', 'HXA13', 'MYBB', 'NKX32', 'MYBA', 'SOX18', 'COT1', 'KLF4', 'BRAC', 'PRGR', 'HXA11', 'AP2C', 'MYOD1', 'SOX9', 'ERR1', 'CEBPG', 'PRGR', 'HXB1', 'MEOX2', 'EHF', 'HME2', 'TBX19', 'BARX1', 'TFCP2', 'EHF', 'EGR2', 'FOXK1', 'GATA2', 'ELK3', 'HXA7', 'CEBPZ', 'MAZ', 'ATF2', 'HXC10', 'SMAD4', 'HXD11', 'PKNX2', 'RXRG', 'XBP1', 'NR2C1', 'ELF5', 'PAX4', 'UBIP1', 'ZN639', 'MAX', 'HXA1', 'HLTF', 'TBX5', 'FOXO6', 'STAT6', 'RARG', 'CENPB', 'FOXA2', 'TBR1', 'GATA1', 'GATA6', 'STAT2', 'CREM', 'EGR1', 'NANOG', 'GATA1', 'HXA10', 'JUNB', 'RAX2', 'EPAS1', 'ZN784', 'SP2', 'E2F6', 'JUND', 'FOXD3', 'ZEP1', 'DRGX', 'MYF6', 'NDF1', 'KLF16', 'ZN333', 'MITF', 'SHOX2', 'CDC5L', 'RFX2', 'HXC6', 'SOX4', 'CUX2', 'MBD2', 'DBP', 'PRRX2', 'GATA6', 'NFIL3', 'ELK4', 'IRF3', 'SMAD1', 'PAX2', 'TBX4', 'ETV4', 'HXB6', 'IRF9', 'PAX2', 'PHX2A', 'ETV6', 'HSF4', 'NGN2', 'BPTF', 'SPIC', 'ETS2', 'HXB13', 'LHX8', 'MAFK', 'ALX3', 'HEY1', 'MAFB', 'ARNT', 'EGR1', 'FEV', 'PAX5', 'MAFK', 'ZBTB4', 'HES5', 'NFYC', 'IRF7', 'ELF2', 'ZBTB4', 'ERR3', 'HAND1', 'CTCF', 'PPARG', 'STF1', 'OLIG3', 'KLF14', 'ASCL2', 'PO4F3', 'EVX1', 'TBX3', 'MYB', 'HXC13', 'DLX6', 'CTCFL', 'PPARA', 'HXD8', 'PTF1A', 'TBX20', 'PPARA', 'HEY2', 'NKX23', 'ZN148', 'NKX25', 'MYCN', 'PO3F1', 'MESP1', 'GATA3', 'PBX1', 'ESR1', 'HME1', 'FOXD2', 'CDX1', 'EOMES', 'PO3F2', 'ESR1', 'DPRX', 'GCM1', 'BRCA1', 'OLIG1', 'BMAL1', 'GRHL1', 'ONEC2', 'NFYA', 'NANOG', 'IRF4', 'ZN652', 'FOXB1', 'BARX2', 'IRF8', 'ZKSC3', 'YBOX1', 'ERG', 'IRF2', 'ISL2', 'ATF1', 'PIT1', 'TGIF1', 'TGIF1', 'KLF1', 'RARA', 'NR4A3', 'BHA15', 'GATA5', 'HIC1', 'RARA', 'KLF6', 'LHX2', 'ZN524']
    # TF2 = []
    # TF1 = ['SP1', 'SP1', 'THA', 'SP2', 'RARG', 'PAX5', 'VDR', 'THA', 'CPEB1', 'PITX2', 'TFDP1', 'FOXO1', 'HMGA1', 'THB', 'RARB', 'IKZF1', 'SP4', 'RXRB', 'TBX2', 'KLF4', 'MAZ', 'SP3', 'PAX5', 'RARA']
    # TF2 = ['FOSL1', 'PAX4', 'LMX1A', 'RX', 'JUND']
    # TF1 = ['LMNB1','UBN1','ASF1A','HMGA2','CABIN1','EP400','TP53','RB1','MYC','RPS6KB1']
    # TF2 = []
    # TF1 = ['SP1','HDAC2','CDKN1A', 'SMAD4', 'SUMO2', 'NFATC1','UBE2I','SUMO1','PIAS1']
    # TF2 = []
    TF1 = ['MMP1','MMP9','PLAUR','EP300','JUN','JUND','MYC','CCND1','JUNB','USF2']
    TF2 = ['CDK4','ATF2','MAPK9','FOS','MAPK10','MAPK8','FOSL1','FOSB','FOSL2','EP300','FOS']
    figuredir = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/GROAnalysis/figures/'
    run(file1,protein_sequences,TF1,TF2,figuredir)

    # folder = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/'
    # find_connectors(folder,file1)
