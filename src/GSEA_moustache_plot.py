__author__ = 'Jonathan Rubin'

import matplotlib.pyplot as plt

def run(GSEA_file_up,GSEA_file_dn):
    path = list()
    x = list()
    y = list()
    upx = list()
    upy = list()
    serx = list()
    sery = list()
    with open(GSEA_file_up) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split(',')
            name = line[0]
            NES = float(line[5])
            FDR = float(line[7])
            if 'SERUM' in name or 'EGF' in name:
                serx.append(NES)
                sery.append(FDR)
            if FDR < 0.1:
                upx.append(NES)
                upy.append(FDR)
            else:
                x.append(NES)
                y.append(FDR)

    dnx = list()
    dny = list()
    with open(GSEA_file_dn) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split(',')
            name = line[0]
            try:
                NES = float(line[5])
            except:
                print line
            FDR = float(line[7])
            if 'OMA' in name or 'CANCER' in name or 'METASTASIS' in name:
                print name
                serx.append(NES)
                sery.append(FDR)
            if FDR < 0.1:
                dnx.append(NES)
                dny.append(FDR)
            else:
                x.append(NES)
                y.append(FDR)
                
    F = plt.figure() 
    ax = F.add_subplot(111)
    plt.scatter(x,y,c='k',edgecolor="",s=40)
    plt.scatter(dnx,dny,c='g',edgecolor="",s=40)
    plt.scatter(upx,upy,c='r',edgecolor="",s=40)
    # plt.scatter(serx,sery,c='y',edgecolor="",s=40,alpha=0.9)
    ax.set_title('GSEA Chemical and Genetic Perturbations gene sets',fontsize=16)
    ax.set_ylabel('FDR q-val',fontsize=18)
    ax.set_xlabel('NES',fontsize=18)
    plt.axhline(0.1, color='black',linestyle='dashed')
    plt.axvline(0, color='black',linestyle='dashed')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.savefig('../figures/GSEA_moustache_DMSOCAt45.svg')
        

    


if __name__ == "__main__":
    GSEA_file_up = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/DMSOCAt45_CGP_GSEA_up.csv'
    GSEA_file_dn = '/Users/jonathanrubin/Google Drive/Colorado University/Jonathan/SerumResponseCA_REP1GROSEQ/DMSOCAt45_CGP_GSEA_dn.csv'

    run(GSEA_file_up,GSEA_file_dn)
