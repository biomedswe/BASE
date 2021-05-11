import pandas as pd
import argparse



def filter_csv(options):
    df = pd.read_csv(options.input)

    minjun_genes = {"DHRS3",
                    "LINC00339",
                    "PTP4A2",
                    "YARS",
                    "TRIT1",
                    "ZMPSTE24",
                    "PDE4B",
                    "NTRK1",
                    "USF1",
                    "PIGC",
                    "ZNF672",
                    "ZNF692",
                    "CMPK2",
                    "LBH",
                    "CEBPZ",
                    "AHSA2",
                    "SFXN5",
                    "TMEM131",
                    "ARHGEF4",
                    "ITGA6",
                    "ATG9A",
                    "ITM2C",
                    "NCL",
                    "THUMPD3",
                    "IQSEC1",
                    "OXNAD1",
                    "NISCH",
                    "MAGEF1",
                    "AC093323.3",
                    "C4orf3",
                    "LSM6",
                    "SLC12A7",
                    "UGT3A2",
                    "WDR41",
                    "SKP1",
                    "SERPINB6",
                    "DTNBP1",
                    "ABT1",
                    "TRIM27",
                    "PBX2",
                    "PEX6",
                    "FAM26F",
                    "TAB2",
                    "EZR",
                    "RNASET2",
                    "CARD11",
                    "NPY",
                    "KBTBD2",
                    "VOPP1",
                    "GTF2I",
                    "NCF1",
                    "POMZP3",
                    "GIGYF1",
                    "TNPO3",
                    "LUC7L2",
                    "PRSS1",
                    "PRSS3P1",
                    "ZYX",
                    "BLK",
                    "SORBS3",
                    "RMDN1",
                    "UQCRB",
                    "MRPL13",
                    "KHDRBS3",
                    "RPL8",
                    "PRSS3",
                    "UBE2R2",
                    "UNC13B",
                    "CCDC107",
                    "SH3GLB2",
                    "RAPGEF1",
                    "RP11-473E2.2",
                    "ASB13",
                    "IL2RA",
                    "APBB1IP",
                    "HNRNPF",
                    "HNRNPH3",
                    "CALHM2",
                    "SMNDC1",
                    "MGMT",
                    "HBB",
                    "TMEM109",
                    "CTSW",
                    "FOLR3",
                    "NCAPD2",
                    "GAPDH",
                    "CD69",
                    "PTGES3",
                    "LYZ",
                    "CHFR",
                    "EXOSC8",
                    "NUDT15",
                    "LPAR6",
                    "NGDN",
                    "IPO4",
                    "TM9SF1",
                    "GZMB",
                    "WARS",
                    "RCOR1",
                    "MTA1",
                    "CHP1",
                    "OIP5-AS1",
                    "SQRDL",
                    "RP11-358M11.3",
                    "ANXA2",
                    "RCN2",
                    "SLC9A3R2",
                    "RPS15A",
                    "CLN3",
                    "CD19",
                    "RFWD3",
                    "EMC8",
                    "C1QBP",
                    "COX11",
                    "SUPT4H1",
                    "UBE2O",
                    "DCXR",
                    "DSC3",
                    "TRAPPC8",
                    "INO80C",
                    "TPGS2",
                    "FAM69C",
                    "PRTN3",
                    "CFD",
                    "GNA15",
                    "APBA3",
                    "STXBP2",
                    "KANK2",
                    "PRKCSH",
                    "CCDC61",
                    "ALDH16A1",
                    "ZNF83",
                    "ZNF444",
                    "ZNF134",
                    "TMEM230",
                    "GINS1",
                    "NINL",
                    "RPN2",
                    "BPI",
                    "CD40",
                    "BIRC7",
                    "SLC2A4RG",
                    "AGPAT3",
                    "SNAP29",
                    "CABIN1",
                    "ASPHD2",
                    "TRABD",
                    "NCAPH2",
                    "SLC25A5",
}


    # create an empty dict
    genes = {}
    jonas_genes1 = []
    jonas_genes2 = []



    # counters for validation
    tot = 0
    sig = 0

    # assign geneName as dict key and 2 empty placeholders for counting total genes and significant genes
    for index, row in df.iterrows():
        genes[row.loc['geneName']] = [0,0]

    # iterate trough rows in df. if conditions are fullfilled for that geneName, add +1 for total count and +1 for count that fullfills condition, i.e., is significant.
    # else just add total count for that geneName

    # dict will show: geneName : {[total counts, matching counts]}

    for index, row in df.iterrows():
        if row.loc['pValue_WGS'] <= float(options.pvalue) and row.loc['pValue_CNV'] <= float(options.pvalue) and not float(options.lower_foldchange) < row.loc['VAF_ratio_WGS'] < float(options.upper_foldchange) and not float(options.lower_foldchange) < row.loc['VAF_ratio_CNV'] <  float(options.upper_foldchange):
            genes[row.loc['geneName']][0] += 1
            genes[row.loc['geneName']][1] += 1
            tot += 1
            sig += 1
        else:
            genes[row.loc['geneName']][0] += 1
            tot += 1

    # print(f"tot: {tot} sig: {sig}")

    for index, row in df.iterrows():
        if genes[row.loc['geneName']][1] > (genes[row.loc['geneName']][0]/2):
            jonas_genes1.append(row.loc['geneName'])
            if row.loc['geneName'] not in jonas_genes2:
                jonas_genes2.append(row.loc['geneName'])


    print("minjun: 184, unique:", len(minjun_genes))
    print("jonas:", len(jonas_genes1), len(jonas_genes2))
    print(len(minjun_genes.intersection(jonas_genes2)))

    # for index, row in df.iterrows():
    #     # if not significant. (if significant gene entries are not bigger than 50% of total entries)
    #     if not genes[row.loc['geneName']][1] > (genes[row.loc['geneName']][0]/2): # if not significant > total/2
    #         df.drop(index, inplace=True)





    # df.drop_duplicates(subset ="geneName", inplace=True)
    # print(df)
    # print filtered file to csv
    # df.to_csv(options.output, sep=',', index=False)


def main():
    # argparse lets ju input arguments to the script before starting it
    parser = argparse.ArgumentParser(description='''This script is used to filter out genes with significant ASE''')
    parser.add_argument("-i", "--input", metavar="", required=True, help="Enter input file")
    parser.add_argument("-o", "--output", metavar="", required=True, help="Enter output file")
    parser.add_argument("-p", "--pvalue", metavar="", required=True, help="Enter threshold pValue")
    parser.add_argument("-l", "--lower_foldchange", metavar="", required=True, help="Enter lower threshold foldchange")
    parser.add_argument("-u", "--upper_foldchange", metavar="", required=True, help="Enter upper threshold foldchange")
    options = parser.parse_args() # all arguments can be called by options. e.g. options.input
    filter_csv(options)

if __name__ == '__main__':
    main()
