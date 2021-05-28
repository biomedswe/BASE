import pandas as pd
import argparse



def filter_csv(options):
    df = pd.read_csv(options.input)




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
    print("tot:", tot, "sig:", sig)


    # Thiw is just used to compare my gene list to minjuns earlier results
    for index, row in df.iterrows():
        if genes[row.loc['geneName']][1] >= (genes[row.loc['geneName']][0]/2):
            jonas_genes1.append(row.loc['geneName'])
            if row.loc['geneName'] not in jonas_genes2:
                jonas_genes2.append(row.loc['geneName'])




    # overlap = 0
    # print("jonas:", len(jonas_genes1), len(jonas_genes2))
    # with open("HeH_2064-01.txt", "r") as minjun:
    #     minjun_list = minjun.read().splitlines()
    #     for gene in minjun_list:
    #         if gene in jonas_genes2:
    #             overlap += 1
    # print(overlap, "/",len(minjun_list))

    for index, row in df.iterrows():
        # if not significant. (if significant gene entries are not bigger than 50% of total entries)
        if not genes[row.loc['geneName']][1] >= (genes[row.loc['geneName']][0]/2): # if not significant > total/2
            df.drop(index, inplace=True)





    # df.drop_duplicates(subset ="geneName", inplace=True)
    # print(df)
    # print filtered file to csv
    df.to_csv(options.output, sep=',', index=False)


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
