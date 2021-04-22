import pandas as pd
from os import sys
import argparse
import numpy as np


def filter_csv(options):
    df = pd.read_csv(options.input)
    variants = {}

    # create a dict with all genes ass key
    for index, row in df.iterrows():
        variants[row.loc['geneName']] = [0,0]

    # iterate trough rows in df. if conditions are fullfilled, add the variantType as key and +1 for total variant and +1 for those that fullfill condition.
    # dict will show: variantType : {[total counts, matching counts]}
    for index, row in df.iterrows():
        if row.loc['pValue_WGS_data'] <= float(options.pvalue) and row.loc['pValue_CNV_data'] <= float(options.pvalue) and not 0.5 < row.loc['VAF_ratio_WGS_data'] < 2 and not 0.5 < row.loc['VAF_ratio_CNV_data'] < 2:
            variants[row.loc['geneName']][0] += 1
            variants[row.loc['geneName']][1] += 1
        else:
            variants[row.loc['geneName']][0] += 1


    for index, row in df.iterrows():
        # if not significant. (significant gene entries are not equal or bigger than the total entries)
        if not variants[row.loc['geneName']][1] >= (variants[row.loc['geneName']][0]/2): # if not significant >= total
            print(index+2, 'not significant:', row.loc['geneName'])
            df.drop(index, inplace=True)
        else:
            print(index+2, 'significant:', row.loc['geneName'])

    # df.dropna(inplace=True)




    # print filtered file to csv
    df.to_csv(options.output, sep=',', index=False)


def main():
    parser = argparse.ArgumentParser(description='''This script is used to filter out genes with significant ASE''')
    parser.add_argument("-i", "--input", metavar="", required=True, help="Enter input file")
    parser.add_argument("-o", "--output", metavar="", required=True, help="Enter output file")
    parser.add_argument("-p", "--pvalue", metavar="", required=True, help="Enter threshold pValue")
    options = parser.parse_args() # all arguments will be passed to the functions
    filter_csv(options)

if __name__ == '__main__':
    main()
