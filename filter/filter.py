import pandas as pd
import argparse



def filter_csv(options):
    '''Final filter of csv after "ASEReadCounter" and "add_wgs_data_to_csv" functions. You can specify the pValue and the lower and upper threshold for RNA/DNA ratio for both WGS and CNV data'''
    df = pd.read_csv(options.input)


    # create a dict with all genes as key
    variants = {}
    for index, row in df.iterrows():
        variants[row.loc['geneName']] = [0,0]

    # iterate trough rows in df. if conditions are fullfilled, add the geneName as key and +1 for total variant and +1 for those that fullfill condition, i.e., significant.
    # dict will show: variantType : {[total counts, matching counts]}
    for index, row in df.iterrows():
        if row.loc['pValue_WGS_VAF'] <= float(options.pvalue) and row.loc['pValue_CNV'] <= float(options.pvalue) and not float(options.lower_foldchange) < row.loc['RNA/DNA_ratio_WGS_VAF'] < float(options.upper_foldchange) and not float(options.lower_foldchange) < row.loc['RNA/DNA_ratio_CNV'] <  float(options.upper_foldchange):
            variants[row.loc['geneName']][0] += 1
            variants[row.loc['geneName']][1] += 1
        else:
            variants[row.loc['geneName']][0] += 1


    for index, row in df.iterrows():
        # if not significant. (if significant gene entries are not bigger than 50% of total entries)
        if not variants[row.loc['geneName']][1] > (variants[row.loc['geneName']][0]/2): # if not significant > total/2
            df.drop(index, inplace=True)

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
