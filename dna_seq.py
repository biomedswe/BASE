from misc import *

from os import listdir, sys, mkdir, getenv
# listdir lets you list the files in a directory, sys lets you exit program with sys.exit()

import subprocess
# lets you run shell commands

import argparse
# shows help sections in shell

import csv
# lets you create a tsv-file

global path_reads
global path_reference
path_reads  = getenv("HOME")+"/sequencing_project/dna_seq/reads/"
path_reference = getenv("HOME")+"/sequencing_project/reference_genome/GRCh38.p13.genome.fa"


global path_exclude_template
path_exclude_template = subprocess.run(["find", ".", "-type", "f", "-name", "human.hg38.excl.tsv"], capture_output=True, text=True).stdout.strip()




'''INSTRUCTIONS
------------------------------------------------------------------------------------------------------------------------------

The script will automatically:

Create the folders: Aligned, Sorted, Merged, Removed_duplicates, Realigned for the output files to be placed in.

Create a library of the reads in a text file which is used by the program to run the analysis. Will be placed in the same folder as this script.


type script.py -h for help

-------------------------------------------------------------------------------------------------------------------------------
Required for the script to work properly:

- A reference genome in fasta format. E.g. human_g1k_v37.fasta

- reads in fastq format

- Input clinical ID of tumor and normal samples for the program to distinguish between the samples. E.g. 2064_01 (t) and 987_02 (n)

- Each line in the library list file should look like this: <clinical_id library_id read> (single-end) or <clinical_id library_id read1 read2> (paired-end)
  Make sure it is correct by checking the library.txt. Check build_library() function if something is wrong.

  Example: clinical_id = 2064-01 (tumor), 987-02 (normal). library_id = whole filename except .fastq


--------------------------------------------------------------------------------------------------------------------------------
Optional:

You may input a .bed file which states the genome interval you want to use for the reference genome in the SNV calling using GATK HaplotypeCaller.
This makes the script go faster compared to searching through the whole reference genome.
'''
############################################################# Functions ################################################################################

def validate_id(options):
    '''This function validates wether the tumor_id and normal_id you entered is present in your reads'''

    if options.tumor_id and options.normal_id not in "".join(listdir(path_reads)):
        print("\nYou have entered a tumor_id and/or normal_id that is not present in your reads")
        input("Press any key to exit program")
        sys.exit()
    else:
        print("\ntumor_id and normal_id correctly validated!")



def print_intro():
    '''Prints an introduction to the program where you can choose either to create a library list file or run analysis'''

    print("\nThis script makes an automated analysis of WGS reads.\nIf you dont have an indexed reference genome or a file containing all the names of the reads,\nplease run option 1 and 2 first\n")
    print("1. Index reference genome")
    print("2. Create library list file")
    print("3. Run analysis")


    return input("(leave blank to exit): ")

def index_reference():
    bwa_index = f"bwa index {path_reference}"
    # subprocess.run(bwa_index, shell=True)
    create_dict = f"samtools dict {path_reference} -o {path_reference[:-2]}dict"
    subprocess.run(create_dict, shell=True)
    print("\nIndexing reference genome completed...\n")
    input("Press any key to return to main menu...")


def build_library(options):
    '''This function lists all WGS files in the directory and writes them to a library-list text file.
       You can create either a file for single-end protocol or paired-end protocol.
       This file is then used in the dna analysis'''

    print("\nsingle end sequencing (1)\nPaired end sequencing (2)\n")

    choice = int(input("Choose (1) or (2): "))


    if  0 < choice < 3:
        files = listdir(path_reads)
        with open("library.txt", 'w') as out_file:
            for line, library_id in enumerate(files, start=1):

                if choice == 1: # Single-end sequencing
                    if options.tumor_id in library_id:
                        out_file.write(f"{options.tumor_id} {library_id[:-6]} {library_id} N/A\n")
                    else:
                        out_file.write(f"{options.normal_id} {library_id[:-6]} {library_id} N/A\n")

                elif choice == 2: # Paired-end sequencing
                    if (line % 2) == 1: # not even
                        if options.tumor_id in library_id:
                            out_file.write(f"{options.tumor_id} {library_id[:-8]} {library_id} ")
                        else:
                            out_file.write(f"{options.normal_id} {library_id[:-8]} {library_id} ")
                    else: # even
                        out_file.write(f"{library_id}\n")
        print("\nLibrary list file created!\nNow that you have created your list library file, you can run the analysis!\n")
        input("Press any key to return to main menu...")



    else:
        print("Invalid choice, try again!")
        build_library()





def validate_bam(align_list):
    '''This function runs picard ValidateSamFile to check if any errors are present in the aligned files.
       returns True if no errors are found or False if errors are found'''

    # print("Validating .bam files...\n")
    # for file in align_list:
    #     validate = f"java -jar {path_picard} ValidateSamFile -I ./Aligned/{file} -MODE SUMMARY"
    #     process = subprocess.run(validate, shell=True, capture_output=True, text=True)
    #     print(process.stdout)
    #     if "No errors found" in process.stdout:
    #         pass
    #     else:
    #         return False

    return True





def alignment():
    '''This function loops through a library list containing either single-end or paired-end protocol, it will automatically detect what protocol it is.

       For every loop, the bwa mem command will be run and the library_name for each run will be saved to a list which is returned for use in the next step.'''

    print("\nAligning reads using Burrows Wheeler aligner...\n")
    try:
        # create output directory
        mkdir("Aligned")
    except FileExistsError:
        pass


    # Creates an empty list to store the names of the output files
    align_list = []

    with open('library.txt', 'r') as fastq_list:
        for line in fastq_list.readlines():
            clinical_id, library_id, read1, read2 = line.split()
            if read2 == 'N/A': # single-end
                run_bwa = f"bwa mem -R '@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}' {path_reference} {path_reads}/{read1} -t 12 | samtools view -bS -o Aligned/align_{library_id}.bam -"
                # subprocess.run(run_bwa, shell=True)
                align_list.append(f"align_{library_name}.bam")
            else: # paired-end
                run_bwa = f"bwa mem -R '@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}' {path_reference} {path_reads}/{read1} {path_reads}/{read2} -t 12 | samtools view -bS -o Aligned/align_{library_id}.bam -"
                # subprocess.run(run_bwa, shell=True)
                align_list.append(f"align_{library_id}.bam")

    print("\nBurrows wheeler aligner completed!\n")

    return align_list



def sort(options, align_list):
    '''This function loops through the list of aligned samples and sort them by executing the picard sortsam
       algoritm as a shell command.
       The function returns a list of string containg the filenames of the sorted tumor samples and normal samples separated.
       The string is needed because you have several inputs in the next function and can therefore not run a for loop'''

    print("\nSorting SAM/BAM files using Picard Sortsam...\n")
    sort_list = []
    tumor_sort_str = ""
    normal_sort_str = ""
    try:
        # create target directory
        mkdir("Sorted")
    except FileExistsError:
        pass

    for sample in align_list:
        cmd = f"picard SortSam -I Aligned/{sample} -O Sorted/sort_{sample} -SORT_ORDER coordinate"
        subprocess.call(cmd, shell=True)
        if options.tumor_id in sample:
            tumor_sort_str += f"-I Sorted/sort_{sample} "
        else:
            normal_sort_str += f"-I Sorted/sort_{sample} "

    sort_list.extend((tumor_sort_str, normal_sort_str))

    print("\nPicard SortSam completed!\n")
    return sort_list


def merge(sort_list, options):
    '''This function merges all the input files in the sort_list to one output file'''

    print("\nMerging SAM/BAM files using Picard MergeSamFiles...\n")
    merge_list = []
    try:
        # create target directory
        mkdir("Merged")
    except FileExistsError:
        pass


    for sample_input in sort_list:
        if options.tumor_id in sample_input:
            cmd = f"picard MergeSamFiles {sample_input} -O Merged/merge_{options.tumor_id}.bam"
            # subprocess.run(cmd, shell=True)
            merge_list.append(f"-I Merged/merge_{options.tumor_id}.bam")

        else:
            cmd = f"picard MergeSamFiles {sample_input} -O Merged/merge_{options.normal_id}.bam"
            # subprocess.run(cmd, shell=True)
            merge_list.append(f"-I Merged/merge_{options.normal_id}.bam")



    print("\nPicard MergeSamFiles completed!\n")
    return merge_list

def remove_duplicate(merge_list):
    '''This function removes duplicates '''

    print("\nRemoving duplicates in SAM/BAM files using Picard MarkDuplicates...\n")
    rd_list = []
    try:
        # create target directory
        mkdir("Removed_duplicates")
    except FileExistsError:
        pass

    for sample_input in merge_list:
        cmd = f"picard MarkDuplicates {sample_input} -O Removed_duplicates/rd_{sample_input[10:]} -M Removed_duplicates/marked_dup_metrics_{sample_input[10:-4]}.txt"
        # subprocess.run(cmd, shell=True)
        rd_list.append(f"Removed_duplicates/rd_{sample_input[10:]}")
    print("\nPicard MarkDuplicates completed!\n")
    print(rd_list)
    return rd_list

def realign(rd_list):
    print("\nRealigning SAM/BAM files using GATK LeftAlignIndels...\n")
    realign_output = []
    try:
        # create target directory
        mkdir("Realigned")
    except FileExistsError:
        pass

    for sample_input in rd_list:
        index = f"samtools index {sample_input}"
        # subprocess.run(index, shell=True)
        cmd = f"{path_gatk} LeftAlignIndels -R {path_reference} -I {sample_input} -O Realigned/realigned_{sample_input[28:]}"
        # subprocess.run(cmd, shell=True)
        realign_output.append(f"Realigned/realigned_{sample_input[28:]}")
    print("\nGATK LeftAlignIndels completed!\n")
    return realign_output

def snv_calling(options, realign_output):
    print("\nLooking for SNV's using GATK HaplotypeCaller...\n")

    try:
        # create target directory
        mkdir("GATK_haplotypecaller")
    except FileExistsError:
        pass

    if options.intervals:
        snv = f"{path_gatk} HaplotypeCaller -R {path_reference} -I {realign_output[0]} -I {realign_output[1]} -O GATK_haplotypecaller/{options.tumor_id}.vcf -L {options.intervals}"
        # subprocess.run(snv, shell=True)
    else:
        snv = f"{path_gatk} HaplotypeCaller -R {path_reference} -I {realign_output[0]} -I {realign_output[1]} -O GATK_haplotypecaller/{options.tumor_id}.vcf"
        # subprocess.run(snv, shell=True)
    print("\nGATK HaplotypeCaller completed!\n")


def delly(options, realign_output):
    '''This function creates an output directory and runs delly to call for somatic SNV's'''

    print("\nLooking for somatic SNV's using delly\n")

    try:
        # create target directory
        mkdir("Delly")
    except FileExistsError:
        pass

    cmd_call = f"delly call -x {path_exclude_template} -g {path_reference} -o Delly/delly.bcf {realign_output[0]} {realign_output[1]}"
    # subprocess.run(cmd_call, shell=True)


    with open('Delly/sample.tsv', 'w', newline='') as tsv:
        tsv_output = csv.writer(tsv, delimiter='\t')
        tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
        tsv_output.writerow([f"{options.normal_id}", 'control'])

    cmd_dos2unix = "dos2unix Delly/sample.tsv"
    subprocess.run(cmd_dos2unix, shell=True)

    cmd_filter = "delly filter -f somatic -o Delly/delly_filter.bcf -s Delly/sample.tsv Delly/delly.bcf"
    subprocess.run(cmd_filter, shell=True)

    cmd_convert = "bcftools view Delly/delly_filter.bcf > Delly/delly_filter.vcf"
    subprocess.run(cmd_convert, shell=True)



def manta():
    pass

# ----------------------------------------------------- Main program starts -----------------------------------------------------------------
def main():

    # metavar makes the help text more tidy
    parser = argparse.ArgumentParser(description='''Automated GDC DNA-seq analysis script for 1 patient sample at a time. Before running, see instructions in dna_seq.py for more information''')
    parser.add_argument("-t", "--tumor_id", metavar="", required=True, help="Input clinical id of tumor samples")
    parser.add_argument("-n", "--normal_id", metavar="", required=True, help="Input clinical id of normal samples")
    parser.add_argument("-i", "--intervals", metavar="", required=False, help="Input path to reference-genome interval if you have any (for use in SNV calling)")
    # parser.add_argument("-R", "--reference_genome", metavar="", required=True, help="Input name of the reference genome, eg. GRCh38.p13.genome.fa")
    # parser.add_argument("-r", "--reads", metavar="", required=True, help="Input name of directory containing sequencing reads")
    options = parser.parse_args() # all arguments will be passed to the functions


    validate_id(options)



    start_choice = print_intro() # return choice
    while start_choice != "":
        if start_choice == '1':
            if confirm_choice():
                index_reference()
                start_choice = print_intro()
            else:
                start_choice = print_intro()


        elif start_choice == '2':
            if confirm_choice():
                build_library(options)
                start_choice = print_intro()
            else:
                start_choice = print_intro()

        elif start_choice == '3':
            # GDC DNA-Seq analysis pipeline
            if confirm_choice():
                align_list = alignment()
                if validate_bam(align_list):
                    print("Validation completed without errors")
                    sort_list = sort(options, align_list)
                    # merge_list = merge(sort_list, options)
                    # rd_list = remove_duplicate(merge_list)
                    # realign_output = realign(rd_list)
                    # snv_calling(options, realign_output)
                    # delly(options, realign_output)
                else:
                    print("An error was found, terminating process!")
                    start_choice = ''
            else:
                start_choice = print_intro()



        else:
            print("Invalid choice, please try again!")
            start_choice = print_intro()

    print("exiting program...")
    sys.exit()






if __name__ == '__main__':
    main()
