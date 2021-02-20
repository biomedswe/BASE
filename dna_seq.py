from misc import *


import multiprocessing # lets you find out how many threads your cpu has



from os import listdir, sys, mkdir, getenv, path
# listdir lets you list the files in a directory, sys lets you exit program with sys.exit()

import subprocess
# lets you run shell commands

import argparse
# shows help sections in shell

import csv
# lets you create a tsv-file


path_reads  = getenv("HOME")+"/sequencing_project/dna_seq/reads/"
path_reference = getenv("HOME")+"/sequencing_project/reference_genome/GRCh38.p13.genome.fa"
path_exclude_template = getenv("HOME")+"/sequencing_project/excludeTemplate/human.hg38.excl.tsv"
configManta_path = getenv("HOME")+"anaconda3/envs/sequencing/bin/configManta.py"
runWorkflow_path = getenv("HOME")+"dna_seq/Manta/runWorkflow.py"

import time

# Global shortcuts to output paths
output_aligned = getenv("HOME")+"/sequencing_project/dna_seq/Aligned/"
output_sorted = getenv("HOME")+"/sequencing_project/dna_seq/Sorted/"
output_merged = getenv("HOME")+"/sequencing_project/dna_seq/Merged/"
output_removedDuplicates = getenv("HOME")+"/sequencing_project/dna_seq/Removed_duplicates/"
output_realigned = getenv("HOME")+"/sequencing_project/dna_seq/Realigned/"
output_haplotypecaller = getenv("HOME")+"/sequencing_project/dna_seq/GATK_haplotypecaller/"
output_delly = getenv("HOME")+"/sequencing_project/dna_seq/Delly/"
output_manta = getenv("HOME")+"/sequencing_project/dna_seq/Manta/"
# print(output_aligned, output_sorted, output_merged, output_removed_duplicates, output_realigned,output_haplotypecaller, output_delly, output_manta)
# time.sleep(3)

# Global variables for output lists
align_list = f"{output_aligned}align.txt"
sort_list = f"{output_sorted}sort.txt"
merge_list = f"{output_merged}merge.txt"
removeDuplicates_list = f"{output_removedDuplicates}remove_duplicate.txt"
haplotypeCaller_list = f"{output_haplotypecaller}haplotypeCaller.txt"
realign_list = f"{output_realigned}realign.txt"
validateBam_file = f"{output_aligned}validateBam.file"
haplotypeCaller_file = f"{output_haplotypecaller}haplotypeCaller.file"
delly_file = f"{output_delly}delly.file"
manta_file = f"{output_manta}manta.file"





'''Optional:

You may input a .bed file which states the genome interval you want to use for the reference genome in the SNV calling using GATK HaplotypeCaller.
This makes the script go faster compared to searching through the whole reference genome.
'''
############################################################# Functions ################################################################################

def validate_id(options):
    '''This function validates wether the tumor_id and normal_id you entered is present in your reads'''

    if options.tumor_id and options.normal_id not in "".join(listdir(path_reads)):
        print("You have entered a tumor_id and/or normal_id that is not present in your reads.\n")
        print("Please restart program and verify that you have typed in the right \"clinical_id\" for both tumor (-t) and normal (-n)! ")
        input("Press any key to exit program")
        sys.exit()
    else:
        print("tumor_id and normal_id correctly validated!\n")



def dna_menu():
    '''Prints an introduction to the program where you can choose either to create a library list file or run analysis'''


    info_script()
    print("This script makes an automated analysis of WGS/WES reads.\nBefore running analysis, you must have indexed the reference genome and also have a list file containing all your reads.\nTo complete these requirements, please run option 1 and 2 first\n")
    print("1. Index reference genome")
    print("2. Create library list file")
    print("3. Run analysis")
    choices = 3
    text = "(leave blank to exit)"
    choice = validate_choice(choices,text)
    return choice


def index_reference():
    print("1. Index reference genome\n")
    bwa_index = f"bwa index {path_reference}"
    # subprocess.run(bwa_index, shell=True)
    create_dict = f"samtools dict {path_reference} -o {path_reference[:-2]}dict"
    # subprocess.run(create_dict, shell=True)
    create_fai = f"samtools faidx {path_reference}"
    # subprocess.run(create_fai, shell=True)
    print("\nIndexing reference genome completed!\n")
    input("Press any key to return to DNA-analysis menu...")


def build_library(options):
    '''This function lists all WGS files in the directory and writes them to a library-list text file.
       You can create either a file for single-end protocol or paired-end protocol.
       This file is then used in the dna analysis'''

    while True:
        print("2. Create library list file.\n\n")
        print("Single end sequencing (1)\n")
        print("Paired end sequencing (2)\n")
        choices = 2
        text = "Choose 1 or 2 (leave blank to return to DNA-analysis menu)"
        choice = validate_choice(choices, text)

        if choice == "":
            return ""

        elif confirm_choice():
            clear_screen()
            files = listdir(path_reads)
            with open("dna_seq/library.txt", 'w') as out_file:
                for line, library_id in enumerate(files, start=1):

                    if choice == '1': # Single-end sequencing
                        if options.tumor_id in library_id:
                            out_file.write(f"{options.tumor_id} {library_id[:-6]} {library_id} N/A\n")
                        else:
                            out_file.write(f"{options.normal_id} {library_id[:-6]} {library_id} N/A\n")

                    elif choice == '2': # Paired-end sequencing
                        if (line % 2) == 1: # not even
                            if options.tumor_id in library_id:
                                out_file.write(f"{options.tumor_id} {library_id[:-8]} {library_id} ")
                            else:
                                out_file.write(f"{options.normal_id} {library_id[:-8]} {library_id} ")
                        else: # even
                            out_file.write(f"{library_id}\n")
            print("Library list file created!\nNow that you have created your list library file, you can run the analysis!\n")
            input("Press any key to return to DNA-analysis menu...")
            return ''

        else:
            clear_screen()
            continue




def validate_bam(): # add align.txt
    '''This function runs picard ValidateSamFile to check if any errors are present in the aligned files.
       returns True if no errors are found or False if errors are found'''

    if path.isfile(f'{output_aligned}validateBam.file'):
        print("Validating bam allready completed, skips step...\n")
        return True

    else:
        print("Validating .bam files...\n")
        for file in align_list:
            validate = f"picard ValidateSamFile -I dna_seq/Aligned/{file} -MODE SUMMARY"
            process = subprocess.run(validate, shell=True, capture_output=True, text=True)
            if "No errors found" in process.stdout:
                pass
            else:
                return False
        with open(validateBam_file, 'a'):
            pass
        print("Validation completed without errors, continuing with next step...")
        return True





def alignment():
    '''This function loops through a library list containing either single-end or paired-end protocol, it will automatically detect what protocol it is.
       For every loop, the bwa mem command will be run and the output name for each run will be saved to a txt-file for use in the next step if the command finnish without errors.'''


    if path.isfile(align_list):
        print("Burrown Wheeler aligner allready completed, skips step... ")

    else:
        while True:
            print("3. Run analysis\n\n\n")
            print("Burrows Wheeler aligner\n")
            # Checks nr of available threads on pc and asks for input
            threads_text = ("Input number of CPU threads you want to use (Avaliable threads: " + str(multiprocessing.cpu_count()) + "): ")
            available_threads = multiprocessing.cpu_count()
            threads = validate_choice(available_threads, threads_text)


            if threads == "":
                return ""

            elif confirm_choice():
                clear_screen()
                print(f"Aligning reads using {threads} CPU threads out of {available_threads}...\n")
                time.sleep(1)

                # create output directory
                try:
                    mkdir(output_aligned)
                except FileExistsError:
                    pass





                with open('dna_seq/library.txt', 'r') as fastq_list:
                    for line in fastq_list.readlines():
                        clinical_id, library_id, read1, read2 = line.split()
                        read_group_header = f'@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}'

                        if read2 == 'N/A': # single-end
                            cmd_bwa = f"1: bwa mem -R {read_group_header} {path_reference} {path_reads}/{read1} -t {threads} -o {output_aligned}"
                        else: # paired-end
                            cmd_bwa = f"bwa mem -R {read_group_header} {path_reference} {path_reads}/{read1} {path_reads}/{read2} -t {threads} -o {output_aligned}"
                        run_command(cmd_bwa, 'Alignment with Burrows Wheeler aligner')
                        create_outputList(align_list, f"{library_id}.bam")
                print("Burrows Wheeler aligner completed!")
                return

            else:
                clear_screen()
                continue



def sort(options):
    '''This function reads the completed_steps.txt to check if the previous step was completed without errors.

       The function returns a list of string containg the filenames of the sorted tumor samples and normal samples separated.
       The string is needed because you have several inputs in the next function and can therefore not run a for loop'''

    if path.isfile(sort_list):
        print("Picard sortsam allready completed, skips step... ")

    else:
        print("\nSorting SAM/BAM files using Picard Sortsam...\n")
        tumor_sort_str = ""
        normal_sort_str = ""
        write_to_file = ""
        try:
            # create target directory
            mkdir(output_sorted)
        except FileExistsError:
            pass

        with open(align_list, 'r') as list:
            for sample in list.readlines():
                cmd_sortsam = f"picard SortSam -I dna_seq/Aligned/{sample} -O dna_seq/Sorted/{sample} -SORT_ORDER coordinate"
                run_command(cmd_sortsam, 'Sorting files with Picard SortSam')
                if options.tumor_id in sample:
                    tumor_sort_str += f" -I dna_seq/Sorted/{sample}".rstrip()
                else:
                    normal_sort_str += f" -I dna_seq/Sorted/{sample}".rstrip()

            write_to_file = tumor_sort_str.lstrip() + '\n' + normal_sort_str.lstrip()
            print(write_to_file)
            create_outputList(sort_list, write_to_file)

            print("\nPicard SortSam completed!\n")



def merge(options):
    '''This function merges all the input files in the sort_list to one output file'''

    if path.isfile(merge_list):
        print("Picard MergeSamFiles allready completed, skips step... ")

    else:
        print("\nMerging SAM/BAM files using Picard MergeSamFiles...\n")

        try:
            # create target directory
            mkdir(output_merged)
        except FileExistsError:
            pass

        with open(sort_list, 'r') as list:
            for sample in list.readlines():
                if options.tumor_id in sample.lstrip():
                    cmd_merge = f"picard MergeSamFiles {sample} -O dna_seq/Merged/{options.tumor_id}.bam"
                    run_command(cmd_merge, f"Merging {sample} with Picard MergeSamFiles")
                    create_outputList(merge_list, f"{options.tumor_id}.bam")
                else:
                    cmd_merge = f"picard MergeSamFiles {sample} -O dna_seq/Merged/{options.normal_id}.bam"
                    run_command(cmd_merge, f"Merging {sample} with Picard MergeSamFiles")
                    create_outputList(merge_list, f"{options.normal_id}.bam")

            print("\nPicard MergeSamFiles completed!\n")


def remove_duplicate():
    '''This function removes duplicates '''

    if path.isfile(removeDuplicates_list):
        print("Picard MarkDuplicates allready completed, skips step... ")

    else:

        print("\nRemoving duplicates in SAM/BAM files using Picard MarkDuplicates...\n")

        try:
            # create target directory
            mkdir(output_removedDuplicates)
        except FileExistsError:
            pass

        with open(merge_list, 'r') as list:
            for sample in list.readlines():
                cmd_rd = f"picard MarkDuplicates {sample} -O {output_removedDuplicates}{sample} -M {output_removedDuplicates}marked_dup_metrics_{sample}.txt"
                run_command(cmd_rd, f'Removing duplicates in {sample.rstrip()} with Picard MarkDuplicates')
                create_outputList(removeDuplicates_list, f"{sample.rstrip()}")

            print("\nPicard MarkDuplicates completed!\n")


def realign():
    '''This function realigns the bam files'''

    if path.isfile(realign_list):
        print("GATK LeftAlignIndels allready completed, skips step... ")

    else:
        print("\nRealigning SAM/BAM files using GATK LeftAlignIndels...\n")

        try:
            # create target directory
            mkdir(output_realigned)
        except FileExistsError:
            pass

        with open(removeDuplicates_list, 'r') as list:
            for sample in list.readlines():
                cmd_index = f"samtools index {sample}"
                run_command(cmd_index, f'Indexing {sample.rstrip()} with samtools index')
                cmd_leftAlignIndels = f"gatk LeftAlignIndels -R {path_reference} -I {sample} -O dna_seq/Realigned/{sample}"
                run_command(cmd_leftAlignIndels, f'Realigning {sample.rstrip()} with GATK LeftAlignIndels')
                create_outputList(realign_list, f"{sample.rstrip()}")

            print("\nGATK LeftAlignIndels completed!\n")


def gatk_snv(options):

    if path.isfile(haplotypeCaller_file):
        print("GATK haplotypeCaller allready completed, skips step... ")

    else:
        print("\nLooking for SNV's using GATK HaplotypeCaller...\n")

        try:
            # create target directory
            mkdir(output_haplotypecaller)
        except FileExistsError:
            pass

        with open(realign_list, 'r') as list:
            sample_1, sample_2 = list.readlines()
            if options.intervals:
                cmd_call = f"gatk HaplotypeCaller -R {path_reference} -I {sample_1.rstrip()} -I {sample_2.rstrip()} -O {output_haplotypecaller}{options.tumor_id}.vcf -L {options.intervals}"
            else:
                cmd_call = f"gatk HaplotypeCaller -R {path_reference} -I {sample_1.rstrip()} -I {sample_2.rstrip()} -O {output_haplotypecaller}{options.tumor_id}.vcf"
            run_command(cmd_call, f'Looking for SNV\'s in {sample_1.rstrip()} and {sample_2.rstrip()} with GATK haplotypeCaller')


            # Remove all reads with read depth less than 10
            cmd_filter_read_depth = f"bcftools view -i 'MIN(FMT/DP)>10' dna_seq/GATK_haplotypecaller/2064-01.vcf > dna_seq/GATK_haplotypecaller/2064-01_filter_RD10.vcf"
            subprocess.run(cmd_filter_read_depth, shell=True)
            with open(haplotypeCaller_file, 'a'):
                pass
            print("\nGATK HaplotypeCaller completed!\n")


def delly(options):
    '''This function creates an output directory and runs delly to call for somatic SNV's'''

    if path.isfile(delly_file):
        print("Delly allready completed, skips step... ")

    else:
        print("\nLooking for somatic SNV's using delly\n")

        try:
            # create target directory
            mkdir(output_delly)
        except FileExistsError:
            pass

        with open(realign_list, 'r') as list:
            sample_1, sample_2 = list.readlines()
            cmd_delly_call = f"delly call -x {path_exclude_template} -g {path_reference} -o {output_delly}delly.bcf {sample_1} {sample_2}"
            run_command(cmd_delly_call, f'Looking for somatic SNV\'s in {sample_1.rstrip()} and {sample_2.rstrip()} using delly')


            with open(f'{output_delly}sample.tsv', 'w', newline='') as tsv:
                tsv_output = csv.writer(tsv, delimiter='\t')
                tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
                tsv_output.writerow([f"{options.normal_id}", 'control'])

            cmd_dos2unix = f"dos2unix {output_delly}sample.tsv"
            run_command(cmd_dos2unix, '')

            cmd_filter = f"delly filter -f somatic -o {output_delly}delly_filter.bcf -s {output_delly}sample.tsv {output_delly}delly.bcf"
            run_command(cmd_filter, 'Filtering variants...')

            cmd_convert = f"bcftools view {output_delly}delly_filter.bcf > {output_delly}delly_filter.vcf"
            run_command(cmd_convert, '')
            with open(delly_file, 'a'):
                pass
            print("\nDelly variant caller completed!\n")


def manta():
    '''This function creates an output directory and runs manta to call for somatic SNV's'''

    if path.isfile(manta_file):
        print("Manta allready completed, skips step... ")

    else:
        print("\nLooking for somatic SNV's using manta\n")
        try:
            # create target directory
            mkdir(output_manta)
        except FileExistsError:
            pass

        with open(realign_list, 'r') as list:
            sample_1, sample_2 = list.readlines()

            cmd_create_config_file = f"{configManta_path} --tumorBam={sample_1.rstrip()} --bam={sample_2.rstrip()} --referenceFasta={path_reference} --runDir={output_manta}"
            run_command(cmd_create_config_file, '')

            cmd_runWorkflow = f"{runWorkflow_path} -m local -j 4"
            run_command(cmd_runWorkflow, '')

        with open(manta_file, 'a'):
            pass
        print("\nManta variant caller completed!\n")
