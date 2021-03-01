# Packages used in script
from os import listdir, sys, mkdir, getenv, path
import subprocess
import argparse
import csv
import multiprocessing # lets you find out how many threads your cpu has
import time
from misc import * # my own package with miscellaneous functions

class DnaSeqAnalysis():

    def __init__(self):

        # Shortcuts to files used in DNA sequencing analysis
        self.reference_genome_file = getenv("HOME")+"/sequencing_project/reference_genome/GRCh38.p13.genome.fa"
        self.reference_genome_exclude_template_file = getenv("HOME")+"/sequencing_project/excludeTemplate/human.hg38.excl.tsv"
        self.configManta_file = getenv("HOME")+"/anaconda3/envs/sequencing/bin/manta-1.6.0.release_src/src/python/bin/configManta.py"
        self.runWorkflow_file = getenv("HOME")+"/sequencing_project/dna_seq/Manta/runWorkflow.py"


        # Shortcuts to input folders
        self.reference_genome_dir = getenv("HOME")+"/sequencing_project/reference_genome/"
        self.dna_reads_dir  = getenv("HOME")+"/sequencing_project/dna_seq/reads/"


        # Shortcuts to general folder
        self.dna_seq_dir = getenv("HOME")+"/sequencing_project/dna_seq/"

        # Shortcuts to output folder
        self.aligned_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Aligned/"
        self.sorted_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Sorted/"
        self.merged_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Merged/"
        self.removed_duplicates_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Removed_duplicates/"
        self.realigned_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Realigned/"
        self.haplotypecaller_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/GATK_haplotypecaller/"
        self.delly_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Delly/"
        self.manta_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/Manta/"
        self.manta_results_dir = getenv("HOME")+"/sequencing_project/dna_seq/Manta/results/variants/"

        # Shortcuts to output lists (used for input in pipeline steps) (and also for validation if pipeline step i allready completed)
        self.alignedFiles_list = f"{self.aligned_output_dir}alignedFiles.txt"
        self.sortedFiles_list = f"{self.sorted_output_dir}sortedFiles.txt"
        self.mergedFiles_list = f"{self.merged_output_dir}mergedFiles.txt"
        self.removeDuplicates_list = f"{self.removed_duplicates_output_dir}remove_duplicate.txt"
        self.realignedFiles_list = f"{self.realigned_output_dir}realignedFiles.txt"

        # Shortcuts to files used to validate if pipeline step is allready completed
        self.index_reference_genome_complete = f"{self.reference_genome_dir}index.complete"
        self.validate_bam_complete = f"{self.aligned_output_dir}validateBam.complete"
        self.haplotypecaller_complete = f"{self.haplotypecaller_output_dir}haplotypeCaller.complete"
        self.delly_complete = f"{self.delly_output_dir}delly.complete"
        self.manta_complete = f"{self.manta_output_dir}manta.complete"


    #---------------------------------------------------------------------------
    def create_directory(path):
        '''This function creates the output directories for the different analysis steps'''
        try:
            mkdir(path)
        except FileExistsError:
            pass

    #---------------------------------------------------------------------------
    def step_completed(self, file, text):
        '''This function checks if an analysis step is completed by looking after a created file *.complete'''
        if path.isfile(file):
            print(text)
            return True

    #---------------------------------------------------------------------------
    def create_trackFile(self, file):
        '''This function creates utput list but also a file after each step as a marker that the step is completed'''
        with open(file, 'w'):
            pass

    #---------------------------------------------------------------------------
    def run_command(self, command, step):
        '''This function executes a command and checks if it was executes without errors'''

        return_code = subprocess.run(command, shell=True)
        if return_code.returncode == 0:
            print(f"{step} without errors!, continuing with next step...\n")
        else:
            print('\nAn error has occured, see shell for more information. Exiting program...')
            sys.exit()

    #---------------------------------------------------------------------------
    def validate_id(self,options):
        '''This function validates wether the tumor_id and normal_id you entered is present in your reads'''

        if options.tumor_id and options.normal_id not in "".join(listdir(reads_dir)):
            print("You have entered a tumor_id and/or normal_id that is not present in your reads.\n")
            print("Please restart program and verify that you have typed in the right \"clinical_id\" for both tumor (-t) and normal (-n)! ")
            input("Press any key to exit program")
            sys.exit()
        else:
            print("tumor_id and normal_id correctly validated!\n")

    #---------------------------------------------------------------------------
    def index_reference(self):
        '''This function indexes the reference genome so it can be used in the analysis'''

        if self.step_completed(index_complete, 'Indexing reference genome'):
            time.sleep(2)
            return True
        else:
            print("1. Index reference genome\n")
            bwa_index = f"bwa index {reference_file}"
            self.run_command(bwa_index, 'Bwa index')
            create_dict = f"samtools dict {reference_file} -o {reference_file[:-2]}dict"
            self.run_command(create_dict, 'Creating .dict with samtools dict')
            create_fai = f"samtools faidx {reference_file}"
            self.run_command(create_fai, 'Creating .fai with samtools faidx')
            with open(index_complete, 'a'):
                pass
            print("\nIndexing reference genome completed!\n")
            input("Press any key to return to DNA-analysis menu...")

    #---------------------------------------------------------------------------
    def create_outputList(self, output_path, write_to_file):
        '''This function creates a list file containing the output name of the files created in the pipeline step where the function is used'''

        with open(output_path, 'a') as c:
            c.write(f'{write_to_file}\n')

    #---------------------------------------------------------------------------
    def validate_bam(self):
        '''This function runs picard ValidateSamFile to check if any errors are present in the aligned files.
           returns True if no errors are found or False if errors are found'''

        if self.step_completed(validateBam_complete, 'Validating files'):
            return True
        else:
            print("Validating .bam files...\n")
            with open(alignedFiles_list, 'r') as list:
                for sample in list.readlines():
                    cmd_validate = f"picard ValidateSamFile -I {output_aligned}{sample.rstrip()} -MODE SUMMARY"
                    self.run_command(cmd_validate, 'Picard ValidateSamFile')
                    self.create_trackFile(validateBam_complete)
                return True

    #---------------------------------------------------------------------------
    def alignment(self):
        '''This function loops through a library list containing either single-end or paired-end protocol, it will automatically detect what protocol it is.
           For every loop, the bwa mem command will be run and the output name for each run will be saved to a txt-file for use in the next step if the command finnish without errors.'''

        if self.step_completed(alignedFiles_list, 'Burrown Wheeler aligner'):
            pass
        else:
            threads = multiprocessing.cpu_count() - 2
            print("3. Run analysis\n\n\n")
            print("Burrows Wheeler aligner\n")
            print(f"Aligning reads using {threads} CPU threads out of {available_threads}...\n")
            self.create_directory(output_aligned)
            with open(f'{dna_seq_dir}library.txt', 'r') as fastq_list:
                for line in fastq_list.readlines():
                    clinical_id, library_id, read1, read2 = line.split()
                    read_group_header = f'\'@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}\''
                    if read2 == 'N/A': # single-end
                        cmd_bwa = f"1: bwa mem -R {read_group_header} {reference_file} {reads_dir}/{read1} -t {threads} | samtools view -bS -o {output_aligned}{library_id}.bam" # samtools view converts SAM to BAM
                    else: # paired-end
                        cmd_bwa = f"bwa mem -R {read_group_header} {reference_file} {reads_dir}/{read1} {reads_dir}/{read2} -t {threads} | samtools view -bS -o {output_aligned}{library_id}.bam"
                    # run_command(cmd_bwa, 'Alignment with Burrows Wheeler aligner')
                    self.create_outputList(alignedFiles_list, f"{library_id}.bam")
            return

    #---------------------------------------------------------------------------
    def sort(self, options):
        '''This function reads the completed_steps.txt to check if the previous step was completed without errors.
           The function returns a list of string containg the filenames of the sorted tumor samples and normal samples separated.
           The string is needed because you have several inputs in the next function and can therefore not run a for loop'''

        if self.step_completed(sortedFiles_list, 'Picard sortsam'):
            pass
        else:
            print("Sorting SAM/BAM files using Picard Sortsam...\n")

            # Empty strings to store the output
            tumor_sort_str = ""
            normal_sort_str = ""
            write_to_file = ""

            self.create_directory(output_sorted)
            with open(alignedFiles_list, 'r') as list:
                for sample in list.readlines():
                    cmd_sortsam = f"picard SortSam -I dna_seq/Aligned/{sample.rstrip()} -O dna_seq/Sorted/{sample.rstrip()} -SORT_ORDER coordinate --TMP_DIR $PWD"
                    self.run_command(cmd_sortsam, 'Sorting files with Picard SortSam')
                    if options.tumor_id in sample:
                        tumor_sort_str += f" -I dna_seq/Sorted/{sample}".rstrip()
                    else:
                        normal_sort_str += f" -I dna_seq/Sorted/{sample}".rstrip()
                write_to_file = tumor_sort_str.lstrip() + '\n' + normal_sort_str.lstrip()
                print(write_to_file)
                self.create_outputList(sortedFiles_list, write_to_file)
                print("\nPicard SortSam completed!\n")

    #---------------------------------------------------------------------------
    def merge(self, options):
        '''This function merges all the input files in the sortedFiles_list to one output file'''

        if self.step_completed(mergedFiles_list, 'Picard MergeSamFiles'):
            pass
        else:
            print("Merging SAM/BAM files using Picard MergeSamFiles...\n")
            self.create_directory(output_merged)
            with open(sortedFiles_list, 'r') as list:
                for sample in list.readlines():
                    if options.tumor_id in sample:
                        cmd_merge = f"picard MergeSamFiles {sample.rstrip()} -O dna_seq/Merged/{options.tumor_id}.bam"
                        self.run_command(cmd_merge, f"Merging {sample.rstrip()} with Picard MergeSamFiles")
                        self.create_outputList(mergedFiles_lis, f"{options.tumor_id}.bam")
                    else:
                        cmd_merge = f"picard MergeSamFiles {sample.rstrip()} -O dna_seq/Merged/{options.normal_id}.bam"
                        self.run_command(cmd_merge, f"Merging {sample.rstrip()} with Picard MergeSamFiles")
                        self.create_outputList(mergedFiles_list, f"{options.normal_id}.bam")
                print("\nPicard MergeSamFiles completed!\n")

    #---------------------------------------------------------------------------
    def remove_duplicate(self):
        '''This function removes duplicates '''

        if self.step_completed(removeDuplicates_list, 'Picard MarkDuplicates'):
            pass
        else:
            print("Removing duplicates in SAM/BAM files using Picard MarkDuplicates...\n")
            self.create_directory(output_removedDuplicates)
            with open(mergedFiles_list, 'r') as list:
                for sample in list.readlines():
                    cmd_rd = f"picard MarkDuplicates -I {output_merged}{sample.rstrip()} -O {output_removedDuplicates}{sample.rstrip()} -M {output_removedDuplicates}marked_dup_metrics_{sample.rstrip()}.txt"
                    print(cmd_rd)
                    self.run_command(cmd_rd, f'Removing duplicates in {sample.rstrip()} with Picard MarkDuplicates')
                    self.create_outputList(removeDuplicates_list, f"{sample.rstrip()}")
                print("\nPicard MarkDuplicates completed!\n")

    #---------------------------------------------------------------------------
    def realign(self):
        '''This function realigns the bam files'''

        if self.step_completed(realignedFiles_list, 'GATK LeftAlignIndels'):
            pass
        else:
            print("Realigning SAM/BAM files using GATK LeftAlignIndels...\n")
            self.create_directory(output_realigned)
            with open(removeDuplicates_list, 'r') as list:
                for sample in list.readlines():
                    cmd_index = f"samtools index {output_removedDuplicates}{sample.rstrip()}"
                    self.run_command(cmd_index, f'Indexing {sample.rstrip()} with samtools index')
                    cmd_leftAlignIndels = f"gatk LeftAlignIndels -R {reference_file} -I {output_removedDuplicates}{sample.rstrip()} -O {output_realigned}{sample.rstrip()}"
                    self.run_command(cmd_leftAlignIndels, f'Realigning {sample.rstrip()} with GATK LeftAlignIndels')
                    self.create_outputList(realignedFiles_list, f"{sample.rstrip()}")
                print("\nGATK LeftAlignIndels completed!\n")

    #---------------------------------------------------------------------------
    def gatk_snv(self, options):

        if self.step_completed(haplotypeCaller_complete, 'GATK haplotypeCaller'):
            pass
        else:
            print("Looking for SNV's using GATK HaplotypeCaller...\n")
            self.create_directory(output_haplotypecaller)
            with open(realignedFiles_list, 'r') as list:
                sample_1, sample_2 = list.readlines()
                if options.intervals:
                    cmd_call = f"gatk HaplotypeCaller -R {reference_file} -I {output_realigned}{sample_1.rstrip()} -I {output_realigned}{sample_2.rstrip()} -O {output_haplotypecaller}{options.tumor_id}.vcf -L {options.intervals}"
                else:
                    cmd_call = f"gatk HaplotypeCaller -R {reference_file} -I {output_realigned}{sample_1.rstrip()} -I {output_realigned}{sample_2.rstrip()} -O {output_haplotypecaller}{options.tumor_id}.vcf"
                self.run_command(cmd_call, f'Looking for SNV\'s with GATK haplotypeCaller')
                # Remove all reads with read depth less than 10
                cmd_filter_read_depth = f"bcftools view -i 'MIN(FMT/DP)>10' -e 'GT=1/1' -e 'GT=1/2' dna_seq/GATK_haplotypecaller/{options.tumor_id}.vcf > dna_seq/GATK_haplotypecaller/{options.tumor_id}_filtered_ReadDepthOver10.vcf"
                self.run_command(cmd_filter_read_depth, '')
                with open(haplotypeCaller_complete, 'a'):
                    pass
                print("\nGATK HaplotypeCaller completed!\n")

    #---------------------------------------------------------------------------
    def delly(self, options):
        '''This function creates an output directory and runs delly to call for somatic SNV's'''

        if self.step_completed(delly_complete, 'Delly'):
            pass
        else:
            print("Looking for somatic SNV's using delly\n")
            self.create_directory(output_delly)
            with open(realignedFiles_list, 'r') as list:
                sample_1, sample_2 = list.readlines()
                cmd_delly_call = f"delly call -x {file_exclude_template} -g {reference_file} -o {output_delly}delly.bcf {output_realigned}{sample_1.rstrip()} {output_realigned}{sample_2.rstrip()}"
                self.run_command(cmd_delly_call, f'Looking for somatic SNV\'s in {sample_1.rstrip()} and {sample_2.rstrip()} using delly')
                with open(f'{output_delly}sample.tsv', 'w', newline='') as tsv:
                    tsv_output = csv.writer(tsv, delimiter='\t')
                    tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
                    tsv_output.writerow([f"{options.normal_id}", 'control'])
                cmd_dos2unix = f"dos2unix {output_delly}sample.tsv"
                self.run_command(cmd_dos2unix, '')
                cmd_filter = f"delly filter -f somatic -o {output_delly}delly_filter.bcf -s {output_delly}sample.tsv {output_delly}delly.bcf"
                self.run_command(cmd_filter, 'Filtering variants...')
                cmd_convert = f"bcftools view {output_delly}delly_filter.bcf > {output_delly}delly_filter.vcf"
                self.run_command(cmd_convert, '')
                with open(delly_complete, 'a'):
                    pass
                print("\nDelly variant caller completed!\n")

    #---------------------------------------------------------------------------
    def manta(self):
        '''This function creates an output directory and runs manta to call for somatic SNV's'''

        if self.step_completed(manta_complete, 'Manta'):
            pass
        else:
            print("\nLooking for somatic SNV's using manta\n")
            self.create_directory(output_manta)
            with open(realignedFiles_list, 'r') as list:
                sample_1, sample_2 = list.readlines()
                cmd_create_config_file = f"{configManta_path} --tumorBam={output_realigned}{sample_1.rstrip()} --bam={output_realigned}{sample_2.rstrip()} --referenceFasta={reference_file} --runDir={output_manta}"
                self.run_command(cmd_create_config_file, '')
                cmd_runWorkflow = f"{runWorkflow_path} -m local -j 4"
                self.run_command(cmd_runWorkflow, '')
                cmd_unzip = f'gunzip {output_mantaVariants}somaticSV.vcf.gz'
                self.run_command(cmd_unzip, '')
                cmd_filter = f'bcftools view -i \'FILTER=="PASS"\' {output_mantaVariants}somaticSV.vcf > {output_mantaVariants}somaticSV_PASS.vcf'
                self.run_command(cmd_filter, 'Filtering of passed SNV\'s')
            with open(manta_complete, 'a'):
                pass
            print("\nManta variant caller completed!\n")
