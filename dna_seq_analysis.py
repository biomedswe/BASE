# Packages used in script
from os import listdir, sys, mkdir, getenv, path
import subprocess
import argparse
import csv
import multiprocessing
import time
import timeit
import logging
logging.basicConfig(filename='Logname.txt', filemode='a', format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)



class DnaSeqAnalysis():

    def __init__(self):
        pass


    #---------------------------------------------------------------------------
    def index_genome_dna(self, misc, shortcuts):
        '''This function indexes the reference genome so it can be used in the analysis'''

        try:
            ref_file = shortcuts.reference_genome_file
            ref_dir = shortcuts.reference_genome_dir
            allready_completed = shortcuts.bwa_index_whole_reference_genome_complete

            if misc.step_allready_completed(allready_completed):
                misc.log_to_file('Burrows Wheeler aligner index allready completed, skips step...')
                time.sleep(2.5)
            else:
                misc.log_to_file('\nStarting: indexing with bwa index')
                misc.clear_screen()
                print("\033[1mIndex whole genome\033[0m\n")
                cmd_bwa_index = f"bwa index {ref_file}"
                misc.run_command(cmd_bwa_index)
                misc.log_to_file('Bwa index completed without errors')
                cmd_create_dict = f"samtools dict {ref_file} -o {ref_file[:-2]}dict"
                misc.run_command(cmd_create_dict)
                misc.log_to_file('Creating .dict with samtools dict completed without errors')
                cmd_create_fai = f"samtools faidx {ref_file} -o {ref_file}.fai"
                misc.run_command(cmd_create_fai)
                misc.log_to_file('Creating .fai with samtools faidx completed without errors')
                misc.create_trackFile(allready_completed)
                misc.log_to_file('Indexing reference genome successfully completed!\n')

        except Exception as e:
            misc.log_to_file(f'Error with index_genome_dna() in dna_seq_analysis.py: {e}')
            input("Press any key to continue")

    #---------------------------------------------------------------------------
    def create_outputList_dna(self, output_path, write_to_file):
        '''This function creates a list file containing the output name of the files created in the pipeline step where the function is used, the list is then used in the next step'''

        try:
            with open(output_path, 'a') as c:
                c.write(f'{write_to_file}\n')
        except Exception as e:
            logging.info(f'Error with create_outputList_dna() in dna_seq_analysis.py: {e}')

    #---------------------------------------------------------------------------
    def validate_bam_dna(self, misc, shortcuts):
        '''This function runs picard ValidateSamFile to check if any errors are present in the aligned files.
           returns True if no errors are found or False if errors are found'''

        try:
            if misc.step_allready_completed(shortcuts.validate_bam_complete):
                logging.info('Validate bam allready completed, skips step...')

                return True
            else:
                print("Validating .bam files...\n")
                with open(shortcuts.alignedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_validate = f"picard ValidateSamFile -I {shortcuts.aligned_output_dir}{sample.rstrip()} -MODE SUMMARY"
                        misc.run_command(cmd_validate)
                        misc.create_trackFile(shortcuts.validate_bam_complete)
                    misc.log_to_file('Picard ValidateSamFile completed - OK!')
                    return True
        except Exception as e:
            misc.log_to_file(f'Error with validate_bam_dna() in dna_seq_analysis.py: {e}')

    #---------------------------------------------------------------------------
    def alignment(self, misc, shortcuts):
        '''This function loops through a library list containing either single-end or paired-end protocol, it will automatically detect what protocol it is.
           For every loop, the bwa mem command will be run and the output name for each run will be saved to a txt-file for use in the next step if the command finnish without errors.'''

        # Whole genome alignment
        try:
            if misc.step_allready_completed(shortcuts.alignedFiles_list):
                misc.log_to_file('Burrown Wheeler aligner allready completed, skips step...')
            else:
                start = timeit.default_timer()
                threads = multiprocessing.cpu_count() - 2
                misc.log_to_file(f'\nStarting: Burrows Wheeler aligner\nUsing {threads} out of {threads+2} available threads')

                with open(f'{shortcuts.dna_seq_dir}library.txt', 'r') as fastq_list:
                    for line in fastq_list.readlines():
                        clinical_id, library_id, read1, read2 = line.split()
                        read_group_header = f'\'@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}\''
                        if read2 == 'N/A': # single-end
                            cmd_bwa = f"1: bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}/{read1} -t {threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam" # samtools view converts SAM to BAM
                        else: # paired-end
                            cmd_bwa = f"bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}/{read1} {shortcuts.dna_reads_dir}/{read2} -t {threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam"
                        # misc.run_command(cmd_bwa)
                        self.create_outputList_dna(shortcuts.alignedFiles_list, f"{library_id}.bam")
                end = timeit.default_timer()
                elapsed_time = end-start
                misc.log_to_file(f'Burrows Wheeler aligner succesfully completed in {elapsed_time/60:.1g} min')
        except Exception as e:
            misc.log_to_file(f'Error with def alignment() in dna_seq_analysis.py: {e}')
            input('press any key to exit')

    #---------------------------------------------------------------------------
    def sort(self, options, misc, shortcuts):
        '''This function reads the completed_steps.txt to check if the previous step was completed without errors.
           The function returns a list of string containg the filenames of the sorted tumor samples and normal samples separated.
           The string is needed because you have several inputs in the next function and can therefore not run a for loop'''

        try:
            if misc.step_allready_completed(shortcuts.sortedFiles_list):
                misc.log_to_file('Picard sortsam allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: sorting SAM/BAM files using Picard Sortsam")


                # Empty strings to store the output
                tumor_sort_str = ""
                normal_sort_str = ""
                write_to_file = ""
                with open(shortcuts.alignedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_sortsam = f"picard SortSam -I {shortcuts.aligned_output_dir}{sample.rstrip()} -O {shortcuts.sorted_output_dir}{sample.rstrip()} -SORT_ORDER coordinate --TMP_DIR $PWD"
                        return_code = misc.run_command(cmd_sortsam)
                        if options.tumor_id in sample:
                            tumor_sort_str += f" -I {shortcuts.sorted_output_dir}{sample}".rstrip()
                        else:
                            normal_sort_str += f" -I {shortcuts.sorted_output_dir}{sample}".rstrip()
                    write_to_file = tumor_sort_str.lstrip() + '\n' + normal_sort_str.lstrip()
                    self.create_outputList_dna(shortcuts.sortedFiles_list, write_to_file)
                misc.log_to_file('Picard SortSam completed - OK!')
                if return_code == 0:
                    cmd_remove_aligned_files = f'rm {shortcuts.aligned_output_dir}*.bam'
                    misc.run_command(cmd_remove_aligned_files)
                    misc.log_to_file('Aligned BAM files removed')
        except Exception as e:
            misc.log_to_file(f'Error with sort() in dna_seq_analysis.py: {e}')
            input('press any key to exit')
    #---------------------------------------------------------------------------
    def merge(self, options, misc, shortcuts):
        '''This function merges all the input files in the sortedFiles_list to one output file'''

        try:
            if misc.step_allready_completed(shortcuts.mergedFiles_list):
                misc.log_to_file('Picard MergeSamFiles allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: merging SAM/BAM files using Picard MergeSamFiles")

                with open(shortcuts.sortedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        if options.tumor_id in sample:
                            cmd_merge = f"picard MergeSamFiles {sample.rstrip()} -O {shortcuts.merged_output_dir}{options.tumor_id}.bam"
                            return_code = misc.run_command(cmd_merge)
                            self.create_outputList_dna(shortcuts.mergedFiles_list, f"{options.tumor_id}.bam")
                        else:
                            cmd_merge = f"picard MergeSamFiles {sample.rstrip()} -O {shortcuts.merged_output_dir}{options.normal_id}.bam"
                            return_code = misc.run_command(cmd_merge)
                            self.create_outputList_dna(shortcuts.mergedFiles_list, f"{options.normal_id}.bam")
                misc.log_to_file('Picard MergeSamFiles completed - OK!')
                if return_code == 0:
                    cmd_remove_sorted_files = f'rm {shortcuts.sorted_output_dir}*.bam'
                    misc.run_command(cmd_remove_sorted_files)
                    misc.log_to_file('Sorted BAM files removed')
        except Exception as e:
            misc.log_to_file(f'Error with merge() in dna_seq_analysis.py: {e}')
            input('press any key to exit')

    #---------------------------------------------------------------------------
    def remove_duplicate(self, misc, shortcuts):
        '''This function removes duplicates '''

        try:
            if misc.step_allready_completed(shortcuts.removeDuplicates_list):
                misc.log_to_file('Picard MarkDuplicates allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: removing duplicates in SAM/BAM files using Picard MarkDuplicates")

                with open(shortcuts.mergedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_rd = f"picard MarkDuplicates -I {shortcuts.merged_output_dir}{sample.rstrip()} -O {shortcuts.removed_duplicates_output_dir}{sample.rstrip()} -M {shortcuts.removed_duplicates_output_dir}marked_dup_metrics_{sample.rstrip()}.txt"
                        return_code = misc.run_command(cmd_rd)
                        self.create_outputList_dna(shortcuts.removeDuplicates_list, f"{sample.rstrip()}")
                misc.log_to_file('Picard MarkDuplicates completed - OK!')
                if return_code == 0:
                    cmd_remove_merged_files = f'rm {shortcuts.merged_output_dir}*.bam'
                    misc.run_command(cmd_remove_merged_files)
                    misc.log_to_file('Merged BAM files removed')
        except Exception as e:
            misc.log_to_file(f'Error with remove_duplicate() in dna_seq_analysis.py: {e}')
            input('press any key to exit')

    #---------------------------------------------------------------------------
    def realign(self, misc, shortcuts):
        '''This function realigns the bam files'''

        try:
            if misc.step_allready_completed(shortcuts.realignedFiles_list):
                misc.log_to_file('GATK LeftAlignIndels allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: realigning SAM/BAM files using GATK LeftAlignIndels")

                with open(shortcuts.removeDuplicates_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_index = f"samtools index {shortcuts.removed_duplicates_output_dir}{sample.rstrip()}"
                        misc.run_command(cmd_index)
                        cmd_leftAlignIndels = f"gatk LeftAlignIndels -R {shortcuts.reference_genome_file} -I {shortcuts.removed_duplicates_output_dir}{sample.rstrip()} -O {shortcuts.realigned_output_dir}{sample.rstrip()}"
                        return_code = misc.run_command(cmd_leftAlignIndels)
                        self.create_outputList_dna(shortcuts.realignedFiles_list, f"{sample.rstrip()}")
                misc.log_to_file('gatk LeftAlignIndels completed - OK!')
        except Exception as e:
            misc.log_to_file(f'Error with realign() in dna_seq_analysis.py: {e}')
            input('press any key to exit')

    #---------------------------------------------------------------------------
    def gatk_haplotype(self, options, misc, shortcuts):

        try:
            if misc.step_allready_completed(shortcuts.haplotypecaller_complete):
                misc.log_to_file('GATK haplotypeCaller allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: looking for SNV's using GATK HaplotypeCaller")

                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.readlines()
                    if options.intervals:
                        cmd_call = f"gatk HaplotypeCaller -R {shortcuts.reference_genome_file} -I {shortcuts.realigned_output_dir}{sample_1.rstrip()} -I {shortcuts.realigned_output_dir}{sample_2.rstrip()} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf -L {options.intervals}"
                    else:
                        cmd_call = f"gatk HaplotypeCaller -R {shortcuts.reference_genome_file} -I {shortcuts.realigned_output_dir}{sample_1.rstrip()} -I {shortcuts.realigned_output_dir}{sample_2.rstrip()} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf"
                    misc.run_command(cmd_call)
                    # Remove all reads with read depth less than 10, selects only snps, exludes normal samples
                    cmd_filter_read_depth = f"bcftools view -i 'MIN(FMT/DP)>10' -v snps -s ^987-02 {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf"
                    misc.run_command(cmd_filter_read_depth)
                    # select heterozygous genotype, excludes GT=1/2
                    cmd_filter_het = f"bcftools view -g het -e 'GT=\"1/2\"' {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf"
                    misc.run_command(cmd_filter_het)
                    cmd_indexFeatureFile = f"gatk IndexFeatureFile -I {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf"
                    misc.run_command(cmd_indexFeatureFile)

                    # Annotates vcf file
                    cmd_annotate = f'''java -Xmx4g -jar $HOME/anaconda3/envs/sequencing/share/snpeff-5.0-0/snpEff.jar \\
                    -v GRCh38.99 -canon -noInteraction -noNextProt -noMotif -strict \\
                    -onlyProtein {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf \\
                    > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf'''
                    misc.run_command(cmd_annotate)
                    misc.create_trackFile(shortcuts.haplotypecaller_complete)
                misc.log_to_file('gatk HaplotypeCaller succesfully completed, vcf file is filtered and annotated!')
        except Exception as e:
            misc.log_to_file(f'Error with gatk_haplotype in dna_seq_analysis.py: {e}')
            input('press any key to exit')

    #---------------------------------------------------------------------------
    def delly(self, options, misc, shortcuts):
        '''This function creates an output directory and runs delly to call for somatic SNV's'''

        try:
            if misc.step_allready_completed(shortcuts.delly_complete):
                misc.log_to_file('Delly allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: looking for somatic SNV's using delly")

                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.readlines()
                    cmd_delly_call = f"delly call -x {shortcuts.reference_genome_exclude_template_file} -g {shortcuts.reference_genome_file} -o {shortcuts.delly_output_dir}delly.bcf {shortcuts.realigned_output_dir}{sample_1.rstrip()} {shortcuts.realigned_output_dir}{sample_2.rstrip()}"
                    misc.run_command(cmd_delly_call)
                with open(f'{shortcuts.delly_output_dir}sample.tsv', 'w', newline='') as tsv:
                    tsv_output = csv.writer(tsv, delimiter='\t')
                    tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
                    tsv_output.writerow([f"{options.normal_id}", 'control'])
                cmd_dos2unix = f"dos2unix {shortcuts.delly_output_dir}sample.tsv"
                misc.run_command(cmd_dos2unix)
                cmd_filter = f"delly filter -f somatic -o {shortcuts.delly_output_dir}delly_filter.bcf -s {shortcuts.delly_output_dir}sample.tsv {shortcuts.delly_output_dir}delly.bcf"
                misc.run_command(cmd_filter)
                cmd_convert = f"bcftools view {shortcuts.delly_output_dir}delly_filter.bcf > {shortcuts.delly_output_dir}delly_filter.vcf"
                misc.run_command(cmd_convert)
                misc.create_trackFile(shortcuts.delly_complete)
                misc.log_to_file('Delly succesfully completed')
        except Exception as e:
            misc.log_to_file(f'Error with delly() in dna_seq_analysis.py: {e}')
            input('press any key to exit')

        #---------------------------------------------------------------------------
    def manta(self, misc, shortcuts):
        '''This function creates an output directory and runs manta to call for somatic SNV's'''

        try:
            if misc.step_allready_completed(shortcuts.manta_complete):
                misc.log_to_file('Manta allready completed, skips step...')
            else:
                misc.log_to_file("\nStarting: looking for somatic SNV's using manta")

                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.readlines()
                    cmd_create_config_file = f"{shortcuts.configManta_file} --tumorBam={shortcuts.realigned_output_dir}{sample_1.rstrip()} --bam={shortcuts.realigned_output_dir}{sample_2.rstrip()} --referenceFasta={shortcuts.reference_genome_file} --runDir={shortcuts.manta_output_dir}"
                    misc.run_command(cmd_create_config_file, '')
                    cmd_runWorkflow = f"{shortcuts.runWorkflow_file} -m local -j 4"
                    misc.run_command(cmd_runWorkflow, '')
                    cmd_unzip = f'gunzip {shortcuts.manta_variants_dir}somaticSV.vcf.gz'
                    misc.run_command(cmd_unzip, '')
                    cmd_filter = f'bcftools view -i \'FILTER=="PASS"\' {shortcuts.manta_variants_dir}somaticSV.vcf > {shortcuts.manta_variants_dir}somaticSV_PASS.vcf'
                    misc.run_command(cmd_filter, 'Filtering of passed SNV\'s')
                misc.create_trackFile(shortcuts.manta_complete)
                misc.log_to_file('Manta succesfully created')
        except Exception as e:
            misc.log_to_file(f'Error with manta() in dna_seq_analysis.py: {e} seconds')
            input('press any key to exit')
