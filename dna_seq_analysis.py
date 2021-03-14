# Packages used in script
from os import listdir, sys, mkdir, getenv, path
import subprocess
import argparse
import csv
import multiprocessing # lets you find out how many threads your cpu has
import time



class DnaSeqAnalysis():

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def create_outputList_dna(self, output_path, write_to_file):
        '''This function creates a list file containing the output name of the files created in the pipeline step where the function is used'''

        try:
            with open(output_path, 'a') as c:
                c.write(f'{write_to_file}\n')
        except Exception as e:
            print(f'Error with create_outputList_dna(): {e}')

    #---------------------------------------------------------------------------
    def validate_bam_dna(self, misc, shortcuts):
        '''This function runs picard ValidateSamFile to check if any errors are present in the aligned files.
           returns True if no errors are found or False if errors are found'''

        try:
            if misc.step_completed(shortcuts.validate_bam_complete, 'Validating files allready completed, skips step...'):
                return True
            else:
                print("Validating .bam files...\n")
                with open(shortcuts.alignedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_validate = f"picard ValidateSamFile -I {shortcuts.aligned_output_dir}{sample.rstrip()} -MODE SUMMARY"
                        misc.run_command(cmd_validate, 'Picard ValidateSamFile')
                        misc.create_trackFile(shortcuts.validate_bam_complete)
                    return True
        except Exception as e:
            print(f'Error with validate_bam_dna(): {e}')

    #---------------------------------------------------------------------------
    def alignment(self, misc, shortcuts):
        '''This function loops through a library list containing either single-end or paired-end protocol, it will automatically detect what protocol it is.
           For every loop, the bwa mem command will be run and the output name for each run will be saved to a txt-file for use in the next step if the command finnish without errors.'''

        try:
            if misc.step_completed(shortcuts.alignedFiles_list, 'Burrown Wheeler aligner allready completed, skips step...'):
                pass
            else:
                threads = multiprocessing.cpu_count() - 2
                print("3. Run analysis\n\n\n")
                print("Burrows Wheeler aligner\n")
                print(f"Aligning reads using {threads} CPU threads out of {available_threads}...\n")
                with open(f'{shortcuts.dna_seq_dir}library.txt', 'r') as fastq_list:
                    for line in fastq_list.readlines():
                        clinical_id, library_id, read1, read2 = line.split()
                        read_group_header = f'\'@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}\''
                        if read2 == 'N/A': # single-end
                            cmd_bwa = f"1: bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}/{read1} -t {threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam" # samtools view converts SAM to BAM
                        else: # paired-end
                            cmd_bwa = f"bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}/{read1} {shortcuts.dna_reads_dir}/{read2} -t {threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam"
                        # run_command(cmd_bwa, 'Alignment with Burrows Wheeler aligner')
                        self.create_outputList(shortcuts.alignedFiles_list, f"{library_id}.bam")
                return
        except Exception as e:
            print(f'Error with alignment(): {e}')

    #---------------------------------------------------------------------------
    def sort(self, options, misc, shortcuts):
        '''This function reads the completed_steps.txt to check if the previous step was completed without errors.
           The function returns a list of string containg the filenames of the sorted tumor samples and normal samples separated.
           The string is needed because you have several inputs in the next function and can therefore not run a for loop'''

        try:
            if misc.step_completed(shortcuts.sortedFiles_list, 'Picard sortsam allready completed, skips step...'):
                pass
            else:
                print("Sorting SAM/BAM files using Picard Sortsam...\n")

                # Empty strings to store the output
                tumor_sort_str = ""
                normal_sort_str = ""
                write_to_file = ""

                with open(shortcuts.alignedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_sortsam = f"picard SortSam -I {shortcuts.aligned_output_dir}{sample.rstrip()} -O {shortcuts.sorted_output_dir}{sample.rstrip()} -SORT_ORDER coordinate --TMP_DIR $PWD"
                        misc.run_command(cmd_sortsam, 'Picard SortSam completed')
                        if options.tumor_id in sample:
                            tumor_sort_str += f" -I {shortcuts.sorted_output_dir}{sample}".rstrip()
                        else:
                            normal_sort_str += f" -I {shortcuts.sorted_output_dir}{sample}".rstrip()
                    write_to_file = tumor_sort_str.lstrip() + '\n' + normal_sort_str.lstrip()
                    print(write_to_file)
                    self.create_outputList_dna(shortcuts.sortedFiles_list, write_to_file)
                    print("\nPicard SortSam completed!\n")
        except Exception as e:
            print(f'Error with sort(): {e}')
    #---------------------------------------------------------------------------
    def merge(self, options, misc, shortcuts):
        '''This function merges all the input files in the sortedFiles_list to one output file'''

        try:
            if misc.step_completed(shortcuts.mergedFiles_list, 'Picard MergeSamFiles allready completed, skips step...'):
                pass
            else:
                print("Merging SAM/BAM files using Picard MergeSamFiles...\n")
                with open(shortcuts.sortedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        if options.tumor_id in sample:
                            cmd_merge = f"picard MergeSamFiles {sample.rstrip()} -O {shortcuts.merged_output_dir}{options.tumor_id}.bam"
                            misc.run_command(cmd_merge, f"Merging {sample.rstrip()} with Picard MergeSamFiles completed")
                            self.create_outputList_dna(shortcuts.mergedFiles_list, f"{options.tumor_id}.bam")
                        else:
                            cmd_merge = f"picard MergeSamFiles {sample.rstrip()} -O {shortcuts.merged_output_dir}{options.normal_id}.bam"
                            misc.run_command(cmd_merge, f"Merging {sample.rstrip()} with Picard MergeSamFiles completed")
                            self.create_outputList_dna(shortcuts.mergedFiles_list, f"{options.normal_id}.bam")
                    print("\nPicard MergeSamFiles completed!\n")
        except Exception as e:
            print(f'Error with merge(): {e}')

    #---------------------------------------------------------------------------
    def remove_duplicate(self, misc, shortcuts):
        '''This function removes duplicates '''

        try:
            if misc.step_completed(shortcuts.removeDuplicates_list, 'Picard MarkDuplicates allready completed, skips step...'):
                pass
            else:
                print("Removing duplicates in SAM/BAM files using Picard MarkDuplicates...\n")
                with open(shortcuts.mergedFiles_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_rd = f"picard MarkDuplicates -I {shortcuts.removed_duplicates_output_dir}{sample.rstrip()} -O {shortcuts.removed_duplicates_output_dir}{sample.rstrip()} -M {shortcuts.removed_duplicates_output_dir}marked_dup_metrics_{sample.rstrip()}.txt"
                        misc.run_command(cmd_rd, f'Removing duplicates in {sample.rstrip()} with Picard MarkDuplicates completed')
                        self.create_outputList_dna(shortcuts.removeDuplicates_list, f"{sample.rstrip()}")
                    print("\nPicard MarkDuplicates completed!\n")
        except Exception as e:
            print(f'Error with remove_duplicate(): {e}')

    #---------------------------------------------------------------------------
    def realign(self, misc, shortcuts):
        '''This function realigns the bam files'''

        try:
            if misc.step_completed(shortcuts.realignedFiles_list, 'GATK LeftAlignIndels allready completed, skips step...'):
                pass
            else:
                print("Realigning SAM/BAM files using GATK LeftAlignIndels...\n")
                with open(shortcuts.removeDuplicates_list, 'r') as list:
                    for sample in list.readlines():
                        cmd_index = f"samtools index {shortcuts.realigned_output_dir}{sample.rstrip()}"
                        misc.run_command(cmd_index, f'Indexing {sample.rstrip()} with samtools index step_completed')
                        cmd_leftAlignIndels = f"gatk LeftAlignIndels -R {shortcuts.reference_genome_file} -I {shortcuts.removed_duplicates_output_dir}{sample.rstrip()} -O {shortcuts.realigned_output_dir}{sample.rstrip()}"
                        misc.run_command(cmd_leftAlignIndels, f'Realigning {sample.rstrip()} with GATK LeftAlignIndels completed')
                        self.create_outputList_dna(shortcuts.realignedFiles_list, f"{sample.rstrip()}")
                    print("\nGATK LeftAlignIndels completed!\n")
        except Exception as e:
            print(f'Error with realign(): {e}')

    #---------------------------------------------------------------------------
    def gatk_haplotype(self, options, misc, shortcuts):

        try:
            if misc.step_completed(shortcuts.haplotypecaller_complete, 'GATK haplotypeCaller allready completed, skips step...'):
                pass
            else:
                print("Looking for SNV's using GATK HaplotypeCaller...\n")
                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.readlines()
                    if options.intervals:
                        cmd_call = f"gatk HaplotypeCaller -R {shortcuts.reference_genome_file} -I {shortcuts.realigned_output_dir}{sample_1.rstrip()} -I {shortcuts.realigned_output_dir}{sample_2.rstrip()} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf -L {options.intervals}"
                    else:
                        cmd_call = f"gatk HaplotypeCaller -R {shortcuts.reference_genome_file} -I {shortcuts.realigned_output_dir}{sample_1.rstrip()} -I {shortcuts.realigned_output_dir}{sample_2.rstrip()} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf"
                    # misc.run_command(cmd_call, f'Looking for SNV\'s with GATK haplotypeCaller')
                    # Remove all reads with read depth less than 10, selects only snps, exludes normal samples
                    cmd_filter_read_depth = f"bcftools view -i 'MIN(FMT/DP)>10' -v snps -s ^987-02 {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf"
                    # misc.run_command(cmd_filter_read_depth, '')
                    # select heterozygous genotype, excludes GT=1/2
                    cmd_filter_het = f"bcftools view -g het -e 'GT=\"1/2\"' {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf"
                    # misc.run_command(cmd_filter_het, None)
                    cmd_indexFeatureFile = f"gatk IndexFeatureFile -I {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf"
                    # misc.run_command(cmd_indexFeatureFile, None)

                    # Annotates vcf file
                    cmd_annotate = f'''java -Xmx4g -jar $HOME/anaconda3/envs/sequencing/share/snpeff-5.0-0/snpEff.jar \\
                    -v GRCh38.99 -canon -noInteraction -noNextProt -noMotif -strict \\
                    -onlyProtein {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf \\
                    > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf'''
                    misc.run_command(cmd_annotate, 'Annotation of vcf file')
                    misc.create_trackFile(shortcuts.haplotypecaller_complete)
                    print("\nGATK HaplotypeCaller completed!\n")
        except Exception as e:
            print(f'Error with gatk_haplotype: {e}')

    #---------------------------------------------------------------------------
    def delly(self, options, misc, shortcuts):
        '''This function creates an output directory and runs delly to call for somatic SNV's'''

        try:
            if misc.step_completed(shortcuts.delly_complete, 'Delly allready completed, skips step...'):
                pass
            else:
                print("Looking for somatic SNV's using delly\n")
                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.readlines()
                    cmd_delly_call = f"delly call -x {shortcuts.reference_genome_exclude_template_file} -g {shortcuts.reference_genome_file} -o {shortcuts.delly_output_dir}delly.bcf {shortcuts.realigned_output_dir}{sample_1.rstrip()} {shortcuts.realigned_output_dir}{sample_2.rstrip()}"
                    misc.run_command(cmd_delly_call, f'Calling somatic SNV\'s in {sample_1.rstrip()} and {sample_2.rstrip()} using delly completed')
                    with open(f'{shortcuts.delly_output_dir}sample.tsv', 'w', newline='') as tsv:
                        tsv_output = csv.writer(tsv, delimiter='\t')
                        tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
                        tsv_output.writerow([f"{options.normal_id}", 'control'])
                    cmd_dos2unix = f"dos2unix {shortcuts.delly_output_dir}sample.tsv"
                    misc.run_command(cmd_dos2unix, '')
                    cmd_filter = f"delly filter -f somatic -o {shortcuts.delly_output_dir}delly_filter.bcf -s {shortcuts.delly_output_dir}sample.tsv {shortcuts.delly_output_dir}delly.bcf"
                    misc.run_command(cmd_filter, 'Filtering variants completed')
                    cmd_convert = f"bcftools view {shortcuts.delly_output_dir}delly_filter.bcf > {shortcuts.delly_output_dir}delly_filter.vcf"
                    misc.run_command(cmd_convert, 'Converting .bcf to .vcf completed')
                    misc.create_trackFile(shortcuts.delly_complete)
                    print("\nDelly variant caller completed!\n")
        except Exception as e:
            print(f'Error with delly(): {e}')

        #---------------------------------------------------------------------------
    def manta(self, misc, shortcuts):
        '''This function creates an output directory and runs manta to call for somatic SNV's'''

        try:
            if misc.step_completed(shortcuts.manta_complete, 'Manta allready completed, skips step...'):
                pass
            else:
                print("\nLooking for somatic SNV's using manta\n")
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
                print("\nManta variant caller completed!\n")
        except Exception as e:
            print(f'Error with manta(): {e}')
