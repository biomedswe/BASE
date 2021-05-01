# Packages used in script
from os import listdir, sys, mkdir, getenv, path, rename, remove
import subprocess
import argparse
import csv
import multiprocessing as mp
import time
import timeit
import re
from shutil import copy
from functools import partial




class DnaSeqAnalysis():

    def __init__(self):
        pass


    #---------------------------------------------------------------------------
    def index_genome_dna(self, misc, shortcuts):
        '''This function indexes the reference genome so it can be used in the analysis'''

        try:
            ref_file = shortcuts.reference_genome_file
            chunks_dir = shortcuts.reference_genome_chunks_dir
            allready_completed = shortcuts.bwa_index_complete

            start = timeit.default_timer()
            if not misc.step_allready_completed(allready_completed, "Indexing GRCh38.p13.genome"):
                misc.clear_screen()
                misc.log_to_file('Starting: indexing GRCh38.p13.genome with bwa index')
                cmd_bwa_index = f"bwa index {ref_file}"
                misc.run_command(cmd_bwa_index, "Bwa index", f"{ref_file}.sa", None)
                cmd_create_dict = f"samtools dict {ref_file} -o {ref_file[:-2]}dict"
                misc.run_command(cmd_create_dict, "Creating .dict with samtools dict", f"{ref_file[:-2]}dict", None)
                cmd_create_fai = f"samtools faidx {ref_file} -o {ref_file}.fai"
                misc.run_command(cmd_create_fai, "Creating .fai with samtools faidx", f"{ref_file}.fai", None)
            cmd_split_fasta = f"bedtools makewindows -w 10000000 -g {ref_file}.fai > {chunks_dir}chunk.bed"
            misc.run_command(cmd_split_fasta, "Spliting fa.fai with bedtools makewindows", f"{chunks_dir}chunk.bed", None)
            cmd_split_bed = f"split -l 2 {chunks_dir}chunk.bed {chunks_dir}chunk.split"
            misc.run_command(cmd_split_bed, "Splitting chunk.bed to one file per line", f"{chunks_dir}chunk.splitaa.bed", allready_completed)
            for file in listdir(chunks_dir):
                if re.search("split*", file) and not file.endswith("bed"): rename(f"{chunks_dir}{file}", f"{chunks_dir}{file}.bed")
                if file == "chunk.bed": remove(f"{chunks_dir}{file}")
            misc.log_to_file('Renaming chunk.split* > chunk.split*.bed succesfully completed - OK!')
            elapsed = timeit.default_timer() - start
            misc.log_to_file(f'Indexing reference genome successfully completed in {misc.elapsed_time(elapsed)} - OK!!')
            input("Press any key to return to DNA analysis menu...")
            return ''
        except Exception as e:
            misc.log_exception(".index_genome_dna() in dna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def validate_bam_dna(self, options, misc, shortcuts):
        '''This function runs picard ValidateSamFile to check if any errors are present in the aligned files.
           returns True if no errors are found or False if errors are found'''

        try:
            start = timeit.default_timer()
            misc.log_to_file("Validating .bam files...")
            with open(shortcuts.alignedFiles_list, 'r') as list:
                for sample in list.read().splitlines():
                    cmd_validate = f"picard ValidateSamFile -I {shortcuts.aligned_output_dir}{options.tumor_id}/{sample} -MODE SUMMARY"
                    validated = misc.run_command(cmd_validate, f"Picard ValidateSamFile {sample[:-4]}", f"{shortcuts.aligned_output_dir}{options.tumor_id}/{sample[:-6]}.validated", f"{shortcuts.aligned_output_dir}{options.tumor_id}/{sample[:-6]}.validated")
                elapsed = timeit.default_timer() - start
                if validated: misc.log_to_file(f'All .bam files succesfully validated in {misc.elapsed_time(elapsed)} - OK!')
            return True
        except Exception as e:
            misc.log_exception(".validate_bam_dna() in dna_seq_analysis.py:", e)
            sys.exit()
    #---------------------------------------------------------------------------
    def alignment(self, options, misc, shortcuts):
        '''This function loops through a library list containing either single-end or paired-end protocol, it will automatically detect what protocol it is.
           For every loop, the bwa mem command will be run and the output name for each run will be saved to a txt-file for use in the next step if the command finnish without errors.'''

        try:
            if not misc.step_allready_completed(shortcuts.alignedFiles_list, "Burrows Wheeler aligner"):
                misc.create_directory([f"{shortcuts.aligned_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                threads = mp.cpu_count() - 2
                misc.log_to_file(f'Starting: Burrows Wheeler aligner Using {threads} out of {threads+2} available threads')
                with open(f'{shortcuts.dna_seq_dir}{options.tumor_id}_library.txt', 'r') as fastq_list:
                    for line in fastq_list.readlines():
                        clinical_id, library_id, read1, read2 = line.split()
                        read_group_header = f'\'@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}\''
                        if read2 == 'N/A': # single-end
                            cmd_bwa = f"1: bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}/{read1} -t {threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{options.tumor_id}/{library_id}.bam" # samtools view converts SAM to BAM
                        else: # paired-end
                            cmd_bwa = f"bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}/{read1} {shortcuts.dna_reads_dir}/{read2} -t {threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{options.tumor_id}/{library_id}.bam"
                        misc.run_command(cmd_bwa, f"Aligning {library_id}.bam", f"{shortcuts.aligned_output_dir}{options.tumor_id}/{library_id}.bam", None)
                        misc.create_outputList_dna(shortcuts.alignedFiles_list, f"{library_id}.bam")
                    elapsed = timeit.default_timer() - start
                    misc.log_to_file(f'Burrows Wheeler aligner succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except OSError as e:
            misc.log_exception("You have to create a library file first", e)
            sys.exit()
        except Exception as e:
            misc.log_exception("alignment() in dna_seq_analysis.py:", e)
            sys.exit()


    #---------------------------------------------------------------------------
    def sort(self, options, misc, shortcuts):
        '''This function reads the completed_steps.txt to check if the previous step was completed without errors.
           The function returns a list of string containg the filenames of the sorted tumor samples and normal samples separated.
           The string is needed because you have several inputs in the next function and can therefore not run a for loop'''

        try:
            if not misc.step_allready_completed(shortcuts.sortedFiles_list, "Picard sortsam"):
                misc.create_directory([f"{shortcuts.sorted_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                misc.log_to_file("Starting: sorting SAM/BAM files using Picard Sortsam")
                # Empty strings to store the output
                cmd_sort = []
                tumor_sort_str = ""
                normal_sort_str = ""
                write_to_file = ""
                with open(shortcuts.alignedFiles_list, 'r') as list:
                    for sample in list.read().splitlines():
                        cmd_sort.append(f"picard SortSam -I {shortcuts.aligned_output_dir}{options.tumor_id}/{sample} -O {shortcuts.sorted_output_dir}{options.tumor_id}/{sample} -SORT_ORDER coordinate --TMP_DIR {shortcuts.sorted_output_dir}{options.tumor_id}/tmp")
                        if options.tumor_id in sample:
                            tumor_sort_str += f" -I {shortcuts.sorted_output_dir}{options.tumor_id}/{sample}".rstrip()
                        else:
                            normal_sort_str += f" -I {shortcuts.sorted_output_dir}{options.tumor_id}/{sample}".rstrip()
                    write_to_file = f"{tumor_sort_str.lstrip()}\n{normal_sort_str.lstrip()}"

                    with mp.Pool(processes=2) as pool:
                        pool.map(partial(self.multi_processing, f"{shortcuts.sorted_output_dir}{options.tumor_id}/{sample}", options, misc, shortcuts),cmd_sort)
                    misc.create_outputList_dna(shortcuts.sortedFiles_list, write_to_file)
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'Picard SortSam succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                misc.run_command(f"rm {shortcuts.aligned_output_dir}{options.tumor_id}/*.bam", 'Removing aligned BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".sort() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def merge(self, options, misc, shortcuts):
        '''This function merges all the input files in the sortedFiles_list to one output file'''

        try:
            if not misc.step_allready_completed(shortcuts.mergedFiles_list, "Picard MergeSamFiles"):
                misc.create_directory([f"{shortcuts.merged_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                misc.log_to_file("Starting: merging SAM/BAM files using Picard MergeSamFiles")
                cmd_merge = []
                with open(shortcuts.sortedFiles_list, 'r') as list:
                    for sample in list.read().splitlines():
                        if f"{options.tumor_id}." in sample: tumor = sample
                        else: normal = sample
                    cmd_merge.extend([f"picard MergeSamFiles {tumor} -O {shortcuts.merged_output_dir}{options.tumor_id}/{options.tumor_id}.bam", f"picard MergeSamFiles {normal} -O {shortcuts.merged_output_dir}{options.tumor_id}/{options.normal_id}.bam"])
                    misc.create_outputList_dna(shortcuts.mergedFiles_list, f"{options.tumor_id}.bam")
                    misc.create_outputList_dna(shortcuts.mergedFiles_list, f"{options.normal_id}.bam")
                    with mp.Pool() as pool:
                        pool.map(partial(self.multi_processing, [f"{shortcuts.merged_output_dir}{options.tumor_id}/{options.tumor_id}.bam", f"{shortcuts.merged_output_dir}{options.tumor_id}/{options.normal_id}.bam" ], options, misc, shortcuts),cmd_merge)
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'Picard MergeSamFiles succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                misc.run_command(f"rm {shortcuts.sorted_output_dir}{options.tumor_id}/*.bam", 'Removing sorted BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".merge() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def remove_duplicate(self, options, misc, shortcuts):
        '''This function removes duplicates '''

        try:
            if not misc.step_allready_completed(shortcuts.removeDuplicates_list, "Picard MarkDuplicates"):
                misc.create_directory([f"{shortcuts.removed_duplicates_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                misc.log_to_file("Starting: removing duplicates in SAM/BAM files using Picard MarkDuplicates")
                cmd_removedup = []
                with open(shortcuts.mergedFiles_list, 'r') as list:
                    tumor, normal = list.read().splitlines()
                    cmd_removedup.extend([f"picard -Xmx70g MarkDuplicates -I {shortcuts.merged_output_dir}{options.tumor_id}/{tumor} -O {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{tumor} -M {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/marked_dup_metrics_{tumor}.txt --TMP_DIR {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/tmp",
                                         f"picard -Xmx70g MarkDuplicates -I {shortcuts.merged_output_dir}{options.tumor_id}/{normal} -O {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{normal} -M {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/marked_dup_metrics_{normal}.txt --TMP_DIR {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/tmp"])
                    # Runs multiprocessing
                    with mp.Pool() as pool:
                        pool.map(partial(self.multi_processing, [f"{shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{tumor}", f"{shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{normal}"], options, misc, shortcuts),cmd_removedup)
                copy(shortcuts.mergedFiles_list, shortcuts.removeDuplicates_list) # just copying because the content will be the same
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'Picard MarkDuplicates succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                misc.run_command(f"rm {shortcuts.merged_output_dir}{options.tumor_id}/*.bam", 'Removing merged BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".remove_duplicate() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def realign(self, options, misc, shortcuts):
        '''This function realigns the bam files'''

        try:
            if not misc.step_allready_completed(shortcuts.realignedFiles_list, "GATK LeftAlignIndels"):
                misc.create_directory([f"{shortcuts.realigned_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                misc.log_to_file("Starting: realigning SAM/BAM files using GATK LeftAlignIndels")
                cmd_index = []
                cmd_leftAlignIndels = []
                with open(shortcuts.removeDuplicates_list, 'r') as list:
                    tumor, normal = list.read().splitlines()
                    cmd_index.extend([f"samtools index {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{tumor}", f"samtools index {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{normal}"])
                    with mp.Pool() as pool:
                        pool.map(partial(self.multi_processing, "N/A", options, misc, shortcuts),cmd_index)
                    cmd_leftAlignIndels.extend([f"gatk LeftAlignIndels -R {shortcuts.reference_genome_file} -I {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{tumor} -O {shortcuts.realigned_output_dir}{options.tumor_id}/{tumor}",
                                                f"gatk LeftAlignIndels -R {shortcuts.reference_genome_file} -I {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/{normal} -O {shortcuts.realigned_output_dir}{options.tumor_id}/{normal}"])
                    with mp.Pool() as pool:
                        pool.map(partial(self.multi_processing, "N/A", options, misc, shortcuts),cmd_leftAlignIndels)
                    copy(shortcuts.removeDuplicates_list, shortcuts.realignedFiles_list)
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'gatk LeftAlignIndels succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                misc.run_command(f"rm {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/*.bam", 'Removing remove_duplicate BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".realign() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def multi_processing(self, file, options, misc, shortcuts, input):
        '''This function runs pipeline steps in parallell'''
        try:
            if "picard MergeSamFiles" or "picard MarkDuplicates" in input:
                file = file[0] if options.tumor_id in input else file[1]
            misc.run_command(input, None, file, None)
        except Exception as e:
            misc.log_exception(".multi_processing() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def gatk_haplotype(self, options, misc, shortcuts):

        try:
            start = timeit.default_timer()
            if not misc.step_allready_completed(shortcuts.haplotypecaller_complete, "GATK haplotypeCaller"):
                misc.create_directory([f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/chunks/"])
                cmd_haplotypecaller = []
                misc.log_to_file("Starting: looking for SNV's using GATK HaplotypeCaller (multiprocessing)")
                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.read().splitlines()
                    for chunk in listdir(shortcuts.reference_genome_chunks_dir):
                        cmd_haplotypecaller.append(f"gatk HaplotypeCaller -R {shortcuts.reference_genome_file} -I {shortcuts.realigned_output_dir}{options.tumor_id}/{sample_1} -I {shortcuts.realigned_output_dir}{options.tumor_id}/{sample_2} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/chunks/{options.tumor_id}_{chunk}.vcf -L {shortcuts.reference_genome_chunks_dir}{chunk}")
                        misc.create_outputList_dna(shortcuts.gatk_chunks_list, f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/chunks/{options.tumor_id}_{chunk}.vcf")
                with mp.Pool(processes=60) as pool:
                    pool.map(partial(self.multi_processing, "gatk", options, misc, shortcuts),cmd_haplotypecaller)
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'gatk haplotypecaller (multiprocessing) succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 1 (snv calling) in dna_seq_analysis.py:", e)

        try:
            # Merge all vcf files
            cmd_merge = f"picard MergeVcfs -I {shortcuts.gatk_chunks_list} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}.vcf"
            misc.run_command(cmd_merge, None, f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}.vcf", None)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 2 (merge vcf) in dna_seq_analysis.py:", e)

        try:
            # Remove all reads with read depth less than 10, selects only snps, exludes normal samples
            cmd_filter_read_depth = f"bcftools view -i 'MIN(FMT/DP)>10' -m2 -M2 -v snps -s ^{options.normal_id} {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor.vcf"
            misc.run_command(cmd_filter_read_depth, "GATK haplotypeCaller step 3 (remove read depth < 10, selects only snps, exludes normal samples)", f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor.vcf", None)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 3 in dna_seq_analysis.py:", e)


        try:
            # select heterozygous genotype, excludes GT=1/2
            cmd_filter_het = f"bcftools view -g het -e 'GT=\"1/2\"' {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf"
            misc.run_command(cmd_filter_het, "GATK haplotypeCaller step 4 (select heterozygous genotype, excludes GT=1/2)", f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf", None)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 4 (select heterozygous genotype, excludes GT=1/2) in dna_seq_analysis.py:", e)


        try:
            # Annotate vcf file
            cmd_annotate = f'''java -Xmx4g -jar $HOME/anaconda3/envs/sequencing/share/snpeff-5.0-1/snpEff.jar \\
            -v GRCh38.99 -canon -noInteraction -noNextProt -noMotif -strict \\
            -onlyProtein {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf \\
            > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf'''
            misc.run_command(cmd_annotate, "GATK haplotypeCaller step 5 (annotate vcf file)", f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf", None)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 5 (annotate vcf file) in dna_seq_analysis.py:", e)


        try:
            # Index feature file
            cmd_indexFeatureFile = f"gatk IndexFeatureFile -I {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf"
            if misc.run_command(cmd_indexFeatureFile, "GATK haplotypeCaller step 6 (index fearure file)", f"{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf.idx", shortcuts.haplotypecaller_complete):
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'All steps in GATK HaplotypeCaller succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 6 (IndexFeatureFile) in dna_seq_analysis.py:", e)


    #---------------------------------------------------------------------------
    def delly(self, options, misc, shortcuts):
        '''This function creates an output directory and runs delly to call for somatic SNV's'''

        try:
            if not misc.step_allready_completed(shortcuts.delly_complete, "Delly SNV calling"):
                misc.create_directory([f"{shortcuts.delly_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                misc.log_to_file("Starting: looking for somatic SNV's using delly")
                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.read().splitlines()
                    cmd_delly_call = f"delly call -x {shortcuts.reference_genome_exclude_template_file} -g {shortcuts.reference_genome_file} -o {shortcuts.delly_output_dir}{options.tumor_id}/delly.bcf {shortcuts.realigned_output_dir}{options.tumor_id}/{sample_1} {shortcuts.realigned_output_dir}{options.tumor_id}/{sample_2}"
                    misc.run_command(cmd_delly_call, "Delly calling (step 1)", f"{shortcuts.delly_output_dir}{options.tumor_id}/delly.bcf", None)
                with open(f'{shortcuts.delly_output_dir}{options.tumor_id}/sample.tsv', 'w', newline='') as tsv:
                    tsv_output = csv.writer(tsv, delimiter='\t')
                    tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
                    tsv_output.writerow([f"{options.normal_id}", 'control'])
                cmd_dos2unix = f"dos2unix {shortcuts.delly_output_dir}{options.tumor_id}/sample.tsv"
                misc.run_command(cmd_dos2unix, None, None, None)

                # Filter bcf file to only show somatic mutations
                cmd_filter = f"delly filter -f somatic -o {shortcuts.delly_output_dir}{options.tumor_id}/delly_filter.bcf -s {shortcuts.delly_output_dir}{options.tumor_id}/sample.tsv {shortcuts.delly_output_dir}{options.tumor_id}/delly.bcf"
                misc.run_command(cmd_filter, "Delly filter (step 2)", f"{shortcuts.delly_output_dir}{options.tumor_id}/delly_filter.bcf", None)

                # Convert bcf file to human readable vcf file
                cmd_convert = f"bcftools view {shortcuts.delly_output_dir}{options.tumor_id}/delly_filter.bcf > {shortcuts.delly_output_dir}{options.tumor_id}/delly_filter.vcf"
                misc.run_command(cmd_convert, "Delly convert bcf > vcf (step 3)", f"{shortcuts.delly_output_dir}{options.tumor_id}/delly_filter.vcf", shortcuts.delly_complete)
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'Delly SNV calling succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".delly() in dna_seq_analysis.py:", e)
            sys.exit()

        #---------------------------------------------------------------------------
    def manta(self, options, misc, shortcuts):
        '''This function creates an output directory and runs manta to call for somatic SNV's'''

        try:
            if not misc.step_allready_completed(shortcuts.manta_complete, "Manta SNV calling"):
                start = timeit.default_timer()
                misc.log_to_file("Starting: looking for somatic SNV's using manta")
                with open(shortcuts.realignedFiles_list, 'r') as list:
                    sample_1, sample_2 = list.read().splitlines()
                    cmd_create_config_file = f"{shortcuts.configManta_file} --tumorBam={shortcuts.realigned_output_dir}{options.tumor_id}/{sample_1} --bam={shortcuts.realigned_output_dir}{options.tumor_id}/{sample_2} --referenceFasta={shortcuts.reference_genome_file} --runDir={shortcuts.manta_output_dir}{options.tumor_id}/"
                    misc.run_command(cmd_create_config_file, "Manta create config file (step 1)", shortcuts.runWorkflow_file, None)
                    cmd_runWorkflow = f"{shortcuts.runWorkflow_file} -m local -j 4"
                    misc.run_command(cmd_runWorkflow, 'Manta running workflow (step 2)', f"{shortcuts.manta_variants_dir}somaticSV.vcf.gz", None)
                    cmd_unzip = f'gunzip {shortcuts.manta_variants_dir}somaticSV.vcf.gz'
                    misc.run_command(cmd_unzip, 'Unzipping somaticSV.vcf.gz (step 3)', f"{shortcuts.manta_variants_dir}somaticSV.vcf", None)
                    cmd_filter = f'bcftools view -i \'FILTER=="PASS"\' {shortcuts.manta_variants_dir}somaticSV.vcf > {shortcuts.manta_variants_dir}somaticSV_PASS.vcf'
                    misc.run_command(cmd_filter, 'Filtering of passed SNV\'s (step 4)', f"{shortcuts.manta_variants_dir}somaticSV_PASS.vcf", shortcuts.manta_complete)
                    elapsed = timeit.default_timer() - start
                    misc.log_to_file(f'Manta SNV calling succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".manta() in dna_seq_analysis.py:", e)
            sys.exit()
