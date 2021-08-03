# Packages used in script
from os import listdir, sys, mkdir, getenv, path, rename, remove
import subprocess
import argparse
import csv
import multiprocessing as mp
from sys import path_importer_cache
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
        '''Runs bwa-mem2 index in the linux shell to index the reference genome'''
        misc.log_to_file("INFO", 'Starting: index_genome_dna')
        start = timeit.default_timer()

        try:
            # Create shortcuts
            misc.log_to_file("DEBUG", "# Create shortcuts")
            ref_file = shortcuts.reference_genome_file
            chunks_dir = shortcuts.reference_genome_chunks_dir
            allready_completed = shortcuts.bwa_index_complete

            # Create chunks directory
            misc.create_directory([chunks_dir])



            # Check if step is allready completed
            misc.log_to_file("DEBUG", "# Check if step is allready completed")
            if not misc.step_allready_completed(allready_completed, f"Indexing {ref_file.split('/')[-1]}"):
                misc.clear_screen()

                cmd_bwa_index = {   'cmd' : f"bwa index {ref_file}",
                                    'program' : 'bwa-mem2',
                                    'text' : 'Indexing genome with bwa',
                                    'file' : f'{ref_file}.sa.complete' }

                misc.run_command(cmd_bwa_index)

                cmd_create_dict = { 'cmd' : f"samtools dict {ref_file} -o {ref_file.split('.')[0]}.dict",
                                    'program' : 'bwa-mem2',
                                    'text' : 'Creating .dict with samtools dict',
                                    'file' : f'{ref_file[:-5]}dict.complete' }

                misc.run_command(cmd_create_dict)

                cmd_create_fai = {  'cmd' : f"samtools faidx {ref_file} -o {ref_file}.fai",
                                    'text' : 'Creating .fai with samtools faidx',
                                    'file' :  allready_completed } # This is the last step in indexing genome so here we can create an index.complete file
                
                misc.run_command(cmd_create_fai)

            # Splits ref_file.fai into chunks with 10 Megabases / chunk
            cmd_split_fasta = { 'cmd' : f"bedtools makewindows -w 10000000 -g {ref_file}.fai > {chunks_dir}chunk.bed",
                                'text' : 'Spliting fa.fai with bedtools makewindows' }



            misc.run_command(cmd_split_fasta)

            cmd_split_bed = {   'cmd' : f"split -l 2 {chunks_dir}chunk.bed {chunks_dir}chunk_",
                                'text' : 'Splitting chunk.bed' }
                              


            misc.run_command(cmd_split_bed)

            # add suffix .bed to split files
            for file in listdir(chunks_dir):
                if not file.endswith("bed"): rename(f"{chunks_dir}{file}", f"{chunks_dir}{file}.bed")
                if file == "chunk.bed": remove(f"{chunks_dir}{file}")
            misc.log_to_file("INFO", 'Renaming chunk_* > chunk_*.bed succesfully completed - OK!')
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Indexing reference genome successfully completed in {misc.elapsed_time(elapsed)} - OK!!')
            input("Press any key to return to DNA analysis menu...")
        except Exception as e:
            misc.log_to_file("ERROR", f"{e}: .index_genome_dna() in dna_seq_analysis.py:")

    #---------------------------------------------------------------------------
    def validate_bam_dna(self, options, misc, shortcuts):
        '''Runs picard ValidateSamFile to check if any errors are present in the aligned files. Returns True if no error were found'''

        try:
            misc.log_to_file("INFO", "Starting: Validating bam files")
            start = timeit.default_timer()
            cmd_validate = []
            target_dir = shortcuts.aligned_output_dir

            for sample in listdir(target_dir):
                path_to_sample = path.join(target_dir, sample)
                
                if path.isfile(path_to_sample) and sample.endswith('.bam'):
                    sample = { 'cmd' : f'java -Xmx60g -jar $HOME/anaconda3/envs/sequencing/share/picard-2.25.2-0/picard.jar ValidateSamFile -I {path_to_sample} --MODE SUMMARY --IGNORE_WARNINGS true --MAX_OPEN_TEMP_FILES 1000 --MAX_RECORDS_IN_RAM 10000000 --TMP_DIR {path_to_sample}/tmp_validate',
                               'program' : 'Picard ValidateSamFile',
                               'text' : f"Validating {sample}",
                               'file' : f'{path_to_sample}.validated'} 
                    cmd_validate.append(sample)
                else:
                    pass
                    
                    
            with mp.Pool(processes=3) as pool:
                pool.map(partial(misc.run_command), cmd_validate)
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'All .bam files succesfully validated in {misc.elapsed_time(elapsed)} - OK!')    
            return True

        except Exception as e:
            misc.log_exception(".validate_bam_dna() in dna_seq_analysis.py:", e)
            sys.exit()
    #---------------------------------------------------------------------------
    def alignment(self, options, misc, shortcuts):
        '''This function align reads to reference genome using Burrows Wheeler aligner (bwa-mem2)'''

        misc.log_to_file("INFO", f'Starting: Burrows Wheeler aligner Using {options.threads} out of {mp.cpu_count()} available threads')
        start = timeit.default_timer()

        try:
            misc.create_directory([f"{shortcuts.aligned_output_dir}"])

            with open(f'{shortcuts.dna_seq_dir}{options.tumor_id}_library.txt', 'r') as fastq_list:
                for line in fastq_list.readlines():
                    clinical_id, library_id, read1, read2 = line.split()
                    read_group_header = f'\'@RG\\tID:{library_id}\\tSM:{clinical_id}\\tLB:{library_id}\\tPL:ILLUMINA\\tPU:{library_id}\''
                    
                    if read2 == 'N/A': # single-end
                        cmd_bwa = { 'cmd' : f'set -o pipefail && bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}{read1} -t {options.threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam -', # samtools view converts SAM to BAM
                                    'program' : 'bwa-mem2',
                                    'text' : f"Aligning {library_id}.bam",
                                    'file' : f'{shortcuts.aligned_output_dir}{library_id}.bam.complete'} 

                    else: # paired-end
                        cmd_bwa = { 'cmd' : f'set -o pipefail && bwa mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}{read1} {shortcuts.dna_reads_dir}{read2} -t {options.threads} | samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam -',
                                    'program' : 'bwa-mem2',
                                    'text' : f"Aligning {library_id}.bam",
                                    'file' : f'{shortcuts.aligned_output_dir}{library_id}.bam.complete'} 
                        
                        # cmd_bwa = { 'cmd_1' : f'bwa-mem2 mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}{read1} {shortcuts.dna_reads_dir}{read2} -t {options.threads}',
                        #             'cmd_2' : f'samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam',
                        #             'text' : library_id,
                        #             'file' : f'{shortcuts.aligned_output_dir}{library_id}.bam.complete'} 

                    misc.run_command(cmd_bwa)       
                elapsed = timeit.default_timer() - start
                misc.log_to_file("INFO", f'Burrows Wheeler aligner succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        
        except OSError as e:
            misc.log_exception("You have to create a library file first", e)
            sys.exit()
        
        except Exception as e:
            misc.log_exception("alignment() in dna_seq_analysis.py:", e)
            sys.exit()


    #---------------------------------------------------------------------------
    def sort(self, options, misc, shortcuts):
        '''This function sorts aligned bam files in parallel with mp.Pool'''

        try:
            misc.create_directory([f"{shortcuts.sorted_output_dir}"])
            start = timeit.default_timer()
            misc.log_to_file("INFO", "Starting: sorting SAM/BAM files using Picard Sortsam")
            
            # Empty strings to store the output
            cmd_sort = []
            tumor_sort_str = ""
            normal_sort_str = ""
            write_to_file = ""
            target_dir = shortcuts.aligned_output_dir

            for sample in listdir(target_dir):
                path_to_sample = path.join(target_dir, sample)
                
                if path.isfile(path_to_sample) and sample.endswith('.bam'):
                    input = { 'cmd' : f'java -Xmx60g -jar $HOME/anaconda3/envs/sequencing/share/picard-2.25.2-0/picard.jar SortSam -I {path_to_sample} -O {shortcuts.sorted_output_dir}{sample} --SORT_ORDER coordinate --MAX_RECORDS_IN_RAM 21000000 --TMP_DIR {shortcuts.sorted_output_dir}tmp',
                               'program' : 'Picard SortSam',
                               'text' : f"Sorting {sample}",
                               'file' : f'{shortcuts.sorted_output_dir}{sample}.complete'} 
                    cmd_sort.append(input)
                    
            #         if options.tumor_id in sample:
            #             tumor_sort_str += f" -I {shortcuts.sorted_output_dir}{sample}".rstrip()
                    
            #         else:
            #             normal_sort_str += f" -I {shortcuts.sorted_output_dir}{sample}".rstrip()
            
            # write_to_file = f"{tumor_sort_str.lstrip()}\n{normal_sort_str.lstrip()}"
               
            with mp.Pool(processes=2) as pool:
                pool.map(partial(misc.run_command),cmd_sort)
            misc.create_outputList_dna(shortcuts.sortedFiles_list, write_to_file)
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Picard SortSam succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
           
            # TODO, doesn't work 
            # misc.run_command({'cmd' : f'rm {shortcuts.aligned_output_dir}*.bam', 'text' : 'Removing aligned BAM files to save space', 'program' : 'Removing aligned BAM files to save space'})
           
        except Exception as e:
            misc.log_exception(".sort() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def merge(self, options, misc, shortcuts):
        '''This function merges all the input files in the sortedFiles_list to one output file'''

        try:
            if not misc.step_allready_completed(shortcuts.mergedFiles_list, "Picard MergeSamFiles"):
                misc.create_directory([f"{shortcuts.merged_output_dir}"])
                start = timeit.default_timer()
                misc.log_to_file("INFO", "Starting: merging SAM/BAM files using Picard MergeSamFiles")  
                cmd_merge = []
                tumor = ""
                normal = ""
                target_dir = shortcuts.sorted_output_dir

                for sample in listdir(target_dir):
                    path_to_sample = path.join(target_dir, sample)
                
                    if path.isfile(path_to_sample) and sample.endswith('.bam'):

                        if f"{options.tumor_id}." in sample:
                            tumor += f"-I {path_to_sample} " 

                        else:
                            normal += f"-I {path_to_sample} " 
     
                cmd_merge.extend([{ 'cmd' : f'picard MergeSamFiles {tumor} -O {shortcuts.merged_output_dir}{options.tumor_id}.bam',
                                    'program' : 'Picard MergeSamFiles',
                                    'text' : f'Merging {options.tumor_id}.bam',
                                    'file' : f'{shortcuts.merged_output_dir}{options.tumor_id}.bam.complete'},

                                 {  'cmd' : f'picard MergeSamFiles {normal} -O {shortcuts.merged_output_dir}{options.normal_id}.bam',
                                    'program' : 'Picard MergeSamFiles',
                                    'text' : f'Merging {options.normal_id}.bam',
                                    'file' : f'{shortcuts.merged_output_dir}{options.normal_id}.bam.complete'}])

                with mp.Pool(processes=2) as pool:
                    pool.map(partial(misc.run_command), cmd_merge)
                elapsed = timeit.default_timer() - start
                misc.log_to_file("INFO", f'Picard MergeSamFiles succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                # TODO
                #misc.run_command(f"rm {shortcuts.sorted_output_dir}{options.tumor_id}/*.bam", 'Removing sorted BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".merge() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def remove_duplicate(self, options, misc, shortcuts):
        '''This function removes duplicates '''

        try:
            if not misc.step_allready_completed(shortcuts.removeDuplicates_list, "Picard MarkDuplicates"):
                misc.create_directory([f"{shortcuts.removed_duplicates_output_dir}"])
                start = timeit.default_timer()
                misc.log_to_file("INFO", "Starting: removing duplicates in SAM/BAM files using Picard MarkDuplicates")
                cmd_removedup = []
                target_dir = shortcuts.merged_output_dir

                for sample in listdir(target_dir):
                    path_to_sample = path.join(target_dir, sample)
                
                    if path.isfile(path_to_sample) and sample.endswith('.bam'):
                        cmd_removedup.append({ 'cmd' : f'picard -Xmx20g MarkDuplicates -I {path_to_sample} -O {shortcuts.removed_duplicates_output_dir}{sample} -M {shortcuts.removed_duplicates_output_dir}marked_dup_metrics_{sample}.txt --TMP_DIR {shortcuts.removed_duplicates_output_dir}tmp',
                                'program' : 'Picard MarkDuplicates',
                                'text' : f'Removing duplicates for {sample}',
                                'file' : f'{shortcuts.removed_duplicates_output_dir}{sample}.complete'})    
     
             

                # Runs multiprocessing
                with mp.Pool(processes=2) as pool:
                    pool.map(partial(misc.run_command), cmd_removedup)
                
                elapsed = timeit.default_timer() - start
                misc.log_to_file("INFO", f'Picard MarkDuplicates succesfully completed in {misc.elapsed_time(elapsed)} - OK!')

                # TODO
                # misc.run_command(f"rm {shortcuts.merged_output_dir}{options.tumor_id}/*.bam", 'Removing merged BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".remove_duplicate() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def realign(self, options, misc, shortcuts):
        '''This function realigns the bam files'''

        try:
            if not misc.step_allready_completed(shortcuts.realignedFiles_list, "GATK LeftAlignIndels"):
                misc.create_directory([f"{shortcuts.realigned_output_dir}"])
                start = timeit.default_timer()
                misc.log_to_file("INFO", "Starting: realigning SAM/BAM files using GATK LeftAlignIndels")
                cmd_index = []
                cmd_leftAlignIndels = []

                target_dir = shortcuts.removed_duplicates_output_dir
                for sample in listdir(target_dir):
                    path_to_sample = path.join(target_dir, sample)
                    if path.isfile(path_to_sample) and sample.endswith('.bam'):
                         cmd_index.append({ 'cmd' : f'samtools index {path_to_sample}',
                                'program' : 'Realign index',
                                'text' : f'Indexing {sample}',
                                'file' : f'{path_to_sample}.bai.complete'})  

                         cmd_leftAlignIndels.append({ 'cmd' : f'gatk LeftAlignIndels -R {shortcuts.reference_genome_file} -I {path_to_sample} -O {shortcuts.realigned_output_dir}{sample}',
                                'program' : 'gatk LeftAlignIndels',
                                'text' : f'Realigning {sample}',
                                'file' : f'{shortcuts.realigned_output_dir}{sample}.complete'})  
                    
                with mp.Pool(processes=2) as pool:
                    pool.map(partial(misc.run_command), cmd_index)
                    
                with mp.Pool(processes=2) as pool:
                    pool.map(partial(misc.run_command), cmd_leftAlignIndels)
                    
                elapsed = timeit.default_timer() - start
                misc.log_to_file("INFO", f'gatk LeftAlignIndels succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                # TODO
                # misc.run_command(f"rm {shortcuts.removed_duplicates_output_dir}{options.tumor_id}/*.bam", 'Removing remove_duplicate BAM files to save space', None, None)
        except Exception as e:
            misc.log_exception(".realign() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def gatk_haplotype(self, options, misc, shortcuts):

        try:
            start = timeit.default_timer()
            if not misc.step_allready_completed(shortcuts.haplotypecaller_complete, "GATK haplotypeCaller"):
                misc.create_directory([f"{shortcuts.haplotypecaller_output_dir}chunks/"])
                cmd_haplotypecaller = []
                misc.log_to_file("INFO", "Starting: looking for SNV's using GATK HaplotypeCaller (multiprocessing)")

                samples = []
                target_dir = shortcuts.realigned_output_dir
                for sample in sorted(listdir(target_dir)):
                    path_to_sample = path.join(target_dir, sample)
                    
                    if path.isfile(path_to_sample) and sample.endswith('.bam'):
                        samples.append(path_to_sample)
                    
                sample_1, sample_2 = samples
                for chunk in sorted(listdir(shortcuts.reference_genome_chunks_dir)):
                    path_to_chunk = path.join(shortcuts.reference_genome_chunks_dir, chunk)
                    cmd_haplotypecaller.append({    'cmd' : f'gatk HaplotypeCaller -R {shortcuts.reference_genome_file} -I {sample_1} -I {sample_2} -O {shortcuts.haplotypecaller_output_dir}chunks/{options.tumor_id}_{chunk}.vcf -L {path_to_chunk}',
                                                    'program' : 'gatk HaplotypeCaller',
                                                    'text' : None,
                                                    'file' : f'{shortcuts.haplotypecaller_output_dir}chunks/{options.tumor_id}_{chunk}.vcf.complete'})
                
                    # Creates a list for all chunks for picard MergeVcfs to use when merging all files in next step
                    misc.create_outputList_dna(shortcuts.gatk_chunks_list, f"{shortcuts.haplotypecaller_output_dir}chunks/{options.tumor_id}_{chunk}.vcf")

                with mp.Pool(processes=20) as pool:
                    pool.map(partial(misc.run_command), cmd_haplotypecaller)

                elapsed = timeit.default_timer() - start
                misc.log_to_file("INFO", f'gatk haplotypecaller step 1 (multiprocessing) succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 1 (snv calling) in dna_seq_analysis.py:", e)

        try:
            # Merge all vcf files
            cmd_merge = {   'cmd' : f'picard MergeVcfs -I {shortcuts.gatk_chunks_list} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf',
                            'program' : 'gatk HaplotypeCaller',
                            'text' : 'Merging vcf files',
                            'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf.complete'} 
            misc.run_command(cmd_merge)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 2 (merge vcf) in dna_seq_analysis.py:", e)

        try:
            # Remove all reads with read depth less than 10, selects only snps, exludes normal samples
            cmd_filter_read_depth = {   'cmd' : f"bcftools view -i 'MIN(FMT/DP)>10' -m2 -M2 -v snps -s ^{options.normal_id} {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf",
                                        'program' : 'gatk HaplotypeCaller',
                                        'text' : 'GATK haplotypeCaller step 3 (remove read depth < 10, selects only snps, exludes normal samples)',
                                        'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf.complete'} 
           
           
            misc.run_command(cmd_filter_read_depth)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 3 in dna_seq_analysis.py:", e)


        try:
            # select heterozygous genotype, excludes GT=1/2
            cmd_filter_het = {  'cmd' : f"bcftools view -g het -e 'GT=\"1/2\"' {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf",
                                'program' : 'gatk HaplotypeCaller',
                                'text' : 'GATK haplotypeCaller step 4 (select heterozygous genotype, excludes GT=1/2)',
                                'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf.complete'}

            misc.run_command(cmd_filter_het)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 4 (select heterozygous genotype, excludes GT=1/2) in dna_seq_analysis.py:", e)


        try:
            # Annotate vcf file
            cmd_annotate = {    'cmd' : f'''java -Xmx4g -jar $HOME/anaconda3/envs/sequencing/share/snpeff-5.0-1/snpEff.jar \\
                                        -v GRCh38.99 -canon -noInteraction -noNextProt -noMotif -strict \\
                                        -onlyProtein {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf \\
                                        > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf''',
                                'program' : 'gatk HaplotypeCaller',
                                'text' : 'GATK haplotypeCaller step 5 (annotate vcf file)',
                                'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf.complete'}           
                                        
            misc.run_command(cmd_annotate)
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 5 (annotate vcf file) in dna_seq_analysis.py:", e)


        try:
            # Index feature file
            cmd_indexFeatureFile = {    'cmd' : f"gatk IndexFeatureFile -I {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf",
                                        'program' : 'gatk HaplotypeCaller',
                                        'text' : 'GATK haplotypeCaller step 6 (index fearure file)',
                                        'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf.idx.complete'}   



            if misc.run_command(cmd_indexFeatureFile):
                elapsed = timeit.default_timer() - start
                misc.log_to_file("INFO", f'All steps in GATK HaplotypeCaller succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".gatk_haplotype step 6 (IndexFeatureFile) in dna_seq_analysis.py:", e)


    #---------------------------------------------------------------------------
    def delly(self, options, misc, shortcuts):
        '''This function creates an output directory and runs delly to call for somatic SNV's'''

        try:
            if not misc.step_allready_completed(shortcuts.delly_complete, "Delly SNV calling"):
                misc.create_directory([f"{shortcuts.delly_output_dir}{options.tumor_id}/"])
                start = timeit.default_timer()
                misc.log_to_file("INFO", "Starting: looking for somatic SNV's using delly")
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
                misc.log_to_file("INFO", f'Delly SNV calling succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".delly() in dna_seq_analysis.py:", e)
            sys.exit()

        #---------------------------------------------------------------------------
    def manta(self, options, misc, shortcuts):
        '''This function creates an output directory and runs manta to call for somatic SNV's'''

        try:
            if not misc.step_allready_completed(shortcuts.manta_complete, "Manta SNV calling"):
                start = timeit.default_timer()
                misc.log_to_file("INFO", "Starting: looking for somatic SNV's using manta")
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
                    misc.log_to_file("INFO", f'Manta SNV calling succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".manta() in dna_seq_analysis.py:", e)
            sys.exit()
