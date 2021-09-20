#!~/anaconda3/envs/sequencing/bin/python3

from os import listdir, sys, mkdir, getenv, path, rename, remove
from subprocess import run
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
    '''All steps in the DNA sequencing analysis pipeline are saved in different methods in this class'''

    def __init__(self):
        pass


    #---------------------------------------------------------------------------
    def index_genome_dna(self, misc, shortcuts):
        '''Runs bwa-mem2 index through the linux shell buy calling the misc.run_command method'''
        
        misc.log_to_file("INFO", 'Starting: index_genome_dna')
        start = timeit.default_timer()

        try:
            # Create shortcuts
            ref_file = shortcuts.reference_genome_file
            chunks_dir = shortcuts.reference_genome_chunks_dir

            # Create chunks directory
            misc.create_directory([chunks_dir])

            # Index reference genome
            cmd_bwa_index = {   'cmd' : f"bwa-mem2 index {ref_file}",
                                'program' : 'bwa-mem2',
                                'text' : 'Indexing genome with bwa',
                                'file' : f'{ref_file}.sa.complete' }
            misc.run_command(cmd_bwa_index)

            # Create .dict file
            cmd_create_dict = { 'cmd' : f"samtools dict {ref_file} -o {ref_file[:-2]}.dict",
                                'program' : 'bwa-mem2',
                                'text' : 'Creating .dict with samtools dict',
                                'file' : f"{ref_file}.dict.complete" }
            misc.run_command(cmd_create_dict)

            # Create .fai file
            cmd_create_fai = {  'cmd' : f"samtools faidx {ref_file} -o {ref_file}.fai",
                                'text' : 'Creating .fai with samtools faidx',
                                'file' :  f"{ref_file}.fai.complete" }
            misc.run_command(cmd_create_fai)

            # Splits ref_file.fai into chunks with 10 Megabases / chunk
            
            cmd_split_fasta = { 'cmd' : f"bedtools makewindows -w 20000000 -g {ref_file}.fai > {chunks_dir}chunk.bed",
                                'text' : 'Spliting fa.fai with bedtools makewindows' }
            misc.run_command(cmd_split_fasta)

           
            cmd_split_bed = {   'cmd' : f"split -l 3 {chunks_dir}chunk.bed {chunks_dir}chunk_",
                                'text' : 'Splitting chunk.bed' }
            misc.run_command(cmd_split_bed)

            # adds suffix .bed to split files
            [rename(f"{chunks_dir}{file}", f"{chunks_dir}{file}.bed") for file in listdir(chunks_dir) if not file.endswith("bed")]
            remove(f"{chunks_dir}chunk.bed")


            # for file in listdir(chunks_dir):
            #     if not file.endswith("bed"): rename(f"{chunks_dir}{file}", f"{chunks_dir}{file}.bed")
            #     if file == "chunk.bed": remove(f"{chunks_dir}{file}")
            misc.log_to_file("INFO", 'Renaming chunk_* > chunk_*.bed succesfully completed - OK!')
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Indexing reference genome successfully completed in {misc.elapsed_time(elapsed)} - OK!!')
            input("Press any key to return to DNA analysis menu...")
        except Exception as e:
            misc.log_to_file("ERROR", f"{e}: .index_genome_dna() in dna_seq_analysis.py:")

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
                        cmd_bwa = { 'cmd' : f'set -o pipefail && bwa-mem2 mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}{read1} -t {options.threads} | samtools view -bS --threads 1 -o {shortcuts.aligned_output_dir}{library_id}.bam -', # samtools view converts SAM to BAM
                                    'program' : 'bwa-mem2',
                                    'text' : f"Aligning {library_id}.bam",
                                    'file' : f'{shortcuts.aligned_output_dir}{library_id}.bam.complete'} 

                    else: # paired-end
                        cmd_bwa = { 'cmd' : f'set -o pipefail && bwa-mem2 mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}{read1} {shortcuts.dna_reads_dir}{read2} -t {options.threads} | samtools view -bS --threads 1 -o {shortcuts.aligned_output_dir}{library_id}.bam -',
                                    'program' : 'bwa-mem2',
                                    'text' : f"Aligning {library_id}.bam",
                                    'file' : f'{shortcuts.aligned_output_dir}{library_id}.bam.complete'}


                        # for running without shell=True
                        # cmd_bwa = { 'cmd_1' : f'bwa-mem2 mem -R {read_group_header} {shortcuts.reference_genome_file} {shortcuts.dna_reads_dir}{read1} {shortcuts.dna_reads_dir}{read2} -t {options.threads}',
                        #             'cmd_2' : f'samtools view -bS -o {shortcuts.aligned_output_dir}{library_id}.bam',
                        #             'text' : library_id,
                        #             'file' : f'{shortcuts.aligned_output_dir}{library_id}.bam.complete'} 

                    misc.run_command(cmd_bwa)       
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Burrows Wheeler aligner succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        
        except OSError as e:
            misc.log_to_file("ERROR", "You have to create a library file first", e)
            sys.exit()
        
        except Exception as e:
            misc.log_to_file("ERROR", "alignment() in dna_seq_analysis.py:", e)
            sys.exit()

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
                    sample = { 'cmd' : f'java -Xmx60g -jar $HOME/anaconda3/envs/sequencing/share/picard-2.25.2-0/picard.jar ValidateSamFile -I {path_to_sample} --MODE SUMMARY --IGNORE_WARNINGS true --MAX_OPEN_TEMP_FILES 2000 --MAX_RECORDS_IN_RAM 20000000 --TMP_DIR {path_to_sample}/tmp_validate',
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
            misc.log_to_file("ERROR", ".validate_bam_dna() in dna_seq_analysis.py:", e)
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
                    input = { 'cmd' : f'sambamba sort {path_to_sample} -o {shortcuts.sorted_output_dir}{sample} -p -t 12 -m 65GB --tmpdir {shortcuts.sorted_output_dir}tmp',
                               'program' : 'Sambamba sort',
                               'text' : f"Sorting {sample}",
                               'file' : f'{shortcuts.sorted_output_dir}{sample}.complete',
                               'old_file': f'{path_to_sample}'} 

                    # input = { 'cmd' : f'java -Xmx60g -jar $HOME/anaconda3/envs/sequencing/share/picard-2.25.2-0/picard.jar SortSam -I {path_to_sample} -O {shortcuts.sorted_output_dir}{sample} --SORT_ORDER coordinate --MAX_RECORDS_IN_RAM 11000000 --TMP_DIR {shortcuts.sorted_output_dir}tmp',
                    #            'program' : 'Picard SortSam',
                    #            'text' : f"Sorting {sample}",
                    #            'file' : f'{shortcuts.sorted_output_dir}{sample}.complete',
                    #            'old_file': f'{path_to_sample}'} 
                    cmd_sort.append(input)
               
            with mp.Pool(processes=3) as pool:
                pool.map(partial(misc.run_command),cmd_sort)
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Picard SortSam succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
           
            
           
        except Exception as e:
            misc.log_to_file("ERROR", ".sort() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def merge(self, options, misc, shortcuts):
        '''This function merges all the input files in the sorted_output_dir to one output file'''

        try:
            misc.create_directory([f"{shortcuts.merged_output_dir}"])
            start = timeit.default_timer()
            misc.log_to_file("INFO", "Starting: merging SAM/BAM files using Picard MergeSamFiles")  
            cmd_merge = []
            tumor = ""
            normal = ""
            target_dir = shortcuts.sorted_output_dir


            # Makes a list of all files in sorted_output_dir and loop through everyone of them
            for sample in listdir(target_dir):
                path_to_sample = path.join(target_dir, sample)
                
                if path.isfile(path_to_sample) and sample.endswith('.bam'):

                    if f"{options.tumor_id}." in sample:
                        tumor += f"{path_to_sample} " 

                        # For picard MergeSamFiles
                        # tumor += f"-I {path_to_sample} " 

                    else:
                        misc.run_command({'cmd' : f'mv {path_to_sample} {shortcuts.merged_output_dir}/{options.normal_id}.bam',
                                          'text' : f'Moving {path_to_sample} to {shortcuts.merged_output_dir}/{options.normal_id}.bam',
                                          'file' : f'{shortcuts.merged_output_dir}{options.normal_id}.bam.complete'})

            
            # Sambamba merge
            input = {'cmd' : f'sambamba merge -t 38 {shortcuts.merged_output_dir}{options.tumor_id}.bam {tumor}',
                     'program' : 'Sambamba merge',
                     'text' : f'Merging {options.tumor_id}.bam',
                     'file' : f'{shortcuts.merged_output_dir}{options.tumor_id}.bam.complete',
                     'old_file' : f'{shortcuts.sorted_output_dir}*.bam'}

            # Picard MergeSamFiles
            # input = {'cmd' : f'picard MergeSamFiles {tumor} -O {shortcuts.merged_output_dir}{options.tumor_id}.bam --USE_THREADING True',
            #          'program' : 'Picard MergeSamFiles',
            #          'text' : f'Merging {options.tumor_id}.bam',
            #          'file' : f'{shortcuts.merged_output_dir}{options.tumor_id}.bam.complete',
            #          'old_file' : f'{shortcuts.sorted_output_dir}*.bam'}


            
            misc.run_command(input)


            # cmd_merge.extend([{ 'cmd' : f'picard MergeSamFiles {tumor} -O {shortcuts.merged_output_dir}{options.tumor_id}.bam',
            #                     'program' : 'Picard MergeSamFiles',
            #                     'text' : f'Merging {options.tumor_id}.bam',
            #                     'file' : f'{shortcuts.merged_output_dir}{options.tumor_id}.bam.complete',
            #                     'old_file' : f'{shortcuts.sorted_output_dir}{options.tumor_id}/*.bam}'},

                            #  {  'cmd' : f'picard MergeSamFiles {normal} -O {shortcuts.merged_output_dir}{options.normal_id}.bam',
                            #     'program' : 'Picard MergeSamFiles',
                            #     'text' : f'Merging {options.normal_id}.bam',
                            #     'file' : f'{shortcuts.merged_output_dir}{options.normal_id}.bam.complete'}])

            # with mp.Pool(processes=2) as pool:
            #     pool.map(partial(misc.run_command), cmd_merge)
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Picard MergeSamFiles succesfully completed in {misc.elapsed_time(elapsed)} - OK!')

           
        except Exception as e:
            misc.log_to_file("ERROR", f".merge() in dna_seq_analysis.py: {e}")
            sys.exit()

    #---------------------------------------------------------------------------
    def remove_duplicate(self, options, misc, shortcuts):
        '''This function removes duplicates '''

        try:
            misc.create_directory([f"{shortcuts.removed_duplicates_output_dir}"])
            start = timeit.default_timer()
            misc.log_to_file("INFO", "Starting: removing duplicates in SAM/BAM files using Picard MarkDuplicates")
            cmd_removedup = []
            target_dir = shortcuts.merged_output_dir
            for sample in listdir(target_dir):
                path_to_sample = path.join(target_dir, sample)
            
                if path.isfile(path_to_sample) and sample.endswith('.bam'):
                    cmd_removedup.append({ 'cmd' : f'sambamba markdup --sort-buffer-size 153600 -t 19 --tmpdir {shortcuts.removed_duplicates_output_dir}tmp --overflow-list-size 2000000 -p {path_to_sample} {shortcuts.removed_duplicates_output_dir}{sample}',
                            'program' : 'Sambamba MarkDuplicates',
                            'text' : f'Removing duplicates for {sample}',
                            'file' : f'{shortcuts.removed_duplicates_output_dir}{sample}.complete',
                            'old_file': f'{path_to_sample}' })    
                    # cmd_removedup.append({ 'cmd' : f'picard -Xmx85g MarkDuplicates -I {path_to_sample} -O {shortcuts.removed_duplicates_output_dir}{sample} --ASSUME_SORTED -M {shortcuts.removed_duplicates_output_dir}marked_dup_metrics_{sample}.txt --TMP_DIR {shortcuts.removed_duplicates_output_dir}tmp',
                    #         'program' : 'Picard MarkDuplicates',
                    #         'text' : f'Removing duplicates for {sample}',
                    #         'file' : f'{shortcuts.removed_duplicates_output_dir}{sample}.complete',
                    #         'old_file': f'{shortcuts.merged_output_dir}*.bam'})    
    
            
            # Runs multiprocessing
            with mp.Pool(processes=2) as pool:
                pool.map(partial(misc.run_command), cmd_removedup)
            
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Picard MarkDuplicates succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
            
        except Exception as e:
            misc.log_to_file("ERROR", ".remove_duplicate() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def realign(self, options, misc, shortcuts):
        '''This function realigns the bam files'''

        try:
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
                            'program' : 'realign index',
                            'text' : f'Indexing {sample}',
                            'file' : f'{path_to_sample}.bai.complete'})  
                     cmd_leftAlignIndels.append({ 'cmd' : f'gatk LeftAlignIndels -R {shortcuts.reference_genome_file} -I {path_to_sample} -O {shortcuts.realigned_output_dir}{sample}',
                            'program' : 'realign',
                            'text' : f'Realigning {sample}',
                            'file' : f'{shortcuts.realigned_output_dir}{sample}.complete',
                            'old_file': f'{shortcuts.removed_duplicates_output_dir}*.bam'})  
                
            with mp.Pool(processes=2) as pool:
                pool.map(partial(misc.run_command), cmd_index)
                
            with mp.Pool(processes=2) as pool:
                pool.map(partial(misc.run_command), cmd_leftAlignIndels)
                
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'gatk LeftAlignIndels succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
            
        except Exception as e:
            misc.log_to_file("ERROR", ".realign() in dna_seq_analysis.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def gatk_haplotype(self, options, misc, shortcuts):

        try:
            start = timeit.default_timer()
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
            with mp.Pool(processes=60) as pool:
                pool.map(partial(misc.run_command), cmd_haplotypecaller)
            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'gatk haplotypecaller step 1 (multiprocessing) succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_to_file("ERROR", ".gatk_haplotype step 1 (snv calling) in dna_seq_analysis.py:", e)

        try:
            # Merge all vcf files
            cmd_merge = {   'cmd' : f'picard MergeVcfs -I {shortcuts.gatk_chunks_list} -O {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf',
                            'program' : 'gatk HaplotypeCaller',
                            'text' : 'Merging vcf files',
                            'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf.complete'} 
            misc.run_command(cmd_merge)
        except Exception as e:
            misc.log_to_file("ERROR", ".gatk_haplotype step 2 (merge vcf) in dna_seq_analysis.py:", e)

        try:
            # Remove all reads with read depth less than 10, selects only snps, exludes normal samples
            cmd_filter_read_depth = {   'cmd' : f"bcftools view -i 'MIN(FMT/DP)>10' -m2 -M2 -v snps -s ^{options.normal_id} {shortcuts.haplotypecaller_output_dir}{options.tumor_id}.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf",
                                        'program' : 'gatk HaplotypeCaller',
                                        'text' : 'GATK haplotypeCaller step 3 (remove read depth < 10, selects only snps, exludes normal samples)',
                                        'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf.complete'} 
           
           
            misc.run_command(cmd_filter_read_depth)
        except Exception as e:
            misc.log_to_file("ERROR", ".gatk_haplotype step 3 in dna_seq_analysis.py:", e)


        try:
            # select heterozygous genotype, excludes GT=1/2
            cmd_filter_het = {  'cmd' : f"bcftools view -g het -e 'GT=\"1/2\"' {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor.vcf > {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf",
                                'program' : 'gatk HaplotypeCaller',
                                'text' : 'GATK haplotypeCaller step 4 (select heterozygous genotype, excludes GT=1/2)',
                                'file' : f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf.complete'}

            misc.run_command(cmd_filter_het)
        except Exception as e:
            misc.log_to_file("ERROR", ".gatk_haplotype step 4 (select heterozygous genotype, excludes GT=1/2) in dna_seq_analysis.py:", e)


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
            misc.log_to_file("ERROR", ".gatk_haplotype step 5 (annotate vcf file) in dna_seq_analysis.py:", e)


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
            misc.log_to_file("ERROR", ".gatk_haplotype step 6 (IndexFeatureFile) in dna_seq_analysis.py:", e)


    #---------------------------------------------------------------------------
    def delly(self, options, misc, shortcuts):
        '''This function creates an output directory and runs delly to call for somatic SNV's'''

        try:
            misc.create_directory([shortcuts.delly_output_dir])
            start = timeit.default_timer()
            misc.log_to_file("INFO", "Starting: looking for somatic SNV's using delly")
           

            # Gets path to the samples from realigned output
            samples = []
            target_dir = shortcuts.realigned_output_dir
            for sample in sorted(listdir(target_dir)):
                path_to_sample = path.join(target_dir, sample)
                
                if path.isfile(path_to_sample) and sample.endswith('.bam'):
                    samples.append(path_to_sample)
            
            sample_1, sample_2 = samples
            
            cmd_delly_call = { 'cmd' : f"delly call -x {shortcuts.reference_genome_exclude_template_file} -g {shortcuts.reference_genome_file} -o {shortcuts.delly_output_dir}delly.bcf {sample_1} {sample_2}",
                               'text' : "Delly calling",
                               'file' : f"{shortcuts.delly_output_dir}delly.bcf.complete" }
            misc.run_command(cmd_delly_call)
            
            with open(f'{shortcuts.delly_output_dir}sample.tsv', 'w', newline='') as tsv:
                tsv_output = csv.writer(tsv, delimiter='\t')
                tsv_output.writerow([f"{options.tumor_id}", 'tumor'])
                tsv_output.writerow([f"{options.normal_id}", 'control'])
            cmd_dos2unix = { 'cmd' : f"dos2unix {shortcuts.delly_output_dir}sample.tsv",
                             'text' : "Converting sample.tsv with dos2unix" }
            misc.run_command(cmd_dos2unix)

            # Filter bcf file to only show somatic mutations
            cmd_filter = { 'cmd' : f"delly filter -f somatic -o {shortcuts.delly_output_dir}delly_filter.bcf -s {shortcuts.delly_output_dir}sample.tsv {shortcuts.delly_output_dir}delly.bcf",
                           'text' : "Filtering delly.bcf to show only somatic mutations",
                           'file' : f"{shortcuts.delly_output_dir}delly_filter.bcf.complete" }
            misc.run_command(cmd_filter)
            
            # Convert bcf file to human readable vcf file
            cmd_convert = { 'cmd' : f"bcftools view {shortcuts.delly_output_dir}delly_filter.bcf > {shortcuts.delly_output_dir}delly_filter.vcf",
                            'text' : "Converting delly.bcf to a human readable vcf file",
                            'file' : f"{shortcuts.delly_output_dir}delly_filter.vcf.complete" }
            misc.run_command(cmd_convert)

            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Delly SNV calling succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        
        except Exception as e:
            misc.log_to_file("ERROR", ".delly() in dna_seq_analysis.py:", e)
            sys.exit()

        #---------------------------------------------------------------------------
    def manta(self, options, misc, shortcuts):
        '''This function creates an output directory and runs manta to call for somatic SNV's'''

        try:
            misc.create_directory([f"{shortcuts.manta_output_dir}"])
            start = timeit.default_timer()
            misc.log_to_file("INFO", "Starting: looking for somatic SNV's using manta")
            

            # Gets path to the samples from realigned output
            samples = []
            target_dir = shortcuts.realigned_output_dir
            for sample in sorted(listdir(target_dir)):
                path_to_sample = path.join(target_dir, sample)
                
                if path.isfile(path_to_sample) and sample.endswith('.bam'):
                    samples.append(path_to_sample)
            
            sample_1, sample_2 = samples

            cmd_create_config_file = { 'cmd' : f"{shortcuts.configManta_file} --tumorBam={sample_1} --bam={sample_2} --referenceFasta={shortcuts.reference_genome_file} --runDir={shortcuts.manta_output_dir}",
                                       'text' : "Manta create config file (step 1)",
                                       'file' : f"{shortcuts.manta_output_dir}configManta_file.complete" }
            misc.run_command(cmd_create_config_file)

            cmd_runWorkflow = { 'cmd' : f"{shortcuts.runWorkflow_file} -m local -j 8",
                                'text' : "Manta running workflow (step 2)",
                                'file' : f"{shortcuts.manta_variants_dir}somaticSV.vcf.gz.complete" }
            misc.run_command(cmd_runWorkflow)

            cmd_unzip = { 'cmd' : f'gunzip {shortcuts.manta_variants_dir}somaticSV.vcf.gz',
                          'text' : "Unzipping somaticSV.vcf.gz (step 3)",
                          'file' : f"{shortcuts.manta_variants_dir}somaticSV.vcf.complete" }
            misc.run_command(cmd_unzip)

            cmd_filter = { 'cmd' : f'bcftools view -i \'FILTER=="PASS"\' {shortcuts.manta_variants_dir}somaticSV.vcf > {shortcuts.manta_variants_dir}somaticSV_PASS.vcf',
                           'text' : "Filtering of passed SNV\'s (step 4)",
                           'file' : f"{shortcuts.manta_variants_dir}somaticSV_PASS.vcf.complete" }
            misc.run_command(cmd_filter)

            elapsed = timeit.default_timer() - start
            misc.log_to_file("INFO", f'Manta SNV calling succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
    
        except Exception as e:
            misc.log_to_file("ERROR", ".manta() in dna_seq_analysis.py:", e)
            sys.exit()
