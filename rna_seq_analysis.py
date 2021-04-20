from os import listdir, getenv, sys
import subprocess
import multiprocessing
import time
import timeit
try:
    import vcfpy
    import pandas as pd
    from Bio import SeqIO
    from scipy.stats import binom
    from scipy.stats import binom_test
except Exception as e:
    print(f'{e}')

class RnaSeqAnalysis():

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def index_genome_rna(self, choice, filename, misc, shortcuts):
        '''This function indexes either the whole genome or the chromosomes entered'''

        try:

            ref_dir = shortcuts.reference_genome_dir

            # Index whole genome
            if choice == 1:
                if not misc.step_allready_completed(f"{shortcuts.star_index_dir_whole_genome}{filename}_starIndex.complete", "Indexing whole genom with STAR genomeGenerate"):
                    start = timeit.default_timer()
                    threads = multiprocessing.cpu_count() - 2
                    misc.log_to_file(f"Starting: indexing whole genome with STAR using {threads} out of {threads+2} available threads")
                    cmd_StarIndex = f'''
                    STAR --runThreadN {threads} \\
                    --runMode genomeGenerate \\
                    --genomeDir {shortcuts.star_index_dir_whole_genome} \\
                    --genomeFastaFiles {shortcuts.reference_genome_file} \\
                    --sjdbGTFfile {shortcuts.annotation_gtf_file}'''
                    if misc.run_command(cmd_StarIndex, None, None, f"{shortcuts.star_index_dir_whole_genome}{filename}_starIndex.complete"):
                        elapsed = timeit.default_timer() - start
                        misc.log_to_file(f'Indexing whole genome with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                        input('press any key to exit')
                        sys.exit()

            # Index parts of genome
            elif choice == 2:
                    if not misc.step_allready_completed(f'{ref_dir}{filename}/star_index/starIndex.complete', f'Genome indexing {filename} with star'):
                        start = timeit.default_timer()
                        threads = multiprocessing.cpu_count() - 2
                        misc.log_to_file(f"Starting: indexing {filename} with STAR using {threads} out of {threads+2} available threads")
                        cmd_StarIndex = f'''
                        STAR --runThreadN {threads} \\
                        --genomeSAindexNbases 12 \\
                        --runMode genomeGenerate \\
                        --genomeDir {ref_dir}{filename}/star_index \\
                        --genomeFastaFiles {ref_dir}{filename}/{filename}.fa \\
                        --sjdbGTFfile {ref_dir}{filename}/{filename}.gtf'''
                        if misc.run_command(cmd_StarIndex,  f'Genome indexing {filename} with star', None, f'{ref_dir}{filename}/star_index/starIndex.complete'):
                            elapsed = timeit.default_timer() - start
                            misc.log_to_file(f'Indexing {filename} with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                            input('press any key to return to previous menu...')

        except Exception as e:
            misc.log_exception(".index_genome_rna in rna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def map_reads(self, options, filename, misc, shortcuts):
        '''This function map reads to the reference genome'''



        try:
            if not filename:
                filename = input("You have restarted the program since you indexed parts of genome, please type in the chromosome/s you want to map to (eg. chr13):")
                filename += "_GRCh38.p13.genome"

            if not misc.step_allready_completed(f"{shortcuts.star_output_dir}{filename}/map.complete", f'Map reads to {filename} genome'):
                reads = []
                for read in listdir(shortcuts.rna_reads_dir):
                    if options.tumor_id in read:
                        reads.append(read)
                    else:
                        misc.log_to_file('Rna reads are incorrectly named')
                        sys.exit()
                threads = multiprocessing.cpu_count() - 2
                misc.log_to_file(f'Starting: mapping reads to {filename} genome with STAR using {threads} out of {threads+2} available threads')

                cmd_mapReads = f'''
                STAR --genomeDir {shortcuts.reference_genome_dir}{filename}/star_index/ \\
                --readFilesIn {shortcuts.rna_reads_dir}{reads[0]} {shortcuts.rna_reads_dir}{reads[1]} \\
                --runThreadN {threads} \\
                --alignIntronMax 1000000 \\
                --alignIntronMin 20 \\
                --alignMatesGapMax 1000000 \\
                --alignSJDBoverhangMin 1 \\
                --alignSJoverhangMin 8 \\
                --alignSoftClipAtReferenceEnds Yes \\
                --chimJunctionOverhangMin 15 \\
                --chimMainSegmentMultNmax 1 \\
                --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \\
                --chimSegmentMin 15 \\
                --genomeLoad NoSharedMemory \\
                --limitSjdbInsertNsj 1200000 \\
                --outFileNamePrefix {shortcuts.star_output_dir}{filename}/{options.tumor_id}_ \\
                --outFilterIntronMotifs None \\
                --outFilterMatchNminOverLread 0.33 \\
                --outFilterMismatchNmax 999 \\
                --outFilterMismatchNoverLmax 0.1 \\
                --outFilterMultimapNmax 20 \\
                --outFilterScoreMinOverLread 0.33 \\
                --outFilterType BySJout \\
                --outSAMattributes NH HI AS nM NM MD XS ch vA vG vW \\
                --outSAMstrandField intronMotif \\
                --outSAMtype BAM Unsorted \\
                --outSAMunmapped Within \\
                --quantMode TranscriptomeSAM GeneCounts \\
                --readFilesCommand zcat \\
                --waspOutputMode SAMtag \\
                --varVCFfile {shortcuts.gatk_vcfFile} \\
                --outSAMattrRGline ID:{reads[0][:16]} SM:{options.tumor_id} LB:{reads[0][:16]} PL:"ILLUMINA" PU:{reads[0][:16]} \\
                --twopassMode Basic'''
                start = timeit.default_timer()
                misc.run_command(cmd_mapReads, f'Mapping reads to {filename} genome', f'{shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned.out.bam', None)
                cmd_sortsam = f"picard SortSam -I {shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned.out.bam -O {shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned_sorted.out.bam -SO coordinate"
                misc.run_command(cmd_sortsam, "Sorting BAM with Picard SortSam", f"{shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned_sorted.out.bam", None)
                cmd_samtools_index = f"samtools index {shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned_sorted.out.bam"
                misc.run_command(cmd_samtools_index, "Indexing BAM with Samtools index", f"{shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned_sorted.out.bam.bai",  f"{shortcuts.star_output_dir}{filename}/map.complete")
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'All steps in mapping reads to {filename} with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".map_reads() in rna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def ASEReadCounter(self, options, filename, misc, shortcuts):

        try:
            if not misc.step_allready_completed(f"{shortcuts.star_output_dir}{filename}/ase.complete", f'ASEReadCounter for {filename}'):
                ref_dir = shortcuts.reference_genome_dir
                misc.log_to_file("Starting: Counting ASE reads using ASEReadCounter")
                start = timeit.default_timer()
                cmd_ase = f'''
                gatk ASEReadCounter -R {ref_dir}{filename}/{filename}.fa \\
                --min-mapping-quality 10 --min-depth-of-non-filtered-base 10 \\
                --min-base-quality 2 \\
                --disable-read-filter NotDuplicateReadFilter \\
                --variant {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf \\
                -I {shortcuts.star_output_dir}{filename}/{options.tumor_id}_Aligned_sorted.out.bam \\
                --output-format CSV \\
                --output {shortcuts.star_output_dir}{filename}/{options.tumor_id}_{filename}_STAR_ASE.csv'''
                misc.run_command(cmd_ase, None, f'{shortcuts.star_output_dir}{filename}/ase.complete', f'{shortcuts.star_output_dir}{filename}/ase.complete')
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'ASEReadCounter for {filename} with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".ASEReadCounter() in rna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def add_wgs_data_to_csv(self, options, filename, misc, shortcuts):
        '''This functions reads the CHROM and POS column of the vcf file and the allele depth 'AD' column. CHROM and POS columns are concatenated
        into a string as key in the dict and the 'AD' matching the CHROM and POS is the value. The csv file created by gatk ASEReadCounter is opened
        and 4 new columns with values are written to the file. df.apply applies the lambda function in every row in the csv file.'''

        try:
            start = timeit.default_timer()
            misc.log_to_file(f"Starting: Creating CSV...")
            dict = {}
            vcf_reader = vcfpy.Reader.from_path(f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf', 'r')

            for i, record in enumerate(vcf_reader):
                contig = record.CHROM
                position = record.POS
                refCount = record.calls[0].data.get('AD')[0]
                altCount = record.calls[0].data.get('AD')[1]
                totalCount = refCount + altCount
                geneName = record.__dict__['INFO']['ANN'][0].split('|')[3]
                variantType = record.__dict__['INFO']['ANN'][0].split('|')[1]
                dict[i] = [f"{contig}-{position}", contig, position, refCount, altCount, totalCount, geneName, variantType]

            df_cn = pd.read_excel(f'{shortcuts.star_output_dir}{filename}/3751-04_CN.xlsx')



            df_vcf = pd.DataFrame.from_dict(dict, orient="index", columns=["key", "contig", "position", "DNA_refCount", "DNA_altCount", "DNA_totalCount", "geneName", "variantType"])
            df = pd.read_csv(f'{shortcuts.star_output_dir}{filename}/{options.tumor_id}_{filename}_STAR_ASE.csv')

            df["clinicalID"] = options.tumor_id
            df_merge = pd.merge(df_vcf, df, on=["contig", "position"])
            df_merge.rename(columns={'refCount': 'RNA_refCount', 'altCount': 'RNA_altCount', 'totalCount': 'RNA_totalCount'}, inplace=True)
            pd.set_option('display.max_columns', 20)
            df_merge = df_merge[['clinicalID', 'contig', 'position', 'variantID', 'refAllele', 'altAllele', 'RNA_refCount', 'RNA_altCount', 'RNA_totalCount', 'DNA_refCount', 'DNA_altCount', 'DNA_totalCount']]
            df_merge['RNA_VAF'] = df_merge.apply(lambda row: f"{row[7]/row[8]:.3f}", axis=1) #row[12]
            df_merge['DNA_VAF'] = df_merge.apply(lambda row: f"{row[10]/row[11]:.3f}", axis=1) #row[13]
            df_merge['pValue_WGS_data'] = df_merge.apply(lambda row: f"{binom_test(row[7], row[8], float(row[13])):.3f}", axis=1) #input: RNA_altCount, RNA_totalCount, DNA_VAF
            df_merge['VAF_ratio_WGS_data'] = df_merge.apply(lambda row: f"{float(row[12])/float(row[13]):.3f}", axis=1) # RNA_VAF/DNA_VAF
            df_merge['CN'] = df_merge.apply(lambda row: self.add_CNV(misc, df_cn, row), axis=1) # row[16]
            df_merge['pValue_CNV_data'] = df_merge.apply(lambda row: f"{binom_test(row[7], row[8], float(1/row[16])):.3f}", axis=1) #input: altCount, totalCount,
            df_merge['VAF_ratio_CNV_data'] = df_merge.apply(lambda row: f"{float(row[12])/float(1/row[16]):.3f}", axis=1) # RNA_VAF/DNA_VAF
            df_merge.to_csv(f'{shortcuts.star_output_dir}{filename}/{options.tumor_id}_{filename}_STAR_ASE_completed.csv', sep=',', index=False)
            elapsed = timeit.default_timer() - start
            misc.log_to_file(f'Creating CSV completed in {misc.elapsed_time(elapsed)} - OK!')
            sys.exit()
        except Exception as e:
            misc.log_exception(".add_wgs_data_to_csv() in rna_seq_analysis.py:", e)

    def add_CNV(self, misc, df_cn, row):
        try:
            for i in df_cn.index:
                if row[1] == df_cn["Chromosome"].iloc[i] and df_cn["Start"].iloc[i] <= row[2] < df_cn["End"].iloc[i]:
                    return df_cn["Cn"].iloc[i]
        except Exception as e:
            misc.log_exception(".add_CNV() in rna_seq_analysis.py:", e)
