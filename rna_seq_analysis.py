from os import listdir, getenv, sys, path
import subprocess
import multiprocessing
import time
import timeit
from miscellaneous import Misc
try:
    import vcfpy
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
    from scipy.stats import binom_test
except Exception as e:
    Misc().log_to_file(f"importing in rna_seq_analysis.py: {e}")


class RnaSeqAnalysis():

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def index_genome_rna(self, misc, shortcuts):
        '''This function indexes either the whole genome or the chromosomes entered'''

        try:
            ref_dir = shortcuts.reference_genome_dir
            if not misc.step_allready_completed(f"{shortcuts.star_index_dir}starIndex.complete", "Indexing genome with STAR genomeGenerate"):
                start = timeit.default_timer()
                threads = multiprocessing.cpu_count() - 2
                misc.log_to_file(f"Starting: indexing genome with STAR using {threads} out of {threads+2} available threads")
                cmd_StarIndex = f'''
                STAR --runThreadN {threads} \\
                --runMode genomeGenerate \\
                --genomeDir {shortcuts.star_index_dir} \\
                --genomeFastaFiles {shortcuts.reference_genome_file} \\
                --sjdbGTFfile {shortcuts.annotation_gtf_file}'''
                if misc.run_command(cmd_StarIndex, None, None, f"{shortcuts.star_index_dir}starIndex.complete"):
                    elapsed = timeit.default_timer() - start
                    misc.log_to_file(f'Indexing whole genome with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
                    input('press any key to exit')
                    sys.exit()

        except Exception as e:
            misc.log_exception(".index_genome_rna in rna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def map_reads(self, options, misc, shortcuts):
        '''This function map reads to the reference genome'''



        try:
            if not misc.step_allready_completed(f"{shortcuts.star_output_dir}{options.tumor_id}/map.complete", f'Map reads to {options.tumor_id}'):
                reads = []
                for read in listdir(shortcuts.rna_reads_dir):
                    if options.tumor_id in read:
                        reads.append(read)
                    else:
                        misc.log_to_file('Rna reads are incorrectly named')
                        sys.exit()
                threads = multiprocessing.cpu_count() - 2
                misc.log_to_file(f'Starting: mapping reads ({options.tumor_id}) to genome with STAR using {threads} out of {threads+2} available threads')

                cmd_mapReads = f'''
                STAR --genomeDir {shortcuts.reference_genome_dir}/star_index/ \\
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
                --outFileNamePrefix {shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_ \\
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
                misc.run_command(cmd_mapReads, f'Mapping reads to genome', f'{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned.out.bam', None)
                cmd_sortsam = f"picard SortSam -I {shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned.out.bam -O {shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned_sorted.out.bam -SO coordinate"
                misc.run_command(cmd_sortsam, "Sorting BAM with Picard SortSam", f"{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned_sorted.out.bam", None)
                cmd_samtools_index = f"samtools index {shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned_sorted.out.bam"
                misc.run_command(cmd_samtools_index, "Indexing BAM with Samtools index", f"{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned_sorted.out.bam.bai",  f"{shortcuts.star_output_dir}{options.tumor_id}/map.complete")
                # add star wasp
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'All steps in mapping reads to gemome with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".map_reads() in rna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def ASEReadCounter(self, options, misc, shortcuts):

        try:
            if not misc.step_allready_completed(f"{shortcuts.star_output_dir}{options.tumor_id}/ase.complete", f'ASEReadCounter for {options.tumor_id}'):
                ref_dir = shortcuts.reference_genome_dir
                misc.log_to_file("Starting: Counting ASE reads using ASEReadCounter")
                start = timeit.default_timer()
                cmd_ase = f'''
                gatk ASEReadCounter -R {shortcuts.reference_genome_file} \\
                --min-mapping-quality 10 --min-depth-of-non-filtered-base 10 \\
                --min-base-quality 2 \\
                --disable-read-filter NotDuplicateReadFilter \\
                --variant {shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf \\
                -I {shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_Aligned_sorted.out.bam \\
                --output-format CSV \\
                --output {shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_STAR_ASE.csv'''
                misc.run_command(cmd_ase, None, f'{shortcuts.star_output_dir}{options.tumor_id}/ase.complete', f'{shortcuts.star_output_dir}{options.tumor_id}/ase.complete')
                elapsed = timeit.default_timer() - start
                misc.log_to_file(f'ASEReadCounter for {options.tumor_id} with STAR succesfully completed in {misc.elapsed_time(elapsed)} - OK!')
        except Exception as e:
            misc.log_exception(".ASEReadCounter() in rna_seq_analysis.py:", e)

    #---------------------------------------------------------------------------
    def add_wgs_data_to_csv(self, options, misc, shortcuts):
        '''This functions reads the CHROM and POS column of the vcf file and the allele depth 'AD' column. CHROM and POS columns are concatenated
        into a string as key in the dict and the 'AD' matching the CHROM and POS is the value. The csv file created by gatk ASEReadCounter is opened
        and 4 new columns with values are written to the file. df.apply applies the lambda function in every row in the csv file.'''

        try:

            start = timeit.default_timer()
            misc.log_to_file(f"Starting: Creating CSV...")
            dict = {}
            vcf_reader = vcfpy.Reader.from_path(f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf', 'r')
            variants_to_exclude = ['downstream_gene_variant', 'intergenic_region', 'intragenic_variant', 'intron_variant', 'splice_region_variant', 'splice_region_variant&intron_variant', 'upstream_gene_variant']

            for i, record in enumerate(vcf_reader):
                contig = record.CHROM
                position = record.POS
                refCount = record.calls[0].data.get('AD')[0]
                altCount = record.calls[0].data.get('AD')[1]
                totalCount = refCount + altCount
                geneName = record.__dict__['INFO']['ANN'][0].split('|')[3]
                variantType = record.__dict__['INFO']['ANN'][0].split('|')[1] if not record.__dict__['INFO']['ANN'][0].split('|')[1] in variants_to_exclude else None
                dict[i] = [contig, position, refCount, altCount, totalCount, geneName, variantType]

            # Make DataFrame out of dict
            df_vcf = pd.DataFrame.from_dict(dict, orient="index", columns=["contig", "position", "DNA_refCount", "DNA_altCount", "DNA_totalCount", "geneName", "variantType"])
            print(df_vcf)
            # read exel file with cnv information
            if not path.isfile(f'{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_CN.xlsx'):
                misc.log_to_file(f"No file specifying copynumber, save {options.tumor_id}_CN.xlsx in {shortcuts.star_output_dir}{options.tumor_id}/")
                sys.exit()
            df_cn = pd.read_excel(f'{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_CN.xlsx')

            # read csv from star output
            df = pd.read_csv(f'{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_STAR_ASE.csv')


            df["clinicalID"] = options.tumor_id.replace('-', '_')
            df_merge = pd.merge(df_vcf, df, on=["contig", "position"])
            df_merge.rename(columns={'refCount': 'RNA_refCount', 'altCount': 'RNA_altCount', 'totalCount': 'RNA_totalCount'}, inplace=True)
            df_merge = df_merge[['clinicalID', 'geneName', 'variantType', 'contig', 'position', 'variantID', 'refAllele', 'altAllele', 'RNA_refCount', 'RNA_altCount', 'RNA_totalCount', 'DNA_refCount', 'DNA_altCount', 'DNA_totalCount']]
            df_merge['RNA_VAF'] = df_merge.apply(lambda row: f"{row[9]/row[10]:.3f}", axis=1) #row[14]
            df_merge['DNA_VAF'] = df_merge.apply(lambda row: f"{row[12]/row[13]:.3f}", axis=1) #row[15]
            df_merge['pValue_WGS'] = df_merge.apply(lambda row: f"{binom_test(row[9], row[10], float(row[15])):.3f}", axis=1) #input: RNA_altCount, RNA_totalCount, DNA_VAF
            df_merge['VAF_ratio_WGS'] = df_merge.apply(lambda row: f"{float(row[14])/float(row[15]):.3f}", axis=1) # RNA_VAF/DNA_VAF
            df_merge['CN'] = df_merge.apply(lambda row: self.add_CNV(misc, df_cn, row), axis=1) # row[18]
            df_merge.dropna(inplace=True)
            df_merge['pValue_CNV'] = df_merge.apply(lambda row: self.calculate_pValue_CNV_data(misc, row), axis=1) #input: RNA_altCount, RNA_totalCount, CN
            df_merge['VAF_ratio_CNV'] = df_merge.apply(lambda row: self.caluculate_VAF_ratio_CNV_data(misc, row), axis=1) # RNA_VAF/CN
            df_merge.to_csv(f'{shortcuts.star_output_dir}{options.tumor_id}/{options.tumor_id}_STAR_ASE_completed.csv', sep=',', index=False)
            elapsed = timeit.default_timer() - start
            misc.log_to_file(f'Creating CSV completed in {misc.elapsed_time(elapsed)} - OK!')
            sys.exit()
        except Exception as e:
            misc.log_exception(".add_wgs_data_to_csv() in rna_seq_analysis.py:", e)

    def add_CNV(self, misc, df_cn, row):
        '''This function fetch the CN data from the CN.xlsx document for the specific position'''
        try:
            for i in df_cn.index:
                if row[3] == df_cn["Chromosome"].iloc[i] and df_cn["Start"].iloc[i] <= row[4] < df_cn["End"].iloc[i]:
                    return df_cn["Cn"].iloc[i]
        except Exception as e:
            misc.log_exception(".add_CNV() in rna_seq_analysis.py:", e)

    def caluculate_VAF_ratio_CNV_data(self, misc, row):
        '''This function determines what CN value that should be used when calculating the VAF_ratio_CNV'''
        try:

            if row[18] <= 2:
                return f"{float(row[14])/float(1/row[18]):.3f}" #1/2

            elif row[18] == 3:
                # if refCount < altCount
                if row[8] < row[9]:
                    return f"{float(row[14])/float(2/row[18]):.3f}" # 2/3
                # if refCount > altCount
                if row[8] > row[9]:
                    return f"{float(row[14])/float(1/row[18]):.3f}" #1/3

            elif row[18] == 4:
                if row[8] < row[9]:
                    return f"{float(row[14])/float(3/row[18]):.3f}" # 3/4
                elif row[8] > row[9]:
                    return f"{float(row[14])/float(1/row[18]):.3f}" # 1/4
                else:
                    return f"{float(row[14])/float(2/row[18]):.3f}" # 2/4
        except Exception as e:
            misc.log_exception(".caluculate_alt_chromosomes() in rna_seq_analysis.py:", e)

    def calculate_pValue_CNV_data(self, misc, row):
        '''This function determines what CN value that should be used when calculating the pValue_CNV'''
        try:
            if row[18] <= 2:
                return f"{binom_test(row[9], row[10], float(1/row[18])):.3f}" # 1/1, 1/2

            if row[18] == 3:
                # if refCount < altCount
                if row[8] < row[9]:
                    return f"{binom_test(row[9], row[10], float(2/row[18])):.3f}" # 2/3
                # if refCount > altCount
                if row[8] > row[9]:
                    return f"{binom_test(row[9], row[10], float(1/row[18])):.3f}" #1/3

            if row[18] == 4:
                if row[8] < row[9]:
                    return f"{binom_test(row[9], row[10], float(3/row[18])):.3f}" # 3/4
                elif row[8] > row[9]:
                    return f"{binom_test(row[9], row[10], float(1/row[18])):.3f}" # 1/4
                else:
                    return f"{binom_test(row[9], row[10], float(2/row[18])):.3f}" # 2/4
        except Exception as e:
            misc.log_exception(".caluculate_alt_chromosomes() in rna_seq_analysis.py:", e)
