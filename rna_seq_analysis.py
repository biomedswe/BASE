from os import listdir, getenv, sys
import subprocess
import multiprocessing
import time
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
    def map_reads(self, options, filename, misc, shortcuts):
        '''This function map reads to the reference genome'''

        try:

            if misc.step_completed(f'{shortcuts.star_output_dir}{options.tumor_id}_{filename}_map.complete', 'Map reads to genome allready completed, skips step...'):
                time.sleep(2.5)
                pass
            else:
                reads = []

                for read in listdir(shortcuts.rna_reads_dir):
                    if options.tumor_id in read:
                        reads.append(read)
                    else:
                        print('Rna reads are incorrectly named')
                threads = multiprocessing.cpu_count() - 2
                print('Map reads to the genome\n\n')


                cmd_mapReads = f'''
                STAR --genomeDir {shortcuts.star_index_dir}{filename}_hg38_index \\
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
                --outFileNamePrefix {shortcuts.star_output_dir}{options.tumor_id} \\
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
                --outSAMattrRGline ID:{reads[0][:-16]} SM:{options.tumor_id} LB:{reads[0][:-16]} PL:"ILLUMINA" PU:{reads[0][:-16]} \\
                --twopassMode Basic'''
                misc.run_command(cmd_mapReads, 'Mapping reads to genome completed')
                cmd_sortsam = f"picard SortSam -I {shortcuts.star_output_dir}{options.tumor_id}Aligned.out.bam -O {shortcuts.star_output_dir}{options.tumor_id}Aligned_sorted.out.bam -SO coordinate"
                misc.run_command(cmd_sortsam, None)
                cmd_samtools_index = f"samtools index {shortcuts.star_output_dir}{options.tumor_id}Aligned_sorted.out.bam"
                misc.run_command(cmd_samtools_index, None)
                misc.create_trackFile(f'{shortcuts.star_output_dir}{options.tumor_id}_{filename}_map.complete')
        except Exception as e:
            print(f'Error with map_reads(): {e}')

    #---------------------------------------------------------------------------
    def ASEReadCounter(self, options, filename, misc, shortcuts):

        try:

            ref_dir = shortcuts.reference_genome_dir


            if misc.step_completed(f'{shortcuts.star_output_dir}{options.tumor_id}_{filename}_ase.complete', 'ASEReadCounter allready completed, skips step...'):
                time.sleep(2.5)
                pass
            else:
                cmd_ase = f'''
                gatk ASEReadCounter -R {ref_dir}{filename}_index/{filename}.fa \\
                --min-mapping-quality 10 --min-depth-of-non-filtered-base 10 \\
                --min-base-quality 2 \\
                --disable-read-filter NotDuplicateReadFilter \\
                --variant {shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het.vcf \\
                -I {shortcuts.star_output_dir}{options.tumor_id}Aligned_sorted.out.bam \\
                --output-format CSV \\
                --output {shortcuts.star_output_dir}{options.tumor_id}_{filename}_STAR_ASE.csv'''
                misc.run_command(cmd_ase, None)
                misc.create_trackFile(f'{shortcuts.star_output_dir}{options.tumor_id}_{filename}_ase.complete')
        except Exception as e:
            print(f'Error with ASEReadCounter(): {e}')

    #---------------------------------------------------------------------------
    def add_wgs_data_to_csv(self, options, filename, misc, shortcuts):
        '''This functions reads the CHROM and POS column of the vcf file and the allele depth 'AD' column. CHROM and POS columns are concatenated
        into a string as key in the dict and the 'AD' matching the CHROM and POS is the value. The csv file created by gatk ASEReadCounter is opened
        and 4 new columns with values are written to the file. df.apply applies the lambda function in every row in the csv file.'''

        try:

            dict = {}
            vcf_reader = vcfpy.Reader.from_path(f'{shortcuts.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf', 'r')

            for record in vcf_reader:
                chrom = record.CHROM
                pos = record.POS
                ad = record.calls[0].data.get('AD')
                geneName = record.__dict__['INFO']['ANN'][0].split('|')[3]
                variantType = record.__dict__['INFO']['ANN'][0].split('|')[1]
                dict[f"{chrom}-{pos}"] = [ad, geneName, variantType]



            df = pd.read_csv(f'{shortcuts.star_output_dir}{options.tumor_id}_chr13_STAR_ASE.csv')
            print(df)
            time.sleep(10)
            df['Dna_refCount'] = df.apply(lambda row: self.ad_tumor_coverage(row, 0, dict), axis=1)
            df['Dna_altCount'] = df.apply(lambda row: self.ad_tumor_coverage(row, 1, dict), axis=1)
            df['Dna_totalCount'] = df.apply(lambda row: self.ad_tumor_coverage(row, 'both', dict), axis=1)
            df['Dna_altAlleleFreq'] = df.apply(lambda row: self.alt_allele_freq(row, dict), axis=1)
            df['pValue_dna_altAlleleFreq'] = df.apply(lambda row: self.binominal_test(row), axis=1)
            df['geneName'] = df.apply(lambda row: self.add_gene_name(row, dict), axis=1)
            df['variantType'] = df.apply(lambda row: self.add_variant_type(row, dict), axis=1)
            # df.drop(['variantID', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs'], inplace=True)
            print(df)
            # df.to_csv('rna_seq/star/2064-01_chr13_STAR_ASE.csv', sep=',', index=False)
            print('RNA analysis completed, exiting program...')
            sys.exit()
        except Exception as e:
            print(f'Error with add_wgs_data_to_csv(): {e}')

        #---------------------------------------------------------------------------
    def add_gene_name(self, row, dict):

        try:
            chr_pos = (str(row[0]) + "-" + str(row[1]))
            if chr_pos in list(dict.keys()):
                return dict[chr_pos][1]

        except Exception as e:
            print(f'Error with add_gene_name: {e}')
            return None

    #---------------------------------------------------------------------------
    def add_variant_type(self, row, dict):

        try:
            chr_pos = (str(row[0]) + "-" + str(row[1]))
            if chr_pos in list(dict.keys()):
                return dict[chr_pos][2]

        except Exception as e:
            print(f'Error with add_variant_type: {e}')
            return None
    #---------------------------------------------------------------------------
    def ad_tumor_coverage(self, row, col, dict):
        '''This function looks for the string (CHROM-POS) in the dict and returns the matching value (list with two values eg. [8,8])'''

        try:
            chr_pos = (str(row[0]) + "-" + str(row[1]))
            if chr_pos in list(dict.keys()):
                return (dict[chr_pos][0][0] + dict[chr_pos][0][1]) if col == 'both' else dict[chr_pos][0][col]

        except Exception as e:
            print(f'Error with ad_tumor_coverage: {e}')
            return None

    #---------------------------------------------------------------------------
    def alt_allele_freq(self, row, dict):
        '''This function uses the chrom-pos keys to retrieve altCount and totalCount values from the dict. Then divides altCount with totalCount
        and returns the alt allele frequency'''

        try:
            chr_pos = (str(row[0]) + "-" + str(row[1]))
            altCount = dict[chr_pos][0][1]
            totalCount = dict[chr_pos][0][0] + dict[chr_pos][0][1]

            if chr_pos in list(dict.keys()):
                alt_allele_freq = altCount/totalCount
                return '{0:.3g}'.format(alt_allele_freq)

        except Exception as e:
            print(f'Error with alt_allele_freq: {e}')
            return None

    #---------------------------------------------------------------------------
    def binominal_test(self, row):
        '''This function reads the columns altCount, totalCount and Dna_altAlleleFreq from current row. The runs binom_test with these inputs
        and returns the answer'''

        try:
            altCount = row[6]
            totalCount = row[7]
            Dna_altAlleleFreq = row[16]
            avg_pvalue = binom_test(altCount, totalCount, float(Dna_altAlleleFreq))
            return avg_pvalue
        except Exception as e:
            print(f'Error with binominal_test function: {e}')
            return None

    #---------------------------------------------------------------------------
