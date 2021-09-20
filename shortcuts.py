from os import getenv, sys

class Shortcuts():
    '''This class contains shortcuts to key files and folders neccessary for the program'''

    def __init__(self, options):

        # Shortcut to main folders
        self.BASE_dir = getenv("HOME")+"/BASE/"
        self.dna_seq_dir = f"{self.BASE_dir}dna_seq/"
        self.rna_seq_dir =  f"{self.BASE_dir}rna_seq/"

        # Shortcuts to input folders
        self.dna_reads_dir  = f"{self.dna_seq_dir}reads/{options.tumor_id}/"
        self.reference_genome_dir = f"{self.BASE_dir}reference_genome/"
        self.reference_genome_chunks_dir = f"{self.reference_genome_dir}chunks/"

        # Shortcuts to output folders in DNA sequencing analysis
        self.aligned_output_dir = f"{self.dna_seq_dir}aligned/{options.tumor_id}/"
        self.sorted_output_dir = f"{self.dna_seq_dir}sorted/{options.tumor_id}/"
        self.merged_output_dir = f"{self.dna_seq_dir}merged/{options.tumor_id}/"
        self.removed_duplicates_output_dir = f"{self.dna_seq_dir}removed_duplicates/{options.tumor_id}/"
        self.realigned_output_dir = f"{self.dna_seq_dir}realigned/{options.tumor_id}/"
        self.haplotypecaller_output_dir = f"{self.dna_seq_dir}gatk_haplotypecaller/{options.tumor_id}/"
        self.delly_output_dir = f"{self.dna_seq_dir}delly/{options.tumor_id}/"
        self.manta_output_dir = f"{self.dna_seq_dir}manta/{options.tumor_id}/"
        self.manta_variants_dir = f"{self.dna_seq_dir}manta/{options.tumor_id}/results/variants/"

        # Shortcuts to files used in DNA sequencing analysis
        self.reference_genome_file = f"{self.reference_genome_dir}GRCh38.p13.genome.fa"
        self.reference_genome_exclude_template_file = f"{self.BASE_dir}excludeTemplate/human.hg38.excl.tsv"
        self.configManta_file = getenv("HOME")+"/anaconda3/envs/sequencing/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py"
        self.runWorkflow_file = getenv("HOME")+f"/BASE/dna_seq/manta/{options.tumor_id}/runWorkflow.py"

        # Shortcuts to folders used in RNA sequencing analysis
        self.rna_reads_dir  = f"{self.rna_seq_dir}reads/{options.tumor_id}/"
        self.star_output_dir = f"{self.rna_seq_dir}star/{options.tumor_id}/"
        self.star_index_dir =  f"{self.reference_genome_dir}star_index/{options.tumor_id}/"



        # Shortcuts to files used in RNA sequencing analysis
        self.annotation_gtf_file = f"{self.reference_genome_dir}gencode.v37.primary_assembly.annotation.gtf"
        self.gatk_vcfFile = f"{self.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf"

        # Shortcuts to output lists (used for input in pipeline steps) (and also for validation if pipeline step i allready completed)  
        self.gatk_chunks_list = f"{self.haplotypecaller_output_dir}/chunks.list"

      