from os import getenv, sys

class Shortcuts():
    '''This class contains shortcuts to key files and folders neccessary for the program'''

    def __init__(self, options):

        # Shortcut to main folders
        self.sequencing_project_dir = getenv("HOME")+"/sequencing_project/"
        self.dna_seq_dir = f"{self.sequencing_project_dir}dna_seq/"
        self.rna_seq_dir =  f"{self.sequencing_project_dir}rna_seq/"

        # Shortcuts to input folders
        self.dna_reads_dir  = f"{self.dna_seq_dir}reads/{options.tumor_id}/"
        self.reference_genome_dir = f"{self.sequencing_project_dir}reference_genome/"
        self.reference_genome_chunks_dir = f"{self.reference_genome_dir}chunks/"

        # Shortcuts to output folders in DNA sequencing analysis
        self.aligned_output_dir = f"{self.dna_seq_dir}aligned/"
        self.sorted_output_dir = f"{self.dna_seq_dir}sorted/"
        self.merged_output_dir = f"{self.dna_seq_dir}merged/"
        self.removed_duplicates_output_dir = f"{self.dna_seq_dir}removed_duplicates/"
        self.realigned_output_dir = f"{self.dna_seq_dir}realigned/"
        self.haplotypecaller_output_dir = f"{self.dna_seq_dir}gatk_haplotypecaller/"
        self.delly_output_dir = f"{self.dna_seq_dir}delly/"
        self.manta_output_dir = f"{self.dna_seq_dir}manta/"
        self.manta_variants_dir = f"{self.dna_seq_dir}manta/{options.tumor_id}/results/variants/"

        # Shortcuts to files used in DNA sequencing analysis
        self.reference_genome_file = f"{self.reference_genome_dir}GRCh38.p13.genome.fa"
        self.reference_genome_exclude_template_file = f"{self.sequencing_project_dir}excludeTemplate/human.hg38.excl.tsv"
        self.configManta_file = getenv("HOME")+"/anaconda3/envs/sequencing/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py"
        self.runWorkflow_file = getenv("HOME")+f"/sequencing_project/dna_seq/manta/{options.tumor_id}/runWorkflow.py"

        # Shortcuts to folders used in RNA sequencing analysis
        self.rna_reads_dir  = f"{self.rna_seq_dir}reads/{options.tumor_id}/"
        self.star_output_dir = f"{self.rna_seq_dir}star/{options.tumor_id}/"
        self.star_index_dir =  f"{self.reference_genome_dir}star_index/{options.tumor_id}/"



        # Shortcuts to files used in RNA sequencing analysis
        self.annotation_gtf_file = f"{self.reference_genome_dir}gencode.v37.primary_assembly.annotation.gtf"
        self.gatk_vcfFile = f"{self.haplotypecaller_output_dir}{options.tumor_id}/{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf"

        # Shortcuts to output lists (used for input in pipeline steps) (and also for validation if pipeline step i allready completed)
        self.alignedFiles_list = f"{self.aligned_output_dir}{options.tumor_id}/alignedFiles.txt"
        self.sortedFiles_list = f"{self.sorted_output_dir}{options.tumor_id}/sortedFiles.txt"
        self.mergedFiles_list = f"{self.merged_output_dir}{options.tumor_id}/mergedFiles.txt"
        self.removeDuplicates_list = f"{self.removed_duplicates_output_dir}{options.tumor_id}/remove_duplicate.txt"
        self.realignedFiles_list = f"{self.realigned_output_dir}{options.tumor_id}/realignedFiles.txt"
        self.gatk_chunks_list = f"{self.haplotypecaller_output_dir}{options.tumor_id}/chunks.list"

        # Shortcuts to files used to validate if pipeline step is allready completed
        self.bwa_index_complete = f"{self.reference_genome_dir}index.complete"
        self.haplotypecaller_complete = f"{self.haplotypecaller_output_dir}{options.tumor_id}/haplotypeCaller.complete"
        self.delly_complete = f"{self.delly_output_dir}{options.tumor_id}/delly.complete"
        self.manta_complete = f"{self.manta_output_dir}{options.tumor_id}/manta.complete"
