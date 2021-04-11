from os import getenv, sys

class Shortcuts():
    '''This class contains shortcuts to key files and folders neccessary for the program'''
    try:
        def __init__(self, options):

            # Shortcut to main folders
            self.sequencing_project_dir = getenv("HOME")+"/sequencing_project/"
            self.dna_seq_dir = f"{self.sequencing_project_dir}dna_seq/"
            self.rna_seq_dir =  f"{self.sequencing_project_dir}rna_seq/"

            # Shortcuts to input folders
            self.dna_reads_dir  = f"{self.dna_seq_dir}reads/"
            self.reference_genome_dir = f"{self.sequencing_project_dir}reference_genome/"
            self.GRCh38_dir = f"{self.reference_genome_dir}GRCh38.p13.genome/"
            self.GRCh38_chunks_dir = f"{self.GRCh38_dir}chunks/"

            # Shortcuts to output folders in DNA sequencing analysis
            self.aligned_output_dir = f"{self.dna_seq_dir}aligned/"
            self.sorted_output_dir = f"{self.dna_seq_dir}sorted/"
            self.merged_output_dir = f"{self.dna_seq_dir}merged/"
            self.removed_duplicates_output_dir = f"{self.dna_seq_dir}removed_duplicates/"
            self.realigned_output_dir = f"{self.dna_seq_dir}realigned/"
            self.haplotypecaller_output_dir = f"{self.dna_seq_dir}gatk_haplotypecaller/"
            self.haplotypecaller_chunks_dir = f"{self.haplotypecaller_output_dir}chunks/"
            self.delly_output_dir = f"{self.dna_seq_dir}delly/"
            self.manta_output_dir = f"{self.dna_seq_dir}manta/"
            self.manta_variants_dir = f"{self.dna_seq_dir}manta/results/variants/"

            # Shortcuts to files used in DNA sequencing analysis
            self.reference_genome_file = f"{self.GRCh38_dir}GRCh38.p13.genome.fa"
            self.reference_genome_exclude_template_file = f"{self.sequencing_project_dir}excludeTemplate/human.hg38.excl.tsv"
            self.configManta_file = getenv("HOME")+"/anaconda3/envs/sequencing/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py"
            self.runWorkflow_file = getenv("HOME")+"/sequencing_project/dna_seq/manta/runWorkflow.py"

            # Shortcuts to folders used in RNA sequencing analysis
            self.rna_reads_dir  = f"{self.rna_seq_dir}reads/"
            self.star_output_dir = f"{self.rna_seq_dir}star/"
            self.star_index_dir_whole_genome =  f"{self.GRCh38_dir}star_index/"



            # Shortcuts to files used in RNA sequencing analysis
            self.annotation_gtf_file = f"{self.GRCh38_dir}gencode.v37.primary_assembly.annotation.gtf"
            self.gatk_vcfFile = f"{self.haplotypecaller_output_dir}{options.tumor_id}_filtered_RD10_snps_tumor_het_annotated.vcf"

            # Shortcuts to output lists (used for input in pipeline steps) (and also for validation if pipeline step i allready completed)
            self.alignedFiles_list = f"{self.aligned_output_dir}alignedFiles.txt"
            self.sortedFiles_list = f"{self.sorted_output_dir}sortedFiles.txt"
            self.mergedFiles_list = f"{self.merged_output_dir}mergedFiles.txt"
            self.removeDuplicates_list = f"{self.removed_duplicates_output_dir}remove_duplicate.txt"
            self.realignedFiles_list = f"{self.realigned_output_dir}realignedFiles.txt"
            self.gatk_chunks_list = f"{self.haplotypecaller_chunks_dir}chunks.list"

            # Shortcuts to files used to validate if pipeline step is allready completed
            self.anaconda_setup_complete = getenv("HOME")+'/anaconda3/install.complete'
            self.bwa_index_whole_reference_genome_complete = f"{self.GRCh38_dir}index.complete"
            self.validate_bam_complete = f"{self.aligned_output_dir}validateBam.complete"
            self.haplotypecaller_complete = f"{self.haplotypecaller_output_dir}haplotypeCaller.complete"
            self.delly_complete = f"{self.delly_output_dir}delly.complete"
            self.manta_complete = f"{self.manta_output_dir}manta.complete"

    except Exception as e:
        self.log_to_file(f'Error with Shortcuts.__init__ in shortcuts.py: {e}. Exiting program...')
        sys.exit()
