import os
import configparser

class Shortcuts:
    def __init__(self, config_file="config.ini"):
        config = configparser.ConfigParser()
        config.read(config_file)

        self.gatk_path = config.get("third_party_programs", "gatk_path")
        self.bwa_path = config.get("third_party_programs", "bwa_path")
        self.samtools_path = config.get("third_party_programs", "samtools_path")
        self.sambamba_path = config.get("third_party_programs", "sambamba_path")
        self.STAR_path = config.get("third_party_programs", "STAR_path")
        self.bedtools_path = config.get("third_party_programs", "bedtools_path")
        self.bcftools_path = config.get("third_party_programs", "bcftools_path")
        self.tabix_path = config.get("third_party_programs", "tabix_path")
        self.Rscript_path = config.get("third_party_programs", "Rscript_path")
        self.bgzip_path = config.get("third_party_programs", "bgzip_path")
        
        # Shortcut to main folders
        self.BASE_dir = os.getenv("HOME") + "/BASE/"
        self.Rcode_dir = f"{self.BASE_dir}Rcode/"
        self.reference_genome_dir = f"{self.BASE_dir}reference_genome/"
        self.reference_cnv_dir = f"{self.BASE_dir}reference_genome/cnv_region/"
        self.star_index_dir = f"{self.BASE_dir}star_index/"
        self.tmp = f"{self.BASE_dir}tmp/"
        
        # Running data folder structure
        self.aligned_output_dir = f"{self.BASE_dir}dna_aligned/"
        self.rna_seq_dir = f"{self.BASE_dir}rna_aligned/"
        
        # Reference hg38.analysisSet folder
        self.reference_genome_file = f"{self.reference_genome_dir}hg38.analysisSet.fa"
        self.annotation_gtf_file = f"{self.reference_genome_dir}gencode.v43.basic.annotation.gtf"
        self.annotation_exon_file = f"{self.reference_genome_dir}gencode.v43.basic.exon.bed"
        
        # CNV calling training set folder
        self.reference_genome_cnv_region = f"{self.reference_cnv_dir}cnv_ref.bed"
        self.reference_cnv_panel = f"{self.reference_cnv_dir}cnv_ref_panel_hg38c.npz"
        self.reference_cnv_trainingset = f"{self.reference_cnv_dir}MeanShiftClustering.joblib"
        self.reference_cnv_trainingset_raw = f"{self.reference_cnv_dir}hg38_cnv.trainset_raw.gz"
        
        # R script for CNV calling
        self.CNV_calling_path = f"{self.Rcode_dir}cov2seg.r"
