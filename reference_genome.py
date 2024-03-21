import os
import re
import shutil
import gzip
from Misc2 import Misc
from shortcuts import Shortcuts

class ReferenceGenome:
    def __init__(self):
        self.shortcuts = Shortcuts()
        self.misc = Misc(dry_run=True)

    def create_folder_structure(self):
        directories = [
            self.shortcuts.BASE_dir,
            self.shortcuts.Rcode_dir,
            self.shortcuts.reference_genome_dir,
            self.shortcuts.reference_cnv_dir,
            self.shortcuts.star_index_dir,
            self.shortcuts.aligned_output_dir,
            self.shortcuts.rna_seq_dir,
            self.shortcuts.tmp
        ]

        for directory in directories:
            os.makedirs(directory, exist_ok=True)
            print(f"Directory created: {directory}")

            
    def download_genome_and_annotation(self):
        # Download genome and annotation files
        try:
            self.download_file(
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
                self.shortcuts.reference_genome_dir,
                "hg38.analysisSet.fa.gz",
                "Downloading hg38.analysisSet.fa from UCSC Genome Browser",
                self.shortcuts.reference_genome_file
            )

            self.download_file(
                "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz",
                self.shortcuts.reference_genome_dir,
                "gencode.v43.basic.annotation.gtf.gz",
                "Downloading gencode v43 basic annotation gtf.gz",
                self.shortcuts.annotation_gtf_file
            )

            self.index_genome_dna()

            self.misc.log_to_file("INFO", "Download, indexing, and GTF parsing completed successfully.")
        except Exception as e:
            self.misc.log_exception(self, e, "Error in download_genome_and_annotation")

    def download_file(self, url, save_dir, file_name, log_message, output_file):
        # Download and decompress files
        self.misc.log_to_file("INFO", log_message)
        download_cmd = f"wget {url} -P {save_dir}"
        self.misc.run_command(download_cmd)

        unzip_cmd = f"gunzip -c {os.path.join(save_dir, file_name)} > {output_file}"
        self.misc.run_command(unzip_cmd)

    def index_genome_dna(self):
        # Index the genome with samtools, bwa, STAR, and perform GTF parsing, sorting, merging, compressing, and indexing
        try:
            # Index with samtools
            self.misc.run_command(f"{self.shortcuts.samtools_path} faidx {self.shortcuts.reference_genome_file}")
            # Create a dictionary for the reference
            dict_file = os.path.splitext(self.shortcuts.reference_genome_file)[0] + ".dict"
            self.misc.run_command(f"{self.shortcuts.gatk_path} CreateSequenceDictionary -R {self.shortcuts.reference_genome_file} -O {dict_file}")
            # Index with bwa
            self.misc.log_to_file("INFO", "Indexing genome with bwa index...")
            self.misc.run_command(f"{self.shortcuts.bwa_path} index {self.shortcuts.reference_genome_file}")
            # STAR genome indexing
            self.misc.log_to_file("INFO", "Generating STAR genome index...")
            star_index_cmd = (
                f"{self.shortcuts.STAR_path} --runThreadN 8 --runMode genomeGenerate "
                f"--genomeDir {self.shortcuts.star_index_dir} "
                f"--genomeFastaFiles {self.shortcuts.reference_genome_file} "
                f"--sjdbGTFfile {self.shortcuts.annotation_gtf_file} "
                f"--sjdbOverhang 100"
            )
            self.misc.run_command(star_index_cmd)
            self.misc.log_to_file("INFO", "Genome indexing with STAR completed.")
            
            # Parse the GTF file
            self.parser_gtf()
            parsed_gtf_path = os.path.join(self.shortcuts.reference_genome_dir, "parsed_gtf.txt")
            sorted_gtf_path = os.path.join(self.shortcuts.reference_genome_dir, "sorted_parsed_gtf.bed")
            merged_gtf_path = self.shortcuts.annotation_exon_file
            self.misc.run_command(f"{self.shortcuts.bedtools_path} sort -i {parsed_gtf_path} > {sorted_gtf_path}")
            self.misc.run_command(f"{self.shortcuts.bedtools_path} merge -i {sorted_gtf_path} > {merged_gtf_path}")
            self.misc.run_command(f"{self.shortcuts.bgzip_path} -c {merged_gtf_path} > {merged_gtf_path}.gz")
            self.misc.run_command(f"{self.shortcuts.tabix_path} -p bed {merged_gtf_path}.gz")
            self.misc.log_to_file("INFO", "Genome indexing, GTF parsing, sorting, merging, compression, and indexing completed.")
            #os.remove(parsed_gtf_path)
            #os.remove(sorted_gtf_path)

        except Exception as e:
            self.misc.log_exception(self, e, "Error in index_genome_dna")

    def parser_gtf(self):
        # Parse GTF to extract relevant exon information
        gtf_path = self.shortcuts.annotation_gtf_file
        output_path = os.path.join(self.shortcuts.reference_genome_dir, "parsed_gtf.txt")
        
        with open(gtf_path, "r") as gtf_fh, open(output_path, "w") as result:
            while True:  # Skip header lines
                line = gtf_fh.readline()
                if not line.startswith('#'):
                    break
            for line in gtf_fh:
                record = line.strip().split("\t")
                if record[2] == "exon":
                    gtf_annot = str(record[8]).split(";")
                    MANE_flag = any(re.search('MANE_Select', item) for item in gtf_annot)
                    Egene_idx = next((i for i, item in enumerate(gtf_annot) if re.search('^gene_id', item)), None)
                    name_idx = next((i for i, item in enumerate(gtf_annot) if re.search('gene_name', item)), None)
                    if MANE_flag and Egene_idx is not None and name_idx is not None:
                        ENSID = re.findall(r'"([^"]*)"', gtf_annot[Egene_idx])[0]
                        GeneSymbol = re.findall(r'"([^"]*)"', gtf_annot[name_idx])[0]
                        result.write(f"{record[0]}\t{int(record[3])}\t{int(record[4])}\t{ENSID}\t{GeneSymbol}\t{record[6]}\n")
        self.misc.log_to_file("INFO", "GTF parsing completed.")


    def move_and_decompress_cnv_files(self):
        current_working_directory = os.getenv("PWD")

        Rsource_file_path = os.path.join(current_working_directory, 'cov2seg.r')
        destination_path = os.path.join(self.shortcuts.Rcode_dir, 'cov2seg.r')
        shutil.move(Rsource_file_path, destination_path)
        print(f"File 'cov2seg.r' copied from {Rsource_file_path} to {destination_path}")

        cnv_data_dir = os.path.join(current_working_directory, 'CNV_data')
        for filename in os.listdir(cnv_data_dir):
            src_file = os.path.join(cnv_data_dir, filename)
            dest_file = os.path.join(self.shortcuts.reference_cnv_dir, filename)
            shutil.move(src_file, dest_file)
            print(f"Moved {filename} to {self.shortcuts.reference_cnv_dir}")

        gz_file = os.path.join(self.shortcuts.reference_cnv_dir, 'cnv_ref.bed.gz')
        if os.path.exists(gz_file):
            with gzip.open(gz_file, 'rb') as f_in:
                with open(gz_file[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"Decompressed {gz_file}")


    
def main():
    reference_genome = ReferenceGenome()
    reference_genome.create_folder_structure()
    reference_genome.move_and_decompress_cnv_files()
    reference_genome.download_genome_and_annotation()
   

if __name__ == "__main__":
    main()