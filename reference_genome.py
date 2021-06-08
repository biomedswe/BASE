import time, multiprocessing
from os import sys

class ReferenceGenome():

    def __init__(self):
        pass

    def download(self, misc, shortcuts):
        '''This function downloads the human reference genome GRCh38.p13.genome.fa and the comprehensive gene annotations gencode.v37.primary_assembly.annotation.gtf
        from https://www.gencodegenes.org/human/'''

        try:
            misc.clear_screen()
            misc.log_to_file("info", "Downloading GRCh38.p13.genome.fa and gencode.v37.primary_assembly.annotation.gtf from https://www.gencodegenes.org/human/...")
            cmd_fasta_download = f"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz -P {shortcuts.reference_genome_dir}"
            if misc.run_command(cmd_fasta_download, 'Downloading GRCh38.p13.genome.fa.gz', shortcuts.reference_genome_file, None):
                misc.log_to_file("info", "Unzipping, please wait...")
                cmd_fasta_unzip = f"gunzip {shortcuts.reference_genome_dir}GRCh38.p13.genome.fa.gz"
                misc.run_command(cmd_fasta_unzip, 'Unzipping GRCh38.p13.genome.fa.gz', None, None)

            cmd_gtf_download = f"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz -P {shortcuts.reference_genome_dir}"
            if misc.run_command(cmd_gtf_download, 'Downloading gencode.v37.primary_assembly.annotation.gtf.gz', shortcuts.annotation_gtf_file, None):
                misc.log_to_file("info", "Unzipping, please wait...")
                cmd_gtf_unzip = f"gunzip {shortcuts.reference_genome_dir}gencode.v37.primary_assembly.annotation.gtf.gz"
                misc.run_command(cmd_gtf_unzip, 'Unzipping gencode.v37.primary_assembly.annotation.gtf.gz', None, None)
                misc.log_to_file("info", "Download completed!\nGRCh38.p13.genome.fa and gencode.v37.primary_assembly.annotation.gtf are saved in the reference_genome/GRCh38.p13.genome folder.\n")
            return input("Press any key to return to previous menu...")
        except Exception as e:
            misc.log_exception(self, e, ".download() in reference_genome.py:")
            input('press any key to exit')
            sys.exit()
