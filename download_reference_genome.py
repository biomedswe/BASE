import subprocess
from misc import *

def download_ref():

    while True:
        if confirm_choice():
            clear_screen()
            print("Download reference genome\n\n")
            print("Downloading:")
            print("GRCh38.p13.genome.fa from https://www.gencodegenes.org/human/...\n")
            print("Comprehensive gene annotation from https://www.gencodegenes.org/human/...\n")

            cmd_fasta_download = "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz -P $HOME/sequencing_project/reference_genome/"
            run_command(cmd__genome_download, 'Download of GRCh38.p13.genome.fa.gz completed')

            cmd_gtf_download = "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz -P $HOME/sequencing_project/reference_genome/"
            run_command(cmd__gtf_download, 'Download of gencode.v37.primary_assembly.annotation.gtf.gz completed')

            cmd_fasta_unzip = "gunzip $HOME/sequencing_project/reference_genome/GRCh38.p13.genome.fa.gz"
            run_command(cmd__fasta_unzip, 'Unzip of GRCh38.p13.genome.fa.gz completed')

            cmd_gtf_unzip = "gunzip $HOME/sequencing_project/reference_genome/gencode.v37.primary_assembly.annotation.gtf.gz"
            run_command(cmd__gtf_unzip, 'Unzip of gencode.v37.primary_assembly.annotation.gtf.gz completed')

            print("Completed!\nGRCh38.p13.genome.fa and gencode.v37.primary_assembly.annotation.gtf are saved in the reference_genome folder.\n")
            return input("Press any key to return to main menu...")
        else:
            break
