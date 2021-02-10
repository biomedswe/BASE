import subprocess
from misc import *

def download_ref():

    while True:
        if confirm_choice():

            print("\nDownloading GRCh38.p13.genome.fa from https://www.gencodegenes.org/human/...\n")

            cmd_download = "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz -P $HOME/sequencing_project/reference_genome/"
            # subprocess.run(cmd_download, shell=True)

            cmd_unzip = "gunzip $HOME/sequencing_project/reference_genome/GRCh38.p13.genome.fa.gz"
            # subprocess.run(cmd_unzip, shell=True)
            print("Completed!\nGRCh38.p13.genome.fa is saved in the reference_genome folder.\n")
            return input("Press any key to return to main menu...")
        else:
            break
