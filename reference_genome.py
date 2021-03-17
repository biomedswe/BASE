import time, multiprocessing

class ReferenceGenome():

    def __init__(self):
        pass

    def download(self, misc, shortcuts):
        '''This function downloads the human reference genome GRCh38.fa and the comprehensive gene annotations gencode.v37.primary_assembly.annotation.gtf
        from https://www.gencodegenes.org/human/'''

        try:
            if misc.step_allready_completed(shortcuts.reference_genome_file) and misc.step_allready_completed(shortcuts.annotation_gtf_file):
                misc.logfile('Reference genome and annotation file allready downloaded.')
            else:
                misc.clear_screen()
                print("\033[1mDownload reference genome\033[0m\n\n")
                print("Downloading:")
                print("GRCh38.p13.genome.fa and from Comprehensive gene annotation https://www.gencodegenes.org/human/...\n")

                cmd_fasta_download = f"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz -P {shortcuts.GRCh38_dir}"
                misc.run_command(cmd_fasta_download, 'Download of GRCh38.p13.genome.fa.gz completed')

                cmd_fasta_unzip = f"gunzip {shortcuts.GRCh38_dir}GRCh38.p13.genome.fa.gz"
                misc.run_command(cmd_fasta_unzip, 'Unzip of GRCh38.p13.genome.fa.gz completed')

                cmd_gtf_download = f"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz -P {shortcuts.GRCh38_dir}"
                misc.run_command(cmd_gtf_download, 'Download of gencode.v37.primary_assembly.annotation.gtf.gz completed')

                cmd_gtf_unzip = f"gunzip {shortcuts.GRCh38_dir}gencode.v37.primary_assembly.annotation.gtf.gz"
                misc.run_command(cmd_gtf_unzip, 'Unzip of gencode.v37.primary_assembly.annotation.gtf.gz completed')

                print("Download completed!\nGRCh38.p13.genome.fa and gencode.v37.primary_assembly.annotation.gtf are saved in the reference_genome/bwa_index/GRCh38 folder.\n")
                return input("Press any key to return to main menu...")
        except Exception as e:
            print(f'Error with download(): {e}')
            input("Press any key to continue...")

    
