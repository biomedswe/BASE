import time, multiprocessing
from os import sys

class ReferenceGenome():

    def __init__(self):
        pass

    def download(self, misc, shortcuts):
        '''This function downloads the human reference genome GRCh38.fa and the comprehensive gene annotations gencode.v37.primary_assembly.annotation.gtf
        from https://www.gencodegenes.org/human/'''

        try:
            misc.clear_screen()
            misc.log_to_file("Downloading GRCh38.p13.genome.fa and Comprehensive gene annotation from https://www.gencodegenes.org/human/...")
            cmd_fasta_download = f"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz -P {shortcuts.GRCh38_dir}"
            if misc.run_command(cmd_fasta_download, 'Downloading GRCh38.p13.genome.fa.gz', shortcuts.reference_genome_file, None):
                cmd_fasta_unzip = f"gunzip {shortcuts.GRCh38_dir}GRCh38.p13.genome.fa.gz"
                misc.run_command(cmd_fasta_unzip, 'Unzipping GRCh38.p13.genome.fa.gz', None, None)

            cmd_gtf_download = f"wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz -P {shortcuts.GRCh38_dir}"
            if misc.run_command(cmd_gtf_download, 'Downloading gencode.v37.primary_assembly.annotation.gtf.gz', shortcuts.annotation_gtf_file, None):
                cmd_gtf_unzip = f"gunzip {shortcuts.GRCh38_dir}gencode.v37.primary_assembly.annotation.gtf.gz"
                misc.run_command(cmd_gtf_unzip, 'Unzipping gencode.v37.primary_assembly.annotation.gtf.gz', None, None)
                misc.log_to_file("Download completed!\nGRCh38.p13.genome.fa and gencode.v37.primary_assembly.annotation.gtf are saved in the reference_genome/bwa_index/GRCh38 folder.\n")
            return input("Press any key to return to previous menu...")
        except Exception as e:
            print(f'Error with ReferenceGenome.download() in reference_genome.py: {e}')
            input('press any key to exit')
            sys.exit()

    #---------------------------------------------------------------------------
    def index_genome_rna(self, choice, filename, misc, shortcuts):
        '''This function indexes either the whole genome or the chromosomes entered'''

        try:

            ref_dir = shortcuts.reference_genome_dir

            # Index whole genome
            if choice == 1:
                if misc.step_allready_completed(shortcuts.star_whole_genome_indexing_complete):
                    misc.log_to_file('Whole genome indexing allready completed, returning...')
                else:
                    threads = multiprocessing.cpu_count() - 2
                    cmd_StarIndex = f'''
                    STAR --runThreadN {threads} \\
                    --runMode genomeGenerate \\
                    --genomeDir {shortcuts.star_index_dir_whole_genome} \\
                    --genomeFastaFiles {shortcuts.reference_genome_file} \\
                    --sjdbGTFfile {shortcuts.annotation_gtf_file}'''
                    print(cmd_StarIndex)
                    misc.run_command(cmd_StarIndex, '\nIndexing whole genom with STAR genomeGenerate...')
                    misc.create_trackFile(shortcuts.star_whole_genome_indexing_complete)
                    print("\nWhole genome indexing completed!\n")
                    time.sleep(5)

            # Index parts of genome
            elif choice == 2:
                    if misc.step_allready_completed(f'{ref_dir}{filename}/star_index/starIndex.complete'):
                        misc.log_to_file(f'{filename} genome indexing with star allready completed...')
                    else:
                        threads = multiprocessing.cpu_count() - 2
                        cmd_StarIndex = f'''
                        STAR --runThreadN {threads} \\
                        --genomeSAindexNbases 12 \\
                        --runMode genomeGenerate \\
                        --genomeDir {ref_dir}{filename}/star_index \\
                        --genomeFastaFiles {ref_dir}{filename}/{filename}.fa \\
                        --sjdbGTFfile {ref_dir}{filename}/{filename}.gtf'''
                        misc.run_command(cmd_StarIndex, '\nIndexing parts of genome completed')
                        misc.create_trackFile(f'{ref_dir}{filename}/star_index/starIndex.complete')
        except Exception as e:
            print(f'Error with ReferenceGenome.index_genome_rna in reference_genome.py: {e}')
            input('press any key to exit')
            sys.exit()
