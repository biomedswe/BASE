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

<<<<<<< HEAD
    
=======
    #---------------------------------------------------------------------------
    def index_genome_rna(self, choice, filename, misc, shortcuts):
        '''This function indexes either the whole genome or the chromosomes entered'''

        try:

            ref_dir = shortcuts.reference_genome_dir

            # Index whole genome
            if choice == 1:
                if misc.step_allready_completed(shortcuts.star_whole_genome_indexing_complete):
                    misc.logfile('Whole genome indexing allready completed, returning...')
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
                        misc.logfile(f'{filename} genome indexing with star allready completed...')
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
            print(f'Error with index_genome_rna: {e}')
            input("Press any key to continue...")
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
