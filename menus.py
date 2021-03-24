from os import path, getenv, listdir, makedirs, sys, remove, kill, getppid
import signal
import subprocess
import time
import logging
logging.basicConfig(filename=getenv("HOME")+'/sequencing_project/Logfile.txt', filemode='a', format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

try:
    from Bio import SeqIO
except:
    pass


class Menus():

    def __init__(self):
        self.main_menu = (['Setup Anaconda3 environment', 'DNA-analysis', 'RNA-analysis'], "\033[1mMain menu\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to exit program)")
        self.reference_genome_menu = (['Download reference genome', 'Index reference genome'], "\033[1mSetup reference genome menu\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to return to previous menu)")
        self.reference_genome_index_menu = (['Index whole genome', 'Index parts of genome'], "\033[1mIndex reference genome menu\033[0m\n" + "-"*28, "(leave blank to return to main menu)")
        self.dna_menu = (['Setup reference genome', 'Create library list file', 'Run analysis'], "\033[1mDNA-analysis menu\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to return to main menu)")
        self.rna_menu = (['Index reference genome', 'Map reads to reference genome'], "\033[1m""RNA-analysis menu""\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to return to main menu)")
        self.map_reads_menu = (['Map reads to whole genome', 'Map reads to parts of genome'], "\033[1mMap reads to reference genome menu\033[0m\n" + "-"*28, "(leave blank to return to main menu)")



    #---------------------------------------------------------------------------
    def info_script(self):
        '''Prints information about the creator of the program'''

        print( "\033[1m""Biomedswe Allele-specific expression analyser pipeline (BASEAP) v.1.0. 2021\n" "\033[0m")
        print("This program was created during a masterproject about B-cell precursor acute lymphoblastic leukemia (BCP-ALL) in 2020-2021.\n")
        print("Created by:\nJonas Andersson\nMaster's programme in biomedicine\nLund University\nLund, Sweden\n")
        print("Github: https://github.com/biomedswe/sequencing_project")
        print("For correspondence please contact jonas870318@gmail.com\n\n")

    #---------------------------------------------------------------------------
    def menu(self, misc, choices):
        misc.clear_screen()
        self.info_script()
        print(choices[1], "\n")
        for i, choice in enumerate(choices[0], start=1):
            print(f'{i}. {choice}')
        return misc.validate_choice(len(choices[0]), choices[2])

    #---------------------------------------------------------------------------
    def build_library_dna_menu(self, options, misc, shortcuts):
        '''This function lists all WGS files in the directory and writes them to a library-list text file.
           You can create either a file for single-end protocol or paired-end protocol.
           This file is then used in the dna analysis'''

        while True:
            print("Create library list file.\n\n")
            print("1. Single end sequencing\n")
            print("2. Paired end sequencing\n")
            choice = misc.validate_choice(2, "(leave blank to return to DNA-analysis menu)")
            if choice == "":
                misc.log_to_file('User input: return to DNA-analysis menu')
                return ""
            elif misc.confirm_choice():
                misc.log_to_file(f'User input: confirmed choice: {choice}')
                misc.clear_screen()
                files = listdir(shortcuts.dna_reads_dir)
                with open(f"{shortcuts.dna_seq_dir}library.txt", 'w') as out_file:
                    for line, library_id in enumerate(files, start=1):
                        if choice == '1': # Single-end sequencing
                            if options.tumor_id in library_id:
                                out_file.write(f"{options.tumor_id} {library_id[:-9]} {library_id} N/A\n")
                            else:
                                out_file.write(f"{options.normal_id} {library_id[:-9]} {library_id} N/A\n")

                        elif choice == '2': # Paired-end sequencing
                            if (line % 2) == 1: # not even
                                if options.tumor_id in library_id:
                                    out_file.write(f"{options.tumor_id} {library_id[:-11]} {library_id} ")
                                else:
                                    out_file.write(f"{options.normal_id} {library_id[:-11]} {library_id} ")
                            else: # even
                                out_file.write(f"{library_id}\n")
                misc.log_to_file('Library list file created!')
                print("\nNow that you have created your list library file, you can run the analysis!\n")
                input("Press any key to return to DNA-analysis menu...")
                return
            else:
                misc.clear_screen()
                continue

    #---------------------------------------------------------------------------
    # this can be use for both dna and rna genome index
    def star_index_menu(self, choice, misc, shortcuts, rna_analysis):
        '''This function creates a STAR index with STAR --runMode genomeGenerate'''

        while True:

            if choice == "":
                return ""
            elif misc.confirm_choice():
                misc.clear_screen()
                # whole genome
                if choice == '1':
                    rna_analysis.index_genome_rna_analysis(1, misc, shortcuts)
                # parts of genome
                elif choice == '2':
                    name = rna_analysis.index_genome_rna_analysis(2, misc, shortcuts)
                    return name
            else:
                continue

#-------------------------------------------------------------------------------
class Misc():
    '''This class contains miscellaneous functions related to general functionality'''

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def log_to_file(self, text):
        print(f'\n{text}')
        logging.info(text)

    #---------------------------------------------------------------------------
    def validate_id(self,options, shortcuts):
        '''This function validates wether the tumor_id and normal_id you entered is present in your reads'''

        try:
            if options.tumor_id and options.normal_id in "".join(listdir(shortcuts.dna_reads_dir)):
                self.log_to_file(f"tumor_id {options.tumor_id} and normal_id {options.normal_id} correctly validated")

            else:
                self.log_to_file(f'You have entered a tumor_id: {options.tumor_id} and/or normal_id: {options.normal_id} that is not present in your reads.\nPlease restart program and verify that you have typed in the right \"clinical_id\" for both tumor (-t) and normal (-n)!')
                input("Press any key to exit program")
                sys.exit()
        except Exception as e:
            self.log_to_file(f'Error with validate_id() in menus.py: {e}')



    #---------------------------------------------------------------------------
    def validate_choice(self, choices, text):
        '''This function checks if the input choice is valid'''

        while True:
            choice = input(text)
            if choice == "":
                return ""
            elif choice not in [str(i) for i in range(1,choices+1)]:
                print("Invalid choice, try again!")
            else:
                return choice

    #---------------------------------------------------------------------------
    def confirm_choice(self):
        '''This function asks if you want to confirm your choice'''

        while True:
            choice = input("Are you sure? y or n: ")
            if choice.lower() == "y":
                return True
            elif choice.lower() == 'n':
                return False
            else:
                print("Invalid choice, type 'n' or 'y'!")

    #---------------------------------------------------------------------------
    def clear_screen(self):
        '''This function clears the terminal window'''
        subprocess.run("clear", shell=True)

    #---------------------------------------------------------------------------
    def close_terminal(self):
        '''This function closes the terminal'''
        kill(getppid(), signal.SIGHUP)

    #---------------------------------------------------------------------------
    def create_directory(self, paths):
        '''This function creates the output directories for the different analysis steps'''
        for path in paths:
            try:
                # create target directory
                makedirs(path)
            except FileExistsError:
                pass

    #---------------------------------------------------------------------------
    def step_allready_completed(self, file):
        '''This function checks if an analysis step is completed by looking after a created file *.complete'''
        if path.isfile(file):
            return True
        else:
            return False

    #---------------------------------------------------------------------------
    def create_trackFile(self, file):
        '''This function creates utput list but also a file after each step as a marker that the step is completed'''
        try:
            with open(file, 'w'):
                pass
        except Exception as e:
            print(f'Error with index_genome_dna: {e}')


    #---------------------------------------------------------------------------
    def run_command(self, command):
        '''This function executes a command and checks if it was executes without errors'''
        try:
            return_code = subprocess.run(command, shell=True)
            if return_code.returncode == 0:
                return 0
            else:
                self.log_to_file('Process ended with returncode != 0, see logfile.txt for more information. Exiting program...')
                self.log_to_file(return_code.stderr)
                self.log_to_file(return_code.stdout)
                sys.exit()
        except Exception as e:
            self.log_to_file(f'Error with misc.run_command() in menus.py: {e}')
            input("Press any key to continue")



    #---------------------------------------------------------------------------
    def choose_chromosomes_to_index(self, menus, shortcuts):
        '''Takes one or more chromosome as input, check if syntax is valid and if so, returns chromosomes as a list'''

        self.clear_screen()
        menus.info_script()

        print('''You can add chromosomes separated by a space.
Use this syntax:

chr1, chr2, chr3, ..., chr22
chrX, chrY, chrM
e.g. "chr11 chr12"
''')

        list = ['chr'+str(i) for i in range(1,23)]
        list.extend(['chrX', 'chrY', 'chrM'])

        while True:
            chromosomes = [chr for chr in input("Enter chromosomes to index (leave blank to return to main menu): ").split()] # Splits all input into a list with different entries
            if not chromosomes:
                return
                # print('You must enter at least on chromosome')
            elif all(chr in list for chr in chromosomes): # Checks if all entered chromosomes are in valid syntax by comparing to entries in list
                return chromosomes
            else:
                print('Invalid syntax, please check spelling!')

    #---------------------------------------------------------------------------
    def create_new_fasta(self, chromosomes, shortcuts):
        ref_file = shortcuts.reference_genome_file
        ref_dir = shortcuts.reference_genome_dir

        filename = "".join(chromosomes) + "_GRCh38"
        if self.step_completed(f'{ref_dir}{filename}/{filename}.fa', f'Fasta for {filename} allready created, skips step...'):
            pass
        else:
            print('Creating a new fasta file...')
            self.create_directory([f'{ref_dir}{filename}'])
            sequences = SeqIO.parse(ref_file, 'fasta')
            with open(f'{ref_dir}{filename}/{filename}.fa', 'w+') as fa:
                for chr in chromosomes:
                    for line in sequences:
                        if line.id == chr:
                            fa.write(">" + str(line.id) + "\n")
                            fa.write(str(line.seq)+ "\n\n")
                            break # I use break here, otherwise it will continue with chr10,12,13 etc. if i choose chr1

        return filename

    #---------------------------------------------------------------------------
    def create_new_gtf(self, chromosomes, filename, shortcuts):
        '''Create a new gtf file from choosed chromosomes'''

        ref_dir = shortcuts.reference_genome_dir

        if self.step_completed(f'{ref_dir}{filename}/{filename}.gtf', f'Gtf for {filename} allready created, skips step...'):
            pass
        else:
            print('Creating a new gtf file...')
            sequences = SeqIO.parse(shortcuts.reference_genome_file, 'fasta')
            with open(f'{shortcuts.reference_genome_dir}{filename}/{filename}.bed', 'w') as bed:
                for chr in chromosomes:
                    for line in sequences:
                        if line.id == chr:
                            bed.write(str(line.id) + "\t")
                            bed.write("0\t")
                            bed.write(str(len(line.seq)))
                            break
            cmd_createGTF = f"bedtools intersect -a {shortcuts.annotation_gtf_file} -b {shortcuts.reference_genome_dir}{filename}/{filename}.bed > {shortcuts.reference_genome_dir}{filename}/{filename}.gtf"
            misc.run_command(cmd_createGTF, '')
            remove(f'{shortcuts.reference_genome_dir}{filename}/{filename}.bed')


class Shortcuts():
    '''This class contains shortcuts to key files and folders neccessary for the program'''

    def __init__(self, options):

        '''skapa funktioner av dessa'''

        # Shortcut to main folders
        self.sequencing_project_dir = getenv("HOME")+"/sequencing_project/"
        self.dna_seq_dir = f"{self.sequencing_project_dir}dna_seq/"
        self.rna_seq_dir =  f"{self.sequencing_project_dir}rna_seq/"

        # Shortcuts to input folders
        self.dna_reads_dir  = f"{self.dna_seq_dir}reads/"
        self.reference_genome_dir = f"{self.sequencing_project_dir}reference_genome/"
        self.GRCh38_dir = f"{self.reference_genome_dir}GRCh38/"

        # Shortcuts to output folders in DNA sequencing analysis
        self.aligned_output_dir = f"{self.dna_seq_dir}aligned/"
        self.sorted_output_dir = f"{self.dna_seq_dir}sorted/"
        self.merged_output_dir = f"{self.dna_seq_dir}merged/"
        self.removed_duplicates_output_dir = f"{self.dna_seq_dir}removed_duplicates/"
        self.realigned_output_dir = f"{self.dna_seq_dir}realigned/"
        self.haplotypecaller_output_dir = f"{self.dna_seq_dir}gatk_haplotypecaller/"
        self.delly_output_dir = f"{self.dna_seq_dir}delly/"
        self.manta_output_dir = f"{self.dna_seq_dir}manta/"
        self.manta_variants_dir = f"{self.dna_seq_dir}manta/results/variants/"

        # Shortcuts to files used in DNA sequencing analysis
        self.reference_genome_file = f"{self.GRCh38_dir}GRCh38.p13.genome.fa"
        self.reference_genome_exclude_template_file = f"{self.sequencing_project_dir}excludeTemplate/human.hg38.excl.tsv"
        self.configManta_file = getenv("HOME")+"/anaconda3/envs/sequencing/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py"
        self.runWorkflow_file = getenv("HOME")+"/sequencing_project/dna_seq/Manta/runWorkflow.py"

        # Shortcuts to folders used in RNA sequencing analysis
        self.rna_reads_dir  = f"{self.rna_seq_dir}reads/"
        self.star_output_dir = f"{self.rna_seq_dir}star/"
        self.star_index_dir_whole_genome =  f"{self.GRCh38_dir}star_index/"


        # Shortcuts to files used in RNA sequencing analysis
        self.annotation_gtf_file = f"{self.GRCh38_dir}gencode.v37.primary_assembly.annotation.gtf"
        self.gatk_vcfFile = f"{self.haplotypecaller_output_dir}{options.tumor_id}_filtered_ReadDepthOver10_het.vcf"

        # Shortcuts to output lists (used for input in pipeline steps) (and also for validation if pipeline step i allready completed)
        self.alignedFiles_list = f"{self.aligned_output_dir}alignedFiles.txt"
        self.sortedFiles_list = f"{self.sorted_output_dir}sortedFiles.txt"
        self.mergedFiles_list = f"{self.merged_output_dir}mergedFiles.txt"
        self.removeDuplicates_list = f"{self.removed_duplicates_output_dir}remove_duplicate.txt"
        self.realignedFiles_list = f"{self.realigned_output_dir}realignedFiles.txt"

        # Shortcuts to files used to validate if pipeline step is allready completed
        self.anaconda_setup_complete = getenv("HOME")+'/anaconda3/install.complete'
        self.bwa_index_whole_reference_genome_complete = f"{self.GRCh38_dir}bwa.complete"
        self.validate_bam_complete = f"{self.aligned_output_dir}validateBam.complete"
        self.haplotypecaller_complete = f"{self.haplotypecaller_output_dir}haplotypeCaller.complete"
        self.delly_complete = f"{self.delly_output_dir}delly.complete"
        self.manta_complete = f"{self.manta_output_dir}manta.complete"
        self.star_whole_genome_indexing_complete = f"{self.star_index_dir_whole_genome}whole_index.complete"
