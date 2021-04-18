from os import path, getenv, listdir, makedirs, sys, remove, kill, getppid
import signal
import subprocess
import logging
import timeit
logging.basicConfig(filename=getenv("HOME")+'/sequencing_project/Logfile.txt', filemode='a', format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

try:
    from Bio import SeqIO
except:
    pass

class Misc():
    '''This class contains miscellaneous functions related to general functionality'''

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def choose_chromosomes_to_index(self, menus, shortcuts):
        '''Takes one or more chromosome as input, check if syntax is valid and if so, returns chromosomes as a list'''

        try:
            self.clear_screen()
            menus.info_script(self)

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
        except Exception as e:
            self.log_exception(".choose_chromosomes_to_index() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def clear_screen(self):
        '''This function clears the terminal window'''
        try:
            subprocess.run("clear", shell=True)
        except Exception as e:
            self.log_exception(".clear_screen() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def close_terminal(self):
        '''This function closes the terminal'''
        try:
            kill(getppid(), signal.SIGHUP)
        except Exception as e:
            self.log_exception(".close_terminal() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def confirm_choice(self):
        '''This function asks if you want to confirm your choice'''

        try:
            while True:
                choice = input("Are you sure? y or n: ")
                if choice.lower() == "y":
                    return True
                elif choice.lower() == 'n':
                    return False
                else:
                    print("Invalid choice, type 'n' or 'y'!")
        except Exception as e:
            self.log_exception(".confirm_choice() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def create_directory(self, paths):
        '''This function creates the output directories for the different analysis steps'''

        for path in paths:
            try:
                # create target directory
                makedirs(path)
                self.log_to_file(f"{path} folder(s) created succesfully - OK!")
            except FileExistsError:
                self.log_to_file(f"{path} allready exists - skips step...")
            except Exception as e:
                self.log_exception(".create_directory() in misc.py:", e)
                sys.exit()

    #---------------------------------------------------------------------------
    def create_new_fasta(self, chromosomes, shortcuts):
        try:
            ref_dir = shortcuts.reference_genome_dir
            ref_file = shortcuts.reference_genome_file
            filename = "".join(chromosomes) + "_GRCh38.p13.genome"
            if not self.step_allready_completed(f'{ref_dir}{filename}/{filename}.fa', f'Creating Fasta for {filename}'):
                start = timeit.default_timer()
                self.log_to_file(f'Starting: creating a new fasta file for {filename}...')
                self.create_directory([f'{ref_dir}{filename}'])
                sequences = SeqIO.parse(ref_file, 'fasta')
                with open(f'{ref_dir}{filename}/{filename}.fa', 'w+') as fa:
                    for chr in chromosomes:
                        for line in sequences:
                            if line.id == chr:
                                fa.write(">" + str(line.id) + "\n")
                                fa.write(str(line.seq)+ "\n\n")
                                break # I use break here, otherwise it will continue with chr10,12,13 etc. if i choose chr1
                end = timeit.default_timer()
                self.log_to_file(f"Creating fasta for {filename} succesfully completed in {self.elapsed_time(end-start)} - OK!)")
            return filename
        except Exception as e:
            self.log_exception(".create_new_fasta() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def create_new_gtf(self, chromosomes, filename, shortcuts):
        '''Create a new gtf file from choosed chromosomes'''
        try:
            ref_dir = shortcuts.reference_genome_dir

            if not self.step_allready_completed(f'{ref_dir}{filename}/{filename}.gtf', f'Creating Gtf for {filename}'):
                start = timeit.default_timer()
                self.log_to_file(f'Starting: creating a new gtf file for {filename}...')
                sequences = SeqIO.parse(shortcuts.reference_genome_file, 'fasta')
                with open(f'{ref_dir}{filename}/{filename}.bed', 'w') as bed:
                    for chr in chromosomes:
                        for line in sequences:
                            if line.id == chr:
                                bed.write(str(line.id) + "\t")
                                bed.write("0\t")
                                bed.write(str(len(line.seq)))
                                break
                cmd_createGTF = f"bedtools intersect -a {shortcuts.annotation_gtf_file} -b {ref_dir}{filename}/{filename}.bed > {ref_dir}{filename}/{filename}.gtf"
                self.run_command(cmd_createGTF, None, None, None)
                end = timeit.default_timer()
                self.log_to_file(f"Creating Gtf for {filename} succesfully completed in {self.elapsed_time(end-start)} - OK!)")
                remove(f'{ref_dir}{filename}/{filename}.bed')
        except Exception as e:
            self.log_exception(".create_new_gtf() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def create_outputList_dna(self, output_path, write_to_file):
        '''This function creates a list file containing the output name of the files created in the pipeline step where the function is used, the list is then used in the next step'''

        try:
            with open(output_path, 'a') as c:
                c.write(f"{write_to_file}\n")

        except Exception as e:
            misc.log_exception(".create_outputList_dna() in misc.py:", e)
    #---------------------------------------------------------------------------
    def create_trackFile(self, file):
        '''This function creates a trackfile that step_allready_completed() function can look after when checking if step is allready completed'''
        try:
            if path.isfile(file):
                self.log_to_file("You are trying to owerwrite an existing file. please check the code.")
                sys.exit()
            with open(file, 'w'):
                pass
        except Exception as e:
            self.log_exception(".create_trackFile() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def elapsed_time(self, elapsed):
        try:
            if elapsed >= 86400:
                return f"{round(elapsed/86400, 2)} days"
            elif  3600 <= elapsed < 86400:
                return f"{round(elapsed/3600, 2)} hours"
            elif 60 <= elapsed < 3600:
                return f"{round(elapsed/60, 2)} minutes"
            elif elapsed < 60:
                return f"{round(elapsed, 2)} seconds"
            else:
                return 'error'
        except Exception as e:
            self.log_exception(".elapsed_time() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def log_exception(self, text, exception):
        try:
            self.log_to_file(f'{exception}: {text}. Exiting program...')
            sys.exit()
        except Exception as e:
            self.log_to_file(f"Error with .log_exception() in misc.py: {e}. Exiting program...")

    #---------------------------------------------------------------------------
    def log_to_file(self, text):
        try:
            print(f'{text}\n')
            logging.info(f"{text}\n")
        except Exception as e:
            logging.info(f'Error with {self}.log_to_file() in misc.py: {e}. Exiting program...')
            sys.exit()

    #---------------------------------------------------------------------------
    def remove_file(self, file):
        '''This function removes incomplete files if the processing of file ended with returncode != 0'''
        try:
            if file:
                remove(file)
                self.log_to_file(f'Incomplete {file} removed - OK!')
        except OSError: pass
        except Exception as e:
            self.log_exception("remove_file() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def run_command(self, command, text, file, trackfile):
        '''This function first calls step_allready_completed() to check if the step i allready completed.
        If not completed; if process is executed without errors, it prints to logfile with time taken and passes return.
        else; if process ends  with errors, it prints error to log, removes incomplete file and then exits program'''

        try:
            if file:
                if self.step_allready_completed(file, text):
                    return False
            return_code = subprocess.run(command, shell=True)
            if return_code.returncode == 0:
                if text: self.log_to_file(f"{text} succesfully completed - OK!")
                if trackfile: self.create_trackFile(trackfile)
                return True
            else:
                self.log_to_file(f'Process ended with returncode != 0: {command} - ERROR!')
                self.remove_file(file)
                sys.exit()
        except Exception as e:
            self.log_exception("run_command() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def step_allready_completed(self, file, text):
        '''This function checks if a step is allready completed by checking if "file" allready exists.
        If file exists: returns True, else return False'''

        try:
            if file:
                if path.isfile(file):
                    if text: self.log_to_file(f"{text} allready completed, skips step...")
                    return True
                else:
                    return False
            else:
                return False
        except Exception as e:
            self.log_exception(".step_allready_completed() in miscellaneous.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def validate_choice(self, choices, text):
        '''This function checks if the input choice is valid'''

        try:
            while True:
                choice = input(text)
                if choice == "":
                    return ""
                elif choice not in [str(i) for i in range(1,choices+1)]:
                    print("Invalid choice, try again!")
                else:
                    return choice
        except Exception as e:
            self.log_exception(".validate_choice() in misc.py:", e)
            sys.exit()

    #---------------------------------------------------------------------------
    def validate_id(self,options, shortcuts):
        '''This function validates wether the tumor_id and normal_id you entered is present in your reads'''

        try:
            if options.tumor_id and options.normal_id in "".join(listdir(shortcuts.dna_reads_dir)):
                self.log_to_file(f"tumor_id {options.tumor_id} and normal_id {options.normal_id} correctly validated")

            else:
                self.log_to_file(f'You have entered a tumor_id: {options.tumor_id} and normal_id: {options.normal_id} that is not present in your reads.\nPlease restart program and verify that you have typed in the right \"clinical_id\" for both tumor (-t) and normal (-n)!')
                input("Press any key to exit program")
                sys.exit()
        except Exception as e:
            self.log_exception(".validate_id() in misc.py:", e)
            sys.exit()
