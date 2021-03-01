import subprocess

class Menus():

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def info_script(self):
        '''Prints information about the creator of the program'''

        print("This program was created by biomedicine master student Jonas Andersson during a masterproject 2020-2021.")
        print("For correspondence please contact jonas870318@gmail.com\n")
        print("Biomedswe 2021 v.1.0")
        print("Automated DNA/RNA-sequencing script\n")

    #---------------------------------------------------------------------------
    def validate_choice(choices,text):
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
    def main_menu(self):
        '''Prints an introduction to the program where you can choose what to do'''

        self.info_script()
        print("1. Anaconda3 setup\n")
        print("2. Download reference genome\n")
        print("3. DNA-analysis\n")
        print("4. RNA-analysis\n")
        choices = 4
        text = "(leave blank to exit)"
        choice = self.validate_choice(choices,text)
        return choice

    #---------------------------------------------------------------------------
    def anaconda_menu(self):
        '''Prints an introduction to the program where you can choose what to do'''

        self.info_script()
        print("This script lets you install Anaconda3 and set up an environment with the required software packages.\n")
        print("If you already have installed anaconda, skip first step and proceed with setting up an environment\n")
        print("1. Download and install Anaconda3 with python 3.8")
        print("2. Set up a new conda environment for DNA and RNA-sequence analysis\n")
        choices = 2
        text = "(leave blank to return to main menu)"
        choice = validate_choice(choices, text)
        return choice

    #---------------------------------------------------------------------------
    def dna_menu(self):
        '''Prints an introduction to the program where you can choose either to create a library list file or run analysis'''

        info_script()
        print("This script makes an automated analysis of WGS/WES reads.\nBefore running analysis, you must have indexed the reference genome and also have a list file containing all your reads.\nTo complete these requirements, please run option 1 and 2 first\n")
        print("1. Index reference genome")
        print("2. Create library list file")
        print("3. Run analysis")
        choice = validate_choice(3, "(leave blank to return to main menu)")
        return choice

    #---------------------------------------------------------------------------
    def build_library_menu(self, options):
        '''This function lists all WGS files in the directory and writes them to a library-list text file.
           You can create either a file for single-end protocol or paired-end protocol.
           This file is then used in the dna analysis'''

        while True:
            print("Create library list file.\n\n")
            print("1. Single end sequencing\n")
            print("2. Paired end sequencing\n")
            choice = validate_choice(2, "(leave blank to return to DNA-analysis menu)")
            if choice == "":
                return ""
            elif self.confirm_choice():
                self.clear_screen()
                files = listdir(reads_dir)
                with open("dna_seq/library.txt", 'w') as out_file:
                    for line, library_id in enumerate(files, start=1):
                        if choice == '1': # Single-end sequencing
                            if options.tumor_id in library_id:
                                out_file.write(f"{options.tumor_id} {library_id[:-6]} {library_id} N/A\n")
                            else:
                                out_file.write(f"{options.normal_id} {library_id[:-6]} {library_id} N/A\n")
                        elif choice == '2': # Paired-end sequencing
                            if (line % 2) == 1: # not even
                                if options.tumor_id in library_id:
                                    out_file.write(f"{options.tumor_id} {library_id[:-8]} {library_id} ")
                                else:
                                    out_file.write(f"{options.normal_id} {library_id[:-8]} {library_id} ")
                            else: # even
                                out_file.write(f"{library_id}\n")
                print("Library list file created!\nNow that you have created your list library file, you can run the analysis!\n")
                input("Press any key to return to DNA-analysis menu...")
                return ''
            else:
                clear_screen()
                continue

    #---------------------------------------------------------------------------
    def rna_menu(self):
        '''Prints an introduction to the program where you can choose either to create a STAR index or run analysis'''

        self.info_script()
        print("This script makes an automated analysis of RNA.\nBefore mapping reads to the genome, you must have created a genome index.\n")
        print("1. Create a genome index")
        print("2. Map reads to the genome")
        choice = self.validate_choice(2, "(leave blank to return to main menu)")
        return choice

    #---------------------------------------------------------------------------
    def star_index_menu(self):
        '''This function creates a STAR index with STAR --runMode genomeGenerate'''

        print('Create genome index\n\n')
        print('1. Index whole genome')
        print('2. Index parts of genome')
        choice = self.validate_choice(2, "(leave blank to return to previous menu)")
        if choice == "":
            return ""
        elif self.confirm_choice():
            self.clear_screen()
            if choice == '1':
                index_genome(1)
            elif choice == '2':
                name = index_genome(2)
        return
