from os import listdir

class Menus():


    def __init__(self, misc):
        try:
            self.main_menu = (['Setup Anaconda3 environment', 'DNA-analysis', 'RNA-analysis'], "\033[1mMain menu\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to exit program)")
            self.reference_genome_menu = (['Download reference genome', 'Index reference genome'], "\033[1mSetup reference genome menu\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to return to previous menu)")
            self.dna_menu = (['Setup reference genome', 'Create library list file', 'Run analysis'], "\033[1mDNA-analysis menu\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to return to main menu)")
            self.rna_menu = (['Index reference genome', 'Map reads to reference genome'], "\033[1m""RNA-analysis menu""\033[0m\n" + "-"*31 + "\nRun the options below in order:", "(leave blank to return to main menu)")
        except Exception as e:
            misc.log_exception('.__init__() in menus.py:', e)



    #---------------------------------------------------------------------------
    def info_script(self, misc):
        '''Prints information about the creator of the program'''
        try:
            print( "\033[1m""Biomedswe Allele-specific expression analyser pipeline (BASE) v.1.0. 2021\n" "\033[0m")
            print("This program was created during a masterproject about B-cell precursor acute lymphoblastic leukemia (BCP-ALL) in 2020-2021.\n")
            print("Created by:\nJonas Andersson\nMaster's programme in biomedicine\nLund University\nLund, Sweden\n")
            print("Github: https://github.com/biomedswe/sequencing_project")
            print("For correspondence please contact jonas870318@gmail.com\n\n")
        except Exception as e:
            misc.log_exception('.info_script() in menus.py:', e)

    #---------------------------------------------------------------------------
    def menu(self, misc, options):
        '''This functions creates a general menu with n choices'''
        try:
            misc.clear_screen()
            self.info_script(misc)
            print(options[1], "\n")
            for i, option in enumerate(options[0], start=1):
                print(f'{i}. {option}')
            return misc.validate_choice(len(options[0]), options[2])
        except Exception as e:
            misc.log_exception('.menu() in menus.py:', e)

    #---------------------------------------------------------------------------
    def build_library_dna_menu(self, options, misc, shortcuts):
        '''This function lists all WGS files in the directory and writes them to a library-list text file.
           You can create either a file for single-end protocol or paired-end protocol.
           This file is then used in the dna analysis'''

        try:
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
                    files = sorted(listdir(shortcuts.dna_reads_dir))
                    with open(f"{shortcuts.dna_seq_dir}{options.tumor_id}_library.txt", 'w') as out_file:
                        for line, library_id in enumerate(files, start=1):
                            if choice == '1': # Single-end sequencing
                                if options.tumor_id in library_id:
                                    out_file.write(f"{options.tumor_id} {library_id.split('_')[0]} {library_id} N/A\n")
                                else:
                                    out_file.write(f"{options.normal_id} {library_id.split('_')[0]} {library_id} N/A\n")

                            elif choice == '2': # Paired-end sequencing
                                if (line % 2) == 1: # not even
                                    if options.tumor_id in library_id:
                                        out_file.write(f"{options.tumor_id} {library_id.split('_')[0]} {library_id} ")
                                    else:
                                        out_file.write(f"{options.normal_id} {library_id.split('_')[0]} {library_id} ")
                                else: # even
                                    out_file.write(f"{library_id}\n")
                    misc.log_to_file('Library list file created!')
                    print("Now that you have created your list library file, you can run the analysis!\n")
                    input("Press any key to return to DNA-analysis menu...")
                    return
                else:
                    misc.clear_screen()
                    continue
        except TypeError as t:
            misc.log_exception("Missing tumor id and/or normal id!\nplease enter: \"-t <tumor-id> -n <normal-id>\" when running the script", t)
        except Exception as e:
            misc.log_exception('.build_library_dna_menu() in menus.py:', e)

#-------------------------------------------------------------------------------
