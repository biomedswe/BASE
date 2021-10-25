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
            print( "\033[1m""Biomedswe Allele-Specific Expression analyser (BASE) v.1.0. 2021\n" "\033[0m")
            print("This program was created during a master thesis project\nabout B-cell precursor acute lymphoblastic leukemia (BCP-ALL) in 2020-2021.\n")
            print("Developer:\nJonas Andersson\nMaster's programme in biomedicine\nDivision of Clinical Genetics\nLund University, BMC C13\nSE-221 84 Lund, Sweden\n")
            print("Acknowledgment:\nSpecial thanks to Prof. Kajsa Paulsson and assistant researcher Minjun Yang\n")
            print("Github: https://github.com/biomedswe/BASE")
            print("For correspondence please contact jonas.andersson@med.lu.se\n\n")
        except Exception as e:
            misc.log_exception('.info_script() in menus.py:', e)

    #---------------------------------------------------------------------------
    def menu(self, misc, meny_options):
        '''This functions creates a general menu with n choices'''
        try:
            misc.clear_screen()
            self.info_script(misc)
            print(meny_options[1], "\n")
            for i, option in enumerate(meny_options[0], start=1):
                print(f'{i}. {option}')
            return misc.validate_choice(len(meny_options[0]), meny_options[2])
        except Exception as e:
            misc.log_exception('.menu() in menus.py:', e)

    #---------------------------------------------------------------------------
    def build_library_dna_menu(self, options, misc, shortcuts):
        '''This function lists all WGS files in the directory and writes them to a library-list text file.
           You can create either a file for single-end protocol or paired-end protocol.
           This file is then used in the dna analysis'''

        try:
            misc.clear_screen()
            files = sorted(listdir(shortcuts.dna_reads_dir))
            with open(f"{shortcuts.dna_seq_dir}{options.tumor_id}_library.txt", 'w') as out_file:
                for line, library_id in enumerate(files, start=1):
                    
                    if (line % 2) == 1: # not even
                        
                        if options.tumor_id in library_id:
                            out_file.write(f"{options.tumor_id} {library_id.split('_')[0]} {library_id} ")
                        
                        else:
                            out_file.write(f"{options.normal_id} {library_id.split('_')[0]} {library_id} ")
                    
                    else: # even
                        out_file.write(f"{library_id}\n")
                        
                misc.log_to_file("info", "Library list file created!")
                print("Now that you have created your list library file, you can run the analysis!\n")
                input("Press any key to return to DNA-analysis menu...")
                return

        except TypeError as t:
            misc.log_exception("Missing tumor id and/or normal id!\nplease enter: \"-t <tumor-id> -n <normal-id>\" when running the script", t)
        except Exception as e:
            misc.log_exception('.build_library_dna_menu() in menus.py:', e)

#-------------------------------------------------------------------------------
