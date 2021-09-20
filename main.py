# coding=utf-8
# Packages used in script
from os import sys
import argparse
from menus import Menus
from miscellaneous import Misc
from shortcuts import Shortcuts
import multiprocessing as mp
from setup_anaconda3 import SetupAnaconda3
from rna_seq_analysis import RnaSeqAnalysis
from dna_seq_analysis import DnaSeqAnalysis
from reference_genome import ReferenceGenome
import time
import timeit
import signal









########################################################### Main program starts #######################################################################
# metavar makes the help text more tidy


def main():

    # Accept arguments for the program in the command line
    parser = argparse.ArgumentParser(description='''This script is used to simplify and make ASE-analysis more time efficient''')
    parser.add_argument("-t", "--tumor_id", metavar="", required=True, help="Input clinical id of tumor samples")
    parser.add_argument("-n", "--normal_id", metavar="", required=True, help="Input clinical id of normal samples")
    parser.add_argument("-sg", "--subgroup", metavar="", required=True, help="Input subgroup of your sample (STR)")
    parser.add_argument("-T", "--threads", metavar="", required=True, help="Input number of CPU threads to use (INT)")
    parser.add_argument("-v", "--validate", metavar="", required=False, help="Choose to validate BAM (True/False), default=False")

    # Assign a variable to all written arguments so that they can be passed to different classes etc
    options = parser.parse_args()
    
    # Assign variables to the imported classes
    misc = Misc()
    all_menus = Menus(misc) # Forward misc to Menus so that all methods in misc can be accessed by the Menus class
    shortcuts = Shortcuts(options) # Forward options to Shortcuts so that all arguments in options can be accessed by the Shortcuts class
    rna_analysis = RnaSeqAnalysis()
    dna_analysis = DnaSeqAnalysis()
    ref_genome = ReferenceGenome()
    setup = SetupAnaconda3()

    # Handle what happens if the user presses ctrl + c
    def signal_handler(sig, frame):
        misc.log_to_file("INFO", "Aborted by user!")
        sys.exit()

    # listens for ctrl + c
    signal.signal(signal.SIGINT, signal_handler)

    misc.log_to_file("INFO", f"--tumor_id: {options.tumor_id} --normal_id: {options.normal_id} --subgroup: {options.subgroup} --thread: {options.threads} --validate: {options.validate}")
   

     # Handle what happens if the user assigns more CPU threads than available
    if int(options.threads) > mp.cpu_count():
        misc.log_to_file("error", f"Threads to use: {options.threads} > Available threads: {mp.cpu_count()}")
        print(f"Threads to use: {options.threads} > Available threads: {mp.cpu_count()}")
        sys.exit()

     # Handle what happens if the user assigns exactly the number of CPU threads that is avaliable
    elif int(options.threads) == mp.cpu_count():
        misc.log_to_file("warning", f"Threads to use: {options.threads} = Available threads: {mp.cpu_count()}\nIt is not recommended to use all available threads\nComputer may be unstable")
        input("Press any key to continue")



    misc.log_to_file("info", "-----Program starts-----\n")
    while True:
        # Main menu
        menu_choice = all_menus.menu(misc, all_menus.main_menu)
        misc.log_to_file("info", "Program at main menu")
        if menu_choice == "":
            misc.log_to_file("info", "User input: exit program\n")
            break

        # Setyp anaconda3 environment
        if menu_choice == '1':
            misc.log_to_file("info","User input: 1, Setup anaconda3 environment")
            misc.clear_screen()
            setup.create_anaconda_environment(misc, shortcuts)

        # Dna analysis menu
        elif menu_choice == '2':
            misc.log_to_file("info", "User input: 2. DNA-analysis")
            while True:
                misc.log_to_file("info", "Program at DNA-analysis menu")
                dna_menu_choice = all_menus.menu(misc, all_menus.dna_menu)
                if dna_menu_choice == '':
                    misc.log_to_file("info", "User input: return to main menu")
                    break

                # Setup reference genome
                elif dna_menu_choice == '1':
                    misc.log_to_file("info", "User input: 1. Setup reference genome")
                    while True:
                        misc.log_to_file("info", "Program at Setup reference genome menu")
                        reference_genome_menu_choice = all_menus.menu(misc, all_menus.reference_genome_menu)

                        if reference_genome_menu_choice == '':
                            misc.log_to_file("info", "User input: return to DNA-analysis menu")
                            break
                        # Download reference genome
                        elif reference_genome_menu_choice == '1':
                            misc.log_to_file("info", "User input: 1. Download reference genome\n")
                            ref_genome.download(misc, shortcuts)
                        # Index reference genome
                        elif reference_genome_menu_choice == '2':
                            misc.log_to_file("info", "User input: 2. Index reference genome\n")
                            dna_analysis.index_genome_dna(misc, shortcuts)
                            break

                # Create library list file
                elif dna_menu_choice == '2':
                    misc.log_to_file("info", "User input: 2. Create library list file\n")
                    misc.clear_screen()
                    misc.validate_id(options, shortcuts)
                    all_menus.build_library_dna_menu(options, misc, shortcuts)

                # Run dna analysis
                elif dna_menu_choice == '3':
                    start = timeit.default_timer()
                    misc.log_to_file("info", "User input: 3. Run analysis\n")
                    misc.clear_screen()
                    misc.validate_id(options, shortcuts)
                    dna_analysis.alignment(options, misc, shortcuts)
        
                    if options.validate == True: dna_analysis.validate_bam_dna(options, misc, shortcuts)
                    dna_analysis.sort(options, misc, shortcuts)
                    dna_analysis.merge(options, misc, shortcuts)
                    dna_analysis.remove_duplicate(options, misc, shortcuts)
                    dna_analysis.realign(options, misc, shortcuts)
                    dna_analysis.gatk_haplotype(options, misc, shortcuts)
                    dna_analysis.delly(options, misc, shortcuts)
                    dna_analysis.manta(options, misc, shortcuts)
                    elapsed = timeit.default_timer() - start
                    misc.log_to_file("info", f'GDC DNA-Seq analysis pipeline successfully completed in {misc.elapsed_time(elapsed)} - OK!')
                    sys.exit()

        # Rna analysis menu
        elif menu_choice == '3':
            misc.log_to_file("info", "User input: RNA-analysis")
            while True:
                misc.log_to_file("info", "Program at RNA-analysis menu")
                misc.clear_screen()
                rna_choice = all_menus.menu(misc, all_menus.rna_menu)
                if rna_choice == '':
                    misc.log_to_file("info", "User input: return to previous menu")
                    break

                # Index reference genome
                elif rna_choice == '1':
                    misc.log_to_file("info", "User input: index reference genome\n")
                    rna_analysis.index_genome_rna(options, misc, shortcuts)

                # Map reads to reference genome
                elif rna_choice == '2':
                    misc.log_to_file("info", "User input: Map reads to reference genome\n")
                    rna_analysis.map_reads(options, misc, shortcuts)
                    rna_analysis.ASEReadCounter(options, misc, shortcuts)
                    rna_analysis.add_wgs_data_to_csv(options, misc, shortcuts)
                    sys.exit()





if __name__ == '__main__':
    main()
