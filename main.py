# Packages used in script
from os import getenv, sys, path, listdir, makedirs
# import subprocess
import argparse
from menus import Menus, Misc, Shortcuts
from rna_seq_analysis import RnaSeqAnalysis
from dna_seq_analysis import DnaSeqAnalysis
from reference_genome import ReferenceGenome
from setup_anaconda3 import SetupAnaconda3
import time
<<<<<<< HEAD
import timeit
=======
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510

########################################################### Main program starts #######################################################################
# metavar makes the help text more tidy


def main():
    parser = argparse.ArgumentParser(description='''This script is used to simplify installation and set up of Anaconda, download reference genome and perform DNA and RNA analysis. Please read the README file first!''')
    parser.add_argument("-t", "--tumor_id", metavar="", required=True, help="Input clinical id of tumor samples")
    parser.add_argument("-n", "--normal_id", metavar="", required=True, help="Input clinical id of normal samples")
    parser.add_argument("-i", "--intervals", metavar="", required=False, help="Input path to reference-genome interval if you have any (for use in SNV calling)")
    options = parser.parse_args() # all arguments will be passed to the functions
    filename = None
    # hur göra här? options måste med i shortcuts
    all_menus = Menus()
    misc = Misc()
    shortcuts = Shortcuts(options)
    rna_analysis = RnaSeqAnalysis()
    dna_analysis = DnaSeqAnalysis()
    ref_genome = ReferenceGenome()
    setup = SetupAnaconda3()
<<<<<<< HEAD

=======
    
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
    # main menu is shown untill choice is left blank
    misc.logfile('Program starts')
    while True:
        # Clears screen and prints main menu
        menu_choice = all_menus.menu(misc, all_menus.main_menu)
        misc.logfile('Program at main menu')
        if menu_choice == "":
            misc.logfile("User input: exit program")
            break

        # Setyp anaconda3 menu
        if menu_choice == '1':
            misc.logfile('User input (main menu): 1. Setup anaconda3')
            anaconda_choice = all_menus.menu(misc, all_menus.anaconda_menu)

<<<<<<< HEAD
            if anaconda_choice == '':
                misc.logfile('User input (Setup anaconda3 menu): (leave blank to return to main menu)')
                break

            elif anaconda_choice == '1':
                misc.logfile('User input (Setup anaconda3 menu): 1. Download and install Anaconda3 with python 3.7.6')
=======
            if anaconda_choice == '1':
                misc.logfile('Starting installation of Anaconda3')
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
                setup.install_anaconda(misc, shortcuts)

            elif anaconda_choice == '2':
                misc.logfile('User input (Setup anaconda3 menu): 2, Set up a new conda environment for DNA and RNA-sequence analysis')
                setup.create_anaconda_environment(misc, shortcuts)

<<<<<<< HEAD
        # DNA-analysis menu
        elif menu_choice == '2':
            misc.logfile('User input (main menu): 2. DNA-analysis')
            while True:
                misc.logfile('Program at DNA-analysis menu')
                dna_menu_choice = all_menus.menu(misc, all_menus.dna_menu)
                if dna_menu_choice == '':
                    misc.logfile('User input: return to previous menu')
=======




        # DNA-analysis menu
        elif menu_choice == '2':
            while True:
                dna_menu_choice = all_menus.menu(misc, all_menus.dna_menu)
                if dna_menu_choice == '':
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
                    break

                # Setup reference genome
                elif dna_menu_choice == '1':
<<<<<<< HEAD
                    misc.logfile('User input (DNA-analysis menu): 1. Setup reference genome ')
                    while True:
                        misc.logfile('Program at Setup reference genome menu')
                        reference_genome_menu_choice = all_menus.menu(misc, all_menus.reference_genome_menu)

                        if reference_genome_menu_choice == '':
                            misc.logfile('User input: return to previous menu')
                            break
                        # Download reference genome
                        elif reference_genome_menu_choice == '1':
                            misc.logfile('User input (setup reference genome menu): 1. Download reference genome')
                            ref_genome.download(misc, shortcuts)
                        # Index reference genome
                        elif reference_genome_menu_choice == '2':
                            misc.logfile('User input (setup reference genome menu): 1.Index reference genome')
=======
                    while True:
                        reference_genome_menu_choice = all_menus.menu(misc, all_menus.reference_genome_menu)

                        if reference_genome_menu_choice == '':
                            break
                        # Download reference genome
                        elif reference_genome_menu_choice == '1':
                            ref_genome.download(misc, shortcuts)
                        # Index reference genome
                        elif reference_genome_menu_choice == '2':
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
                            dna_analysis.index_genome_dna(misc, shortcuts)

                # Create library list file
                elif dna_menu_choice == '2':
                    misc.clear_screen()
                    dna_choice = all_menus.build_library_dna_menu(options, misc, shortcuts)

                # Run dna analysis
                elif dna_menu_choice == '3':
<<<<<<< HEAD
                    start = timeit.default_timer()
=======
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
                    misc.clear_screen()
                    misc.validate_id(options, shortcuts)
                    dna_analysis.alignment(misc, shortcuts)
                    if dna_analysis.validate_bam_dna(misc, shortcuts):
                        dna_analysis.sort(options, misc, shortcuts)
                        dna_analysis.merge(options, misc, shortcuts)
                        dna_analysis.remove_duplicate(misc, shortcuts)
                        dna_analysis.realign(misc, shortcuts)
                        dna_analysis.gatk_haplotype(options, misc, shortcuts)
                        dna_analysis.delly(options, misc, shortcuts)
                        dna_analysis.manta(misc, shortcuts)
<<<<<<< HEAD
                        end = timeit.default_timer()
                        elapsed_time = end-start
                        misc.logfile("All steps in DNA-Seq analysis pipeline completed succesfully in: {0:.1g} min \n".format(elapsed_time/60))
=======
                        print("\nAll steps in DNA-Seq analysis pipeline completed succesfully\n")
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
                        sys.exit()

        #Rna analysis menu
        elif menu_choice == '3':
<<<<<<< HEAD
            misc.logfile('User input: RNA-analysis')
=======
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510
            while True:
                misc.logfile('Program at RNA-analysis menu')
                misc.clear_screen()
                rna_choice = all_menus.menu(misc, all_menus.rna_menu)
                if rna_choice == '':
                    misc.logfile('User input: return to previous menu')
                    break

                # Index reference genome
                elif rna_choice == '1':
                    misc.logfile('User input: index reference genome')
                    while True:
                        misc.logfile('Program at reference genome index menu')
                        misc.clear_screen()
                        reference_genome_index_choice = all_menus.menu(misc, all_menus.reference_genome_index_menu)
                        if reference_genome_index_choice == '':
                            misc.logfile('User input: return to previous menu')
                            break

                        elif reference_genome_index_choice == '1':
                            misc.logfile('User input: index whole genome')
                            rna_analysis.index_genome_rna(1, None, misc, shortcuts)

                        elif reference_genome_index_choice == '2':
                            misc.logfile('User input: index parts of genome')
                            chromosomes = misc.choose_chromosomes_to_index(all_menus, shortcuts)
                            filename = misc.create_new_fasta(chromosomes, shortcuts)
                            misc.create_new_gtf(chromosomes, filename, shortcuts)
                            ref_genome.index_genome_rna(2, filename, misc, shortcuts)

#                             # If the program is closed between running index_chromosomes and map_reads, the name variable will be forgot and must be typed in again
#                             if filename is None:
#                                 filename = "".join(misc.choose_chromosomes_to_index(all_menus, shortcuts))
#                                 rna_analysis.map_reads(options, filename, misc, shortcuts)
#                                 rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
#                                 rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)
#                             else:
#                                 rna_analysis.map_reads(options, filename, misc, shortcuts)
#                                 rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
#                                 rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)
#
#                 # Index whole genome
#                 if index_reference_genome_choice == '1':
#
#                 filename = "".join(misc.choose_chromosomes_to_index(all_menus, shortcuts))
#
#
#                 # Index parts of genome
#                 elif index_reference_genome_choice == '2':
#                     chromosomes = misc.choose_chromosomes_to_index(all_menus, shortcuts)
#                     if not chromosomes:
#                         continue
#                     else:
#
#
#                         ref_genome.index_genome_rna(2, filename, misc, shortcuts)
#                         input("Press any key to return to main menu")
#
#
# rna_analysis.map_reads(options, filename, misc, shortcuts)
# rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
# rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)
#

                    #
                    # else:
                    #     input("An error was found, see the output in the shell for more information!\n\nPress any key to exit program")
                    #     sys.exit()


<<<<<<< HEAD
=======
                    # If the program is closed between running index_chromosomes and map_reads, the name variable will be forgot and must be typed in again
                    if filename is None:
                        filename = "".join(misc.choose_chromosomes_to_index(all_menus, shortcuts))
                        rna_analysis.map_reads(options, filename, misc, shortcuts)
                        rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
                        rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)
                    else:
                        rna_analysis.map_reads(options, filename, misc, shortcuts)
                        rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
                        rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)
>>>>>>> 5473b9774558377c381065f090ac58c8cb730510

                # Index whole genome
                if index_reference_genome_choice == '1':
                    ref_genome.index_genome_dna(1, None, misc, shortcuts)
                    ref_genome.index_genome_rna(1, None, misc, shortcuts)


                # Index parts of genome
                elif index_reference_genome_choice == '2':
                    chromosomes = misc.choose_chromosomes_to_index(all_menus, shortcuts)
                    if not chromosomes:
                        continue
                    else:
                        filename = misc.create_new_fasta(chromosomes, shortcuts)
                        misc.create_new_gtf(chromosomes, filename, shortcuts)
                        ref_genome.index_genome_dna(2, filename, misc, shortcuts)
                        ref_genome.index_genome_rna(2, filename, misc, shortcuts)
                        input("Press any key to return to main menu")




                    #
                    # else:
                    #     input("An error was found, see the output in the shell for more information!\n\nPress any key to exit program")
                    #     sys.exit()



if __name__ == '__main__':
    main()
