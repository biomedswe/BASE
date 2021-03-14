# Packages used in script
from os import getenv, sys, path, listdir, makedirs
# import subprocess
import argparse
from menus import Menus, Misc, Shortcuts
import rna_seq_analysis
import dna_seq_analysis
from reference_genome import ReferenceGenome
import setup_anaconda3
import time





# Shortcuts

############################################################# Functions ################################################################################








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
    rna_analysis = rna_seq_analysis.RnaSeqAnalysis()
    dna_analysis = dna_seq_analysis.DnaSeqAnalysis()
    ref_genome = ReferenceGenome()
    setup = setup_anaconda3.SetupAnaconda3()
    # main menu is shown untill choice is left blank
    while True:
        # Clears screen and prints main menu
        menu_choice = all_menus.menu(misc, all_menus.main_menu)

        if menu_choice == "":
            print("exiting program...")
            break

        # Setyp anaconda3 menu
        if menu_choice == '1':
            anaconda_choice = all_menus.menu(misc, all_menus.anaconda_menu)

            if anaconda_choice == '1':
                setup.install_anaconda(misc, shortcuts)

            elif anaconda_choice == '2':
                setup.create_anaconda_environment(misc, shortcuts)





        # Download reference genome
        elif menu_choice == '2':
            reference_genome_menu_choice = all_menus.menu(misc, all_menus.reference_genome_menu)
            if reference_genome_menu_choice == '1':
                download_ref_genome.download(misc, shortcuts)
            elif reference_genome_menu_choice == '2':
                index_reference_genome_choice = all_menus.menu(misc, all_menus.reference_genome_index_menu)

                # Index whole genome
                if index_reference_genome_choice == '1':
                    if misc.step_completed(shortcuts.index_reference_genome_complete, 'Indexing whole reference genome allready completed, returning to main menu...'):
                        time.sleep(2.5)
                        continue
                    else:
                        ref_genome.index_genome_dna(1, None, misc, shortcuts)
                        ref_genome.index_genome_rna(1, None, misc, shortcuts)


                # Index parts of genome
                elif index_reference_genome_choice == '2':
                    chromosomes = misc.choose_chromosomes_to_index(misc, shortcuts)
                    if not chromosomes:
                        continue
                    else:
                        filename = misc.create_new_fasta(chromosomes, shortcuts)
                        misc.create_new_gtf(chromosomes, filename, misc, shortcuts)
                        dna_analysis.index_reference_dna(2, filename, misc, shortcuts)
                        rna_analysis.index_genome_rna_analysis(2, filename, misc, shortcuts)
                        input("Press any key to continue")

        # Dna analysis menu
        elif menu_choice == '3':
            while True:
                misc.clear_screen()
                dna_choice = all_menus.menu(misc, all_menus.dna_menu)

                if dna_choice == '':
                    break

                elif dna_choice == '1':
                    while True:
                        misc.clear_screen()
                        dna_choice = all_menus.build_library_dna_menu(options, misc, shortcuts)
                        if dna_choice == "":
                            break


                elif dna_choice == '2':
                    misc.clear_screen()
                    misc.validate_id(options, shortcuts)
                    dna_choice = dna_analysis.alignment(misc, shortcuts)
                    if dna_choice == "":
                        continue
                    else:
                        if dna_analysis.validate_bam_dna(misc, shortcuts):
                            dna_analysis.sort(options, misc, shortcuts)
                            dna_analysis.merge(options, misc, shortcuts)
                            dna_analysis.remove_duplicate(misc, shortcuts)
                            dna_analysis.realign(misc, shortcuts)
                            dna_analysis.gatk_haplotype(options, misc, shortcuts)
                            dna_analysis.delly(options, misc, shortcuts)
                            dna_analysis.manta(misc, shortcuts)
                            print("\nAll steps in DNA-Seq analysis pipeline completed succesfully\n")
                            sys.exit()


                        else:
                            input("An error was found, see the output in the shell for more information!\n\nPress any key to exit program")
                            sys.exit()

        # Rna analysis menu
        elif menu_choice == '4':
            while True:
                misc.clear_screen()
                rna_choice = all_menus.menu(misc, all_menus.rna_menu)
                if rna_choice == '':
                    break

                # Map reads to genome
                elif rna_choice == '1':
                    misc.clear_screen()

                    # If the program is closed between running index_chromosomes and map_reads, the name variable will be forgot and must be typed in again
                    if filename is None:
                        filename = "".join(misc.choose_chromosomes_to_index(misc, shortcuts))
                        rna_analysis.map_reads(options, filename, misc, shortcuts)
                        rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
                        rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)
                    else:
                        rna_analysis.map_reads(options, filename, misc, shortcuts)
                        rna_analysis.ASEReadCounter(options, filename, misc, shortcuts)
                        rna_analysis.add_wgs_data_to_csv(options, filename, misc, shortcuts)

if __name__ == '__main__':
    main()
