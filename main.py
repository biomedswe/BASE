
# coding=utf-8
# Packages used in script
from os import getenv, sys, path, listdir, makedirs
import argparse
from menus import Menus, Misc, Shortcuts
from rna_seq_analysis import RnaSeqAnalysis
from dna_seq_analysis import DnaSeqAnalysis
from reference_genome import ReferenceGenome
from setup_anaconda3 import SetupAnaconda3
import time
import timeit


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

    misc.log_to_file('-----Program starts-----')
    while True:
        # Main menu
        menu_choice = all_menus.menu(misc, all_menus.main_menu)
        misc.log_to_file('Program at main menu')
        if menu_choice == "":
            misc.log_to_file("User input: exit program")
            break

        # Setyp anaconda3 environment
        if menu_choice == '1':
            misc.log_to_file('User input: 1, Setup anaconda3 environment')
            setup.create_anaconda_environment(misc, shortcuts)

        # Dna analysis menu
        elif menu_choice == '2':
            misc.log_to_file('User input: 2. DNA-analysis')
            while True:
                misc.log_to_file('Program at DNA-analysis menu')
                dna_menu_choice = all_menus.menu(misc, all_menus.dna_menu)
                if dna_menu_choice == '':
                    misc.log_to_file('User input: return to main menu')
                    break

                # Setup reference genome
                elif dna_menu_choice == '1':
                    misc.log_to_file('User input: 1. Setup reference genome ')
                    while True:
                        misc.log_to_file('Program at Setup reference genome menu')
                        reference_genome_menu_choice = all_menus.menu(misc, all_menus.reference_genome_menu)

                        if reference_genome_menu_choice == '':
                            misc.log_to_file('User input: return to DNA-analysis menu')
                            break
                        # Download reference genome
                        elif reference_genome_menu_choice == '1':
                            misc.log_to_file('User input: 1. Download reference genome')
                            ref_genome.download(misc, shortcuts)
                        # Index reference genome
                        elif reference_genome_menu_choice == '2':
                            misc.log_to_file('User input: 2. Index reference genome')
                            dna_analysis.index_genome_dna(misc, shortcuts)

                # Create library list file
                elif dna_menu_choice == '2':
                    misc.log_to_file('User input: 2. Create library list file')
                    misc.clear_screen()
                    all_menus.build_library_dna_menu(options, misc, shortcuts)

                # Run dna analysis
                elif dna_menu_choice == '3':
                    misc.log_to_file('User input: 3. Run analysis')
                    start = timeit.default_timer()
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
                        end = timeit.default_timer()
                        elapsed_time = end-start
                        misc.log_to_file(f"All steps in DNA-Seq analysis pipeline completed succesfully in: {elapsed_time/60:.1g} min")
                        sys.exit()

        # Rna analysis menu
        elif menu_choice == '3':
            misc.log_to_file('User input: RNA-analysis')
            while True:
                misc.log_to_file('Program at RNA-analysis menu')
                misc.clear_screen()
                rna_choice = all_menus.menu(misc, all_menus.rna_menu)
                if rna_choice == '':
                    misc.log_to_file('User input: return to previous menu')
                    break

                # Index reference genome
                elif rna_choice == '1':
                    misc.log_to_file('User input: index reference genome')
                    while True:
                        misc.clear_screen()
                        index_reference_genome_choice = all_menus.menu(misc, all_menus.reference_genome_index_menu)
                        misc.log_to_file('Program at reference genome index menu')
                        if index_reference_genome_choice == '':
                            misc.log_to_file('User input: return to previous menu')
                            break

                        # Index whole genome
                        elif index_reference_genome_choice == '1':
                            misc.log_to_file('User input: Index whole genome')
                            rna_analysis.index_genome_rna(1, None, misc, shortcuts)


                        # Index parts of genome
                        elif index_reference_genome_choice == '2':
                            misc.log_to_file('User input: Index parts of genome')
                            chromosomes = misc.choose_chromosomes_to_index(all_menus, shortcuts)
                            if not chromosomes:
                                continue
                            else:
                                filename = misc.create_new_fasta(chromosomes, shortcuts)
                                misc.create_new_gtf(chromosomes, filename, shortcuts)
                                ref_genome.index_genome_dna(2, filename, misc, shortcuts)
                                ref_genome.index_genome_rna(2, filename, misc, shortcuts)
                                input("Press any key to return to main menu")

                elif rna_choice == '2':
                    misc.log_to_file('User input: Map reads to reference genome')
                    map_reads_choice = all_menus.menu(misc, all_menus.map_reads_menu)
                    while True:
                        misc.log_to_file('Program at Map reads to reference genome menu')
                        reference_genome_menu_choice = all_menus.menu(misc, all_menus.reference_genome_menu)

                        if reference_genome_menu_choice == '':
                            misc.log_to_file('User input: return to previous menu')
                            break

                        # Map to whole genome
                        elif reference_genome_menu_choice == '1':
                            misc.log_to_file('User input: Map reads to whole genome')


                        # Map to parts of genome
                        elif reference_genome_menu_choice == '2':
                            misc.log_to_file('User input (setup reference genome menu): 1.Index reference genome')





if __name__ == '__main__':
    main()
