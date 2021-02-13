import subprocess
import argparse
from setup_anaconda3 import *
from dna_seq import *
from misc import *
from download_reference_genome import *
from os import sys


############################################################# Functions ################################################################################










def main_menu():
    '''Prints an introduction to the program where you can choose what to do'''

    info_script()
    print("1. Anaconda3 setup\n")
    print("2. Download reference genome\n")
    print("3. DNA-analysis\n")
    print("4. RNA-analysis\n")
    choices = 4
    text = "(leave blank to exit)"
    choice = validate_choice(choices,text)
    return choice




# metavar makes the help text more tidy
parser = argparse.ArgumentParser(description='''This script is used to simplify installation and set up of Anaconda, download reference genome and perform DNA and RNA analysis. Please read the README file first!''')
parser.add_argument("-t", "--tumor_id", metavar="", required=True, help="Input clinical id of tumor samples")
parser.add_argument("-n", "--normal_id", metavar="", required=True, help="Input clinical id of normal samples")
parser.add_argument("-i", "--intervals", metavar="", required=False, help="Input path to reference-genome interval if you have any (for use in SNV calling)")
options = parser.parse_args() # all arguments will be passed to the functions


def main():



    # main menu is shown untill choice is left blank
    while True:
        clear_screen()
        menu_choice = main_menu() # Shows main menu and returns choice

        if menu_choice == "":
            break

        # Anaconda menu option
        if menu_choice == '1':
            clear_screen()
            anaconda_choice = anaconda_menu()

            if anaconda_choice == '1':
                install_anaconda()

            elif anaconda_choice == '2':
                create_anaconda_environment()
                continue

            elif anaconda_choice == '':
                continue



        # reference genome option
        elif menu_choice == '2':
            download_ref()
        elif menu_choice == '3':
            while True:
                clear_screen()
                dna_choice = dna_menu()

                if dna_choice == '':
                    break

                if dna_choice == '1':
                    clear_screen()
                    index_reference()

                elif dna_choice == '2':
                    while True:
                        clear_screen()
                        dna_choice = build_library(options)
                        if dna_choice == "":
                            break


                elif dna_choice == '3':
                    clear_screen()
                    validate_id(options)
                    align_list = alignment()
                    if validate_bam(align_list):
                        print("Validation completed without errors")
                        time.sleep(1)
                        sort_list = sort(options, align_list)
                        merge_list = merge(sort_list, options)
                        rd_list = remove_duplicate(merge_list)
                        realign_output = realign(rd_list)
                        snv_calling(options, realign_output)
                        delly(options, realign_output)
                        print("All completed succesfully")
                    else:
                        input("An error was found, see the output in the shell for more information!\n\nPress any key to exit program")
                        sys.exit()


        else:
            continue





    print("exiting program...")
    sys.exit()







if __name__ == '__main__':
    main()
