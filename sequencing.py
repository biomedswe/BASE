import subprocess
import argparse
from setup_anaconda3 import *
from dna_seq import *
from misc import *
from download_reference_genome import *


############################################################# Functions ################################################################################










def main_menu():
    '''Prints an introduction to the program where you can choose what to do'''

    while True:
        info_script()
        print("1. Anaconda3 setup\n")
        print("2. Download reference genome\n")
        print("3. DNA-analysis\n")
        print("4. RNA-analysis\n")


        choice = input("(leave blank to exit): ")
        if choice == "":
            return ""

        elif choice not in ['1','2','3','4']:
            print("Invalid choice, try again!")
            time.sleep(1)
            clear_screen()

        else:
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

        if menu_choice == '1':
            clear_screen()
            anaconda_choice = anaconda_menu() # Shows anaconda menu and returns choice

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







    print("exiting program...")
    sys.exit()







if __name__ == '__main__':
    main()
