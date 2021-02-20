import subprocess
import argparse
from os import sys
from misc import *



############################################################# Functions ################################################################################


def anaconda_menu():
    '''Prints an introduction to the program where you can choose what to do'''

    info_script()
    print("This script lets you install Anaconda3 and set up an environment with the required software packages.\n")
    print("If you already have installed anaconda, skip first step and proceed with setting up an environment\n")
    print("1. Download and install Anaconda3 with python 3.8")
    print("2. Set up a new conda environment for DNA and RNA-sequence analysis\n")
    choices = 2
    text = "(leave blank to exit)"
    choice = validate_choice(choices, text)
    return choice



def install_anaconda():
    '''This function downloads anaconda via wget and the link-adress to the linux installer from anaconda.com.
       It then installs anaconda and creates an environment with the required software packages.
       Finally creates the directory three in the sequencing_project folder'''

    clear_screen()
    print("""Download and install Anaconda3 with python 3.8.


Downloading and installing anaconda from https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh...""")

    # downloads anaconda
    cmd_download = "wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -P $HOME"
    # subprocess.run(cmd_download, shell=True)

    # applies read, write and execute permission for all users to the file
    cmd_chmod = "chmod a=xrw $HOME/Anaconda3-2020.11-Linux-x86_64.sh"
    # subprocess.run(cmd_chmod, shell=True)

    # installs anaconda
    cmd_install = "bash $HOME/Anaconda3-2020.11-Linux-x86_64.sh"
    # subprocess.run(cmd_install, shell=True)

    # create directory three
    create_dir()

    print("\nAnaconda is now successfully installed\n\n")
    input("The terminal must be restarted for anaconda3 to initialize correctly\n\nPress any key to close the terminal\nThen start it again manually")
    close_terminal()

def create_anaconda_environment():
    clear_screen()

    print("""Set up a new conda environment for DNA and RNA-sequence analysis.


Creating a new environment named \"Secuencing\"...

Packages to be installed:
Burrows-Wheeler aligner (bwa): (aligns reads to reference genome)
Picard tools,
Genome analysis tool kit 4 (gatk),
Bcftools
Samtools
Manta (structural variant caller)
Delly (structural variant caller)""")

    cmd_env = "conda create --name sequencing -c bioconda bwa picard gatk4 manta bcftools samtools star"
    subprocess.run(cmd_env, shell=True)

    # delly in bioconda didn't work so I had to do a workaround
    cmd_delly = "wget https://github.com/dellytools/delly/releases/download/v0.8.7/delly_v0.8.7_linux_x86_64bit -P $HOME/anaconda3/envs/sequencing/bin"
    subprocess.call(cmd_delly, shell=True)
    cmd_chmod = "chmod a=xrw $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit"
    subprocess.run(cmd_chmod, shell=True)
    cmd_mv = "mv $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit $HOME/anaconda3/envs/sequencing/bin/delly"
    subprocess.run(cmd_mv, shell=True)


    print("\n\n\nEnvironment named \"sequencing\" was successfully created\n")
    print("You must manually activate the created environment.\nTo do this, type \"conda activate sequencing\" in the shell before proceeding to the next step\n")
    return input("Press any key to return to main menu...")

# ----------------------------------------------------- Main program starts -----------------------------------------------------------------
def main():


    # metavar makes the help text more tidy
    parser = argparse.ArgumentParser(description='''This script is used to simplify installation and set up of the python 3.8 distribution Anaconda.\n Please read the README file first!''')
    options = parser.parse_args() # all arguments will be passed to the functions

    subprocess.run("clear", shell=True)
    menu_choice = anaconda_menu()

    while menu_choice != "":
        subprocess.run("clear", shell=True)
        if menu_choice == '1':
            install_anaconda()
            menu_choice = anaconda_menu()

        elif menu_choice == '2':
            create_environment()
            menu_choice = anaconda_menu()


        else:
            print("Invalid choice, please try again!")
            menu_choice = anaconda_menu()

    print("exiting program...")
    sys.exit()

if __name__ == '__main__':
    main()
