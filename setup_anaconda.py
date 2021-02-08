import subprocess
import argparse
from os import sys


############################################################# Functions ################################################################################


def main_menu():
    '''Prints an introduction to the program where you can choose what to do'''

    print("\nThis script lets you install Anaconda3 and set up an environment with the required software packages.\n")
    print("1. Download and install Anaconda3 with python 3.8 (if it's not allready installed)")
    print("2. Set up a new environment for DNA and RNA-sequence analysis\n")


    return input("(leave blank to exit): ")

def install_anaconda():
    '''This function downloads anaconda via wget and the link-adress to the linux installer from anaconda.com.
       It then installs anaconda and creates an environment with the required software packages'''

    print("Downloading and installing anaconda...\n")

    # downloads anaconda
    cmd_download = "wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -P $HOME"
    subprocess.run(cmd_download, shell=True)

    # applies read, write and execute permission for all users to the file
    cmd_chmod = "chmod a=xrw $HOME/Anaconda3-2020.11-Linux-x86_64.sh"
    subprocess.run(cmd_chmod, shell=True)

    # installs anaconda
    cmd_install = "bash $HOME/Anaconda3-2020.11-Linux-x86_64.sh"
    subprocess.run(cmd_install, shell=True)

    print("Anaconda is now successfully installed\n")
    return input("Now please close your shell and start a new shell for anaconda3 to work")

def create_environment():
    cmd_env = "conda create --name sequencing -c bioconda bwa picard gatk4 manta bcftools"
    # subprocess.run(cmd_env, shell=True)

    cmd_delly = "wget https://github.com/dellytools/delly/releases/download/v0.8.7/delly_v0.8.7_linux_x86_64bit -P $HOME/anaconda3/envs/sequencing/bin"
    subprocess.call(cmd_delly, shell=True)
    cmd_chmod = "chmod a=xrw $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit"
    subprocess.run(cmd_chmod, shell=True)
    cmd_mv = "mv $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit $HOME/anaconda3/envs/sequencing/bin/delly"
    subprocess.run(cmd_mv, shell=True)


    print("Environment is successfully created\n")
    print("Type \"conda activate sequencing\" in the shell to activate the environment before continuing to the next step\n")
    return input("Press any key to return to main menu...")

# ----------------------------------------------------- Main program starts -----------------------------------------------------------------
def main():


    # metavar makes the help text more tidy
    parser = argparse.ArgumentParser(description='''This script is used to simplify installation and set up of the python 3.8 distribution Anaconda.\n Please read the README file first!''')
    options = parser.parse_args() # all arguments will be passed to the functions

    subprocess.run("clear", shell=True)
    menu_choice = main_menu()

    while menu_choice != "":
        subprocess.run("clear", shell=True)
        if menu_choice == '1':
            install_anaconda()
            menu_choice = main_menu()

        elif menu_choice == '2':
            create_environment()
            menu_choice = main_menu()


        else:
            print("Invalid choice, please try again!")
            menu_choice = main_menu()

    print("exiting program...")
    sys.exit()







if __name__ == '__main__':
    main()
