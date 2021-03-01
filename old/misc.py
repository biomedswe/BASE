from os import listdir, sys, mkdir, makedirs, getenv, kill, getppid, path
import subprocess
import signal
import time
import shlex


#################################################### MISCELLANEOUS FUNCTIONS #################################################

def confirm_choice():
    while True:
        choice = input("Are you sure? y or n: ")
        if choice.lower() == "y":
            return True
        elif choice.lower() == 'n':
            return False
        else:
            print("Invalid choice, type 'n' or 'y'!")

def clear_screen():
    subprocess.run("clear", shell=True)

def close_terminal():
    kill(getppid(), signal.SIGHUP)

def info_script():
    print("This program was created by biomedicine master student Jonas Andersson during a masterproject 2020-2021.")
    print("For correspondence please contact jonas870318@gmail.com\n")
    print("Biomedswe 2021 v.1.0")
    print("Automated DNA/RNA-sequencing script\n")

def create_dir():
    try:
        # create directory three
        makedirs(getenv("HOME")+"/sequencing_project/dna_seq/reads")
        mkdir(getenv("HOME")+"/sequencing_project/reference_genome")
        makedirs(getenv("HOME")+"/sequencing_project/rna_seq/reads")

    except FileExistsError:
        pass

def validate_choice(choices,text):

    while True:
        choice = input(text)

        if choice == "":
            return ""

        elif choice not in [str(i) for i in range(1,choices+1)]:
            print("Invalid choice, try again!")

        else:
            return choice

def run_command(command, step):
    return_code = subprocess.run(command, shell=True)
    if return_code.returncode == 0:
        print(f"{step} completed without errors!, continuing with next step...")
    else:
        print('\nAn error has occured, see shell for more information. Exiting program...')
        sys.exit()

def create_outputList(output_path, write_to_file):
    with open(output_path, 'a') as c:
        c.write(f'{write_to_file}\n')

def create_directory(path):
    try:
        # create target directory
        mkdir(path)
    except FileExistsError:
        pass

def step_completed(file, text):
    if path.isfile(file):
        print(f"{text}, allready completed, skips step...\n")
        time.sleep(2)
        return True

def setup_performance():
    pass
