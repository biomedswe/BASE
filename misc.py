from os import listdir, sys, mkdir, makedirs, getenv, kill, getppid
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
    return_code = 0 #subprocess.run(command, shell=True)
    if return_code == 0:
        print(f"{step} completed without errors!, continuing with next step...")
    else:
        print('An error has occured, see shell for more information. Exiting program...')
        sys.exit()

def create_outputList(output_path, write_to_file):
    with open(output_path, 'a') as c:
        c.write(f'{write_to_file}\n')




# def completed_steps(output):
#     with open('dna_seq/completed_steps.txt', 'a') as c:
#         c.write(output)
#         print("Completed! saved to completed_steps.txt.")

def setup_performance():
    pass







# def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
#     """
#     Call in a loop to create terminal progress bar
#     @params:
#         iteration   - Required  : current iteration (Int)
#         total       - Required  : total iterations (Int)
#         prefix      - Optional  : prefix string (Str)
#         suffix      - Optional  : suffix string (Str)
#         decimals    - Optional  : positive number of decimals in percent complete (Int)
#         length      - Optional  : character length of bar (Int)
#         fill        - Optional  : bar fill character (Str)
#         printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
#     """
#     percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
#     filledLength = int(length * iteration // total)
#     bar = fill * filledLength + '-' * (length - filledLength)
#     print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
#     # Print New Line on Complete
#     if iteration == total:
#         print()
#
# import time
#
# # A List of Items
# items = list(range(0, 57))
# l = len(items)
#
# # Initial call to print 0% progress
# printProgressBar(0, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
# for i, item in enumerate(items):
#     # Do stuff...
#     time.sleep(1)
#     # Update Progress Bar
#     printProgressBar(i + 1, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
