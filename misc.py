import subprocess
import os
import signal
import time
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
    os.kill(os.getppid(), signal.SIGHUP)

def info_script():
    print("This program was created by biomedicine master student Jonas Andersson during a masterproject 2020-2021.")
    print("For correspondence please contact jonas870318@gmail.com\n")
    print("Biomedswe 2021 v.1.0")
    print("Automated DNA/RNA-sequencing script\n")
