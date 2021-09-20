from os import path, getenv, listdir, makedirs, sys, remove, kill, getppid
import signal
from subprocess import Popen, PIPE, STDOUT, run
import logging
import timeit
import time
import shlex
import re
logging.basicConfig(filename = getenv("HOME")+'/BASE/Logfile.txt',
                    format = '%(levelname)s     %(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                    level = logging.DEBUG,
                    filemode = 'a')

try:
    from Bio import SeqIO
except:
    pass

class Misc():
    '''This class contains miscellaneous functions related to general functionality'''

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def clear_screen(self):
        '''This function clears the terminal window'''

        run("clear", shell=True)
        

    #---------------------------------------------------------------------------
    def close_terminal(self):
        '''This function closes the terminal'''
        
        kill(getppid(), signal.SIGHUP)
        
    #---------------------------------------------------------------------------
    def confirm_choice(self):
        '''This function asks if you want to confirm your choice'''

        try:
            while True:
                choice = input("Are you sure? y or n: ")
                if choice.lower() == "y":
                    return True
                elif choice.lower() == 'n':
                    return False
                else:
                    print("Invalid choice, type 'n' or 'y'!")
        except Exception as e:
            logging.error("Error with confirm_choice() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def create_directory(self, paths):
        '''This function creates the output directories for the different analysis steps'''

        for path in paths:
            try:
                # create target directory
                makedirs(path)
                self.log_to_file("INFO", f"{path} folder(s) created succesfully - OK!")
            except FileExistsError:
                self.log_to_file("INFO", f"{path} folder(s) allready exists - skips step...")
            except Exception as e:
                logging.error("Error with create_directory() in miscellaneous.py: {0}. Exiting program...".format(e))
                sys.exit()

    #---------------------------------------------------------------------------
    def create_outputList_dna(self, output_path, write_to_file):
        '''This function creates a list file containing the output name of the files created in the pipeline step where the function is used, the list is then used in the next step'''

        try:
            with open(output_path, 'a') as c:
                c.write(f"{write_to_file}\n")

        except Exception as e:
            logging.error("Error with create_outputlist_dna() in miscellaneous.py: {0}. Exiting program...".format(e))
    #---------------------------------------------------------------------------
    def create_trackFile(self, file):
        '''This function creates a trackfile that step_allready_completed() function can look after when checking if step is allready completed'''
        try:
            if path.isfile(file):
                self.log_to_file("WARNING", f"You are trying to overwrite the existing file: {file}")
                if not self.confirm_choice(): 
                    return
                    
            with open(file, 'w'):
                pass
        except Exception as e:
            logging.error("Error with create_trackfile() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def elapsed_time(self, elapsed):
        try:
            if elapsed >= 86400:
                return f"{round(elapsed/86400, 2)} days"
            elif  3600 <= elapsed < 86400:
                return f"{round(elapsed/3600, 2)} hours"
            elif 60 <= elapsed < 3600:
                return f"{round(elapsed/60, 2)} minutes"
            elif elapsed < 60:
                return f"{round(elapsed, 2)} seconds"
            else:
                return 'error'
        except Exception as e:
            logging.error("Error with elapsed_time() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def log_to_file(self, level, text):
        '''This function handles the logging'''
        # .format and not f-strings are used here in order to give python2 compatibility


        try:
            print("{0}: {1}\n".format(level, text))
            if level == "DEBUG": logging.debug("{0}".format(text))
            elif level == "INFO": logging.info("{0}".format(text))
            elif level == "WARNING": logging.warning("{0}".format(text))
            elif level == "ERROR": logging.error("{0}".format(text))
            elif level == "CRITICAL": logging.critical("{0}".format(text))

        except Exception as e:
            logging.error("Error with {0}.log_to_file() in miscellaneous.py: {1}. Exiting program...".format(self, e))
            sys.exit()

    #---------------------------------------------------------------------------
    def remove_file(self, file):
        '''This function removes incomplete files if the processing of file ended with returncode != 0'''
        try:
            if file:
                remove(file)
                self.log_to_file("INFO", f'Incomplete {file} removed - OK!')
        except OSError: pass
        except Exception as e:
            logging.error("Error with remove_file() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def run_command(self, input):
        '''This function first calls step_allready_completed() to check if the step i allready completed.
        If not completed; if process is executed without errors, it prints to logfile with time taken and passes return.
        else; if process ends with errors, it prints error to log, removes incomplete file and then exits program'''

        try:
            
            cmd = input.get('cmd')
            program = input.get('program') if input.get('program') else None
            text = input.get('text')
            trackfile = input.get('file')
            old_file = input.get('old_file') if input.get('old_file') else None
            
            # self.log_to_file("DEBUG", f"run_command(cmd: {cmd},\ntext: {text},\nfile: {trackfile})")


            # Check if step is allready completed
            if self.step_allready_completed(trackfile, text):
                return True

            # Invoke process for easier commands during setup
            if input.get('program') == "setup":
                process = run(cmd, shell=True)
                return True

            # invoke process
            # bwa-mem2 requires set -o pipefail which in turn requires executable='/bin/bash'
            else:
                process = Popen(cmd, executable='/bin/bash', shell=True, stdout=PIPE, stderr=STDOUT, text=True)

            # print stdout during execution
            while process.poll() is None:
                output = process.stdout.readline()
                
                if output:
                    # For some reason, output is empty after process has terminated and any errors can't be logged to file
                    # this workaround (latest = output) solves this issue
                    latest = output
                    print(output.strip())

                    if "No errors found" in output:
                        self.log_to_file("INFO", f"{output.strip()} in {text.split()[1]}")

               

                              
            # Checks the returncode after process has terminated        
            rc = process.poll()
            
            # Do this if process ended with returncode = 0
            if rc == 0:
                
                # If trackfile is specified, create a ".completed" file telling python that this step is allready completed
                if trackfile: self.create_trackFile(trackfile)
                if old_file: 
                    run(f"rm -f {old_file}", shell=True)
                    self.log_to_file("INFO", f"{old_file} was removed to save space")
                return True
           
            # If process didn't end with returncode = 0
            else:
                raise Exception(f"Error message: {latest.strip()}")
                
                
                

            # option2 without shell = True

                # if cmd == "bwa-mem2":
                #     p1 = Popen(shlex.split(cmd_1), stdout=PIPE, stderr=PIPE, text=True)
                #     p2 = Popen(shlex.split(cmd_2), stdin=p1.stdout, stdout=PIPE, stderr=PIPE, text=True)
                #     p1.stdout.close()
                #     output, error = p2.communicate() # let p2 finish running
                #     p1.wait()                        # ensure p1 has properly exited

                #     rc_1 = p1.poll()
                #     rc_2 = p2.poll()
                #     print("rc_1:", rc_1)
                #     print("rc_2:", rc_2)

                #     if rc_1 == 0 and rc_2 == 0:
                #         self.log_to_file("DEBUG", "# Process ended with returncode = 0")
                #         if text: self.log_to_file("INFO", f"{text} succesfully completed")
                    
                #         # If trackfile is specified, create a ".completed" file telling python that this step is allready completed
                #         if trackfile: self.create_trackFile(trackfile)
                #         return True
            
                #     # If process didn't end with returncode = 0
                #     else:
                #         print("Raise exception")
                #         raise Exception(f"stderr1: {p1.stderr.readline()}, stdout2: {output}, stderr2: {error}")

 

            
        
                

        except Exception as e:
            self.log_to_file("ERROR", f"{input.get('program')} Process ended with returncode != 0, {e}")
            sys.exit()

    #---------------------------------------------------------------------------
    def step_allready_completed(self, file, text):
        '''This function checks if a step is allready completed by checking if "file" allready exists.
        If file exists: returns True, else return False'''

        try:
            if file:
                if path.isfile(file):
                    if text: self.log_to_file("INFO", f"{text} allready completed, skips step...")
                    time.sleep(2)
                    return True
                else:
                    return False
            else:
                return False
        except Exception as e:
            logging.error("Error with step_allready_completed() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def validate_choice(self, choices, text):
        '''This function checks if the input choice is valid'''

        try:
            while True:
                choice = input(text)
                if choice == "":
                    return ""
                elif choice not in [str(i) for i in range(1,choices+1)]:
                    print("Invalid choice, try again!")
                else:
                    return choice
        except Exception as e:
            logging.error("Error with validate_choice() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def validate_id(self,options, shortcuts):
        '''This function validates wether the tumor_id and normal_id you entered is present in your reads'''

        try:
            if options.tumor_id and options.normal_id in "".join(listdir(shortcuts.dna_reads_dir)):
                self.log_to_file("INFO", f"tumor_id {options.tumor_id} and normal_id {options.normal_id} correctly validated")

            else:
                self.log_to_file("ERROR", f'You have entered a tumor_id: {options.tumor_id} and normal_id: {options.normal_id} that is not present in your reads.\nPlease restart program and verify that you have typed in the right \"clinical_id\" for both tumor (-t) and normal (-n)!')
                input("Press any key to exit program")
                sys.exit()

        except OSError as e:
            self.log_to_file("ERROR", f"{e} No read files. You have to save you reads in BASE/dna_seq/reads first")
            sys.exit()
        except Exception as e:
            logging.error("Error with validate_id() in miscellaneous.py: {0}. Exiting program...".format(e))
            sys.exit()
