# Supports python 2 and >
from os import sys, path, getenv
import subprocess
import logging
logging.basicConfig(filename='Logfile.txt', filemode='a', format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)

class SetupAnaconda3():

    def __init__(self):
        pass

    #---------------------------------------------------------------------------
    def log_to_file(self, text):
        try:
            print('\n{0}'.format(text))
            logging.info('\n{0}'.format(text))
        except Exception as e:
            logging.info('Error with self.log_to_file() in setup_anaconda3.py: {0}. Exiting program...'.format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def run_command(self, command, text, file, trackfile):
        '''This function first calls _step_allready_completed() to check if the step i allready completed.
        If not completed; if process is executed without errors, it prints to logfile with time taken and passes return.
        else; if process ends  with errors, it prints error to log, removes incomplete file and then exits program'''

        try:
            if not self._step_allready_completed(file, text):
                return_code = subprocess.call(command, shell=True)
                if return_code == 0:
                    self.log_to_file("{0} succesfully completed - OK!".format(text))
                    if trackfile:
                        self._create_trackFile(trackfile)
                        return
                else:
                    self.log_to_file('Process ended with returncode != 0, {0} - ERROR!'.format(command))
                    subprocess.call("rm {0}".format(file), shell=True)
                    self.log_to_file('Incomplete {0} removed - OK!'.format(file))
                    sys.exit()
        except Exception as e:
            self.log_to_file('Error with SetupAnaconda3.run_command() in setup_anaconda3.py: {0}. Exiting program...'.format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def _create_trackFile(self, file):
        '''This function creates a trackfile that _step_allready_completed() function can look after when checking if step is allready completed'''
        try:
            with open(file, 'w'):
                self.log_to_file('Trackfile {0} created - OK!'.format(file))
        except Exception as e:
            self.log_to_file('Error with SetupAnaconda3._create_trackFile() in setup_anaconda3.py: {0}'.format(e))

    #---------------------------------------------------------------------------
    def _step_allready_completed(self, file, text):
        '''This function checks if a step is allready completed by checking if "file" allready exists, if so, returns True, else return False'''
        try:
            if file:
                if path.isfile(file):
                    self.log_to_file("{0} allready completed, skips step...".format(text))
                    return True
                else:
                    return False
            else:
                False
        except Exception as e:
            self.log_to_file('Error with SetupAnaconda3._step_allready_completed() in setup_anaconda3.py: {0}. Exiting program...'.format(e))
            sys.exit()

    #---------------------------------------------------------------------------
    def install_anaconda(self):
        '''This function downloads anaconda via wget and the link-adress to the linux installer from anaconda.com.
           It then installs anaconda and removes the installation file after completion.'''

        self.log_to_file("Downloading anaconda from https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh...")
        # Download Anaconda3
        self.run_command("wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -P $HOME", "Downloading Anaconda3-2020.11-Linux-x86_64.sh", getenv("HOME")+'/anaconda3/install.complete', None)
        # Install Anaconda3
        self.run_command("bash $HOME/Anaconda3-2020.11-Linux-x86_64.sh", "Installation of Anaconda3", getenv("HOME")+'/anaconda3/install.complete', None)
        # remove installation file
        self.run_command("rm $HOME/Anaconda3-2020.11-Linux-x86_64.sh", 'Removal of installation file', getenv("HOME")+'/anaconda3/install.complete', getenv("HOME")+'/anaconda3/install.complete')
        print("\nInstallation of Anaconda3 is now completed,\nBefore you are done, you have to Copy and paste the following in the terminal: \"source $HOME/.bashrc\" to initialize anaconda3")

    #---------------------------------------------------------------------------
    def create_anaconda_environment(self, misc, shortcuts):
        '''This function setups a new anaconda environment with all the packages required for the program to run'''

        try:
            # create directory three
            misc.create_directory([
            shortcuts.aligned_output_dir,
            shortcuts.sorted_output_dir,
            shortcuts.merged_output_dir,
            shortcuts.removed_duplicates_output_dir,
            shortcuts.realigned_output_dir,
            shortcuts.haplotypecaller_output_dir,
            shortcuts.delly_output_dir,
            shortcuts.manta_output_dir,
            shortcuts.rna_reads_dir,
            shortcuts.dna_reads_dir,
            shortcuts.star_output_dir,
            shortcuts.GRCh38_dir,
            shortcuts.star_index_dir_whole_genome,
            shortcuts.GRCh38_chunks_dir
            ], "Directory tree")

            if path.isdir(getenv("HOME")+'/anaconda3/envs/sequencing'):
                misc.log_to_file("Environment allready exists, do you want to owerwrite?")
                if not misc.confirm_choice():
                    misc.log_to_file('User input: No')
                    return

            misc.clear_screen()
            misc.log_to_file('User input: Yes')
            misc.log_to_file('''Set up a new conda environment for DNA and RNA-sequence analysis.

Creating a new environment named \"Secuencing\"...

Packages to be installed:
Bedtools
Bcftools
Biopython
Burrows-Wheeler aligner (bwa)
Delly (structural variant caller)
Genome analysis tool kit 4 (gatk)
Manta (structural variant caller)
Pandas
Picard tools
Python 3.7.6
Samtools 1.9
Scipy
Snpeff
STAR
Vcfpy
''')

            cmd_env = "conda create -n sequencing -c bioconda bedtools bcftools biopython bwa gatk4 picard python=3.7.6 samtools=1.9 star pandas vcfpy scipy snpeff"
            misc.run_command(cmd_env, "Installing bedtools bcftools biopython bwa gatk4 picard python=3.7.6 samtools=1.9 star pandas vcfpy scipy snpeff", None, None)

            # Delly in bioconda didn't work so I had to do a workaround
            cmd_download_delly = "wget https://github.com/dellytools/delly/releases/download/v0.8.7/delly_v0.8.7_linux_x86_64bit -P $HOME/anaconda3/envs/sequencing/bin"
            misc.run_command(cmd_download_delly, 'Downloading delly', None, None)

            cmd_chmod_delly = "chmod a=xrw $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit"
            misc.run_command(cmd_chmod_delly, 'Granting permission to move delly', None, None)
            cmd_mv = "mv $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit $HOME/anaconda3/envs/sequencing/bin/delly"
            misc.run_command(cmd_mv, "Moving Delly to $HOME/anaconda3/envs/sequencing/bin/delly", None, None)

            # Manta in bioconda didn't work so I had to do a workaround
            cmd_download_manta = "wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 -P $HOME/anaconda3/envs/sequencing/bin"
            misc.run_command(cmd_download_manta, "Downloading Manta", None, None)
            cmd_unpack = "tar -xf $HOME/anaconda3/envs/sequencing/bin/manta-1.6.0.centos6_x86_64.tar.bz2 -C $HOME/anaconda3/envs/sequencing/bin/"
            misc.run_command(cmd_unpack, "Unpacking Manta", None, None)



            print("\n\n\nEnvironment named \"sequencing\" was successfully created\n")
            print("You must manually activate the created environment.\nTo do this, type \"conda activate sequencing\" in the shell before restarting the program\n")
            input("Press any key to exit program...")
            sys.exit()
        except Exception as e:
            misc.log_to_file('Error with create_anaconda_environment(): {0}'.format(e))


def main():
    setup = SetupAnaconda3()
    setup.install_anaconda()
if __name__ == '__main__':
    main()
