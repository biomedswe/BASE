from os import sys, path, getenv
import time

class SetupAnaconda3():

    def __init__(self):
        self.anaconda_setup_complete = getenv("HOME")+'/anaconda3_installComplete.file'

    def install_anaconda(self, misc, shortcuts):
        '''This function downloads anaconda via wget and the link-adress to the linux installer from anaconda.com.
           It then installs anaconda and creates an environment with the required software packages.
           Finally creates the directory three in the sequencing_project folder'''



        if misc.step_completed(self.anaconda_setup_complete, 'Installation of Anaconda3 allready completed.'):
            time.sleep(2.5)
            pass
        else:
            clear_screen()
            print("Download and install Anaconda3\n\n\n")
            print("Downloading and installing anaconda from https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh...")

            # downloads anaconda
            cmd_download = "wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -P $HOME"
            run_command(cmd_download, "")

            # applies read, write and execute permission for all users to the file
            cmd_chmod = "chmod a=xrw $HOME/Anaconda3-2020.11-Linux-x86_64.sh"
            run_command(cmd_chmod, "")

            # installs anaconda
            cmd_install = "bash $HOME/Anaconda3-2020.11-Linux-x86_64.sh"
            run_command(cmd_install, "")


            # create file that shows if anaconda3 is allready installed
            create_trackFile(getenv("HOME")+'/anaconda3_installComplete.file')
            print("\nAnaconda is now successfully installed\n\n")
            input("The terminal must be restarted for anaconda3 to initialize correctly\n\nPress any key to close the terminal\nThen start it again manually")
            misc.close_terminal()

    def create_anaconda_environment(self, misc, shortcuts):

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
        shortcuts.star_output_dir,
        shortcuts.star_index_dir,
        shortcuts.bwa_index_dir
        ])

        if path.isdir(getenv("HOME")+'/anaconda3/envs/sequencing'):
            print("Environment allready exists, do you want to owerwrite?\n")
            if not misc.confirm_choice():
                return

        misc.clear_screen()
        print("Set up a new conda environment for DNA and RNA-sequence analysis.\n")
        print("""


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
    Pyvcf
    Samtools 1.9
    Scipy
    Snpeff
    STAR



    """)

        cmd_env = "conda create -n sequencing bedtools bcftools biopython bwa gatk4 picard python=3.7.6 samtools=1.9 star pandas vcfpy scipy snpeff"
        misc.run_command(cmd_env, "")

        # Delly in bioconda didn't work so I had to do a workaround
        cmd_download_delly = "wget https://github.com/dellytools/delly/releases/download/v0.8.7/delly_v0.8.7_linux_x86_64bit -P $HOME/anaconda3/envs/sequencing/bin"
        misc.run_command(cmd_download_delly, "Download Delly")
        cmd_chmod = "chmod a=xrw $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit"
        misc.run_command(cmd_chmod, "Change permission")
        cmd_mv = "mv $HOME/anaconda3/envs/sequencing/bin/delly_v0.8.7_linux_x86_64bit $HOME/anaconda3/envs/sequencing/bin/delly"
        misc.run_command(cmd_mv, "Move Delly")

        # Manta in bioconda didn't work so I had to do a workaround
        cmd_download_manta = "wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.release_src.tar.bz2 -P $HOME/anaconda3/envs/sequencing/bin"
        misc.run_command(cmd_download_manta, "Download Manta")
        cmd_unpack = "tar -xf $HOME/anaconda3/envs/sequencing/bin/manta-1.6.0.release_src.tar.bz2 -C $HOME/anaconda3/envs/sequencing/bin/"
        misc.run_command(cmd_unpack, "Unpacking Manta")

        # Download annotations for GRCh38 genome
        # cmd_download_annotations = f"java -jar $HOME/anaconda3/envs/sequencing/share/snpeff-5.0-0/snpEff.jar download GRCh38.99"
        # misc.run_command(cmd_download_annotations, "Downloading annotations")


        print("\n\n\nEnvironment named \"sequencing\" was successfully created\n")
        print("You must manually activate the created environment.\nTo do this, type \"conda activate sequencing\" in the shell before restarting the program\n")
        input("Press any key to exit program...")
        sys.exit()
