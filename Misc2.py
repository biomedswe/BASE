import os
import sys
import subprocess

class Misc:
    def __init__(self, dry_run=False):
        self.dry_run = dry_run

    def log_to_file(self, level, message):
        print(f"{level}: {message}")

    def log_exception(self, context, exception, message):
        print(f"Exception in {context}: {message} - {str(exception)}")

    def run_command(self, command):
        if self.dry_run:
            print(f"[DRY-RUN] Would execute: {command}")
            self.log_to_file("INFO", f"[DRY-RUN] Would execute: {command}")
        else:
            try:
                completed_process = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"Executed successfully: {command}")
                self.log_to_file("INFO", f"Executed successfully: {command}")
            except subprocess.CalledProcessError as e:
                error_message = e.stderr.decode()
                self.log_exception("run_command", error_message, f"Command failed: {e.cmd}")
                sys.exit(1)

# Example of using the Misc class with dry_run
# misc = Misc(dry_run=True)  # Set dry_run=True for a dry run, False for actual execution
# misc.run_command("ls -l")  # This will log the command without executing it if dry_run is True
