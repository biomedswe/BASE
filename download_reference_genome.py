import subprocess

print("Downloading GRCh38.p13.genome.fa...")

cdm_download = "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz -P $HOME"
