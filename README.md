# sequencing_project



Requirements
Platform:
linux-64
Software:
python >= 3.6
conda >= 4.8.2
Installation and setup instructions
From source code
Clone the git archive:

git clone https://github.com/bioinformatics-polito/PhyliCS.git

This package will simplify your DNA and RNA-sequening.

Follow stese steps:

1. Set-up this git repository in your $HOME folder (eg. Home/biomedswe/) 
	you can for example type "git clone <link to git repository>" in your $HOME folder


2. Navigate to the sequencing_project folder and run setup_anaconda.py to install all neccessary software packages

	a) Install anaconda and the restart your shell. (If you allready have installed anaconda, proceed with step c)
	b) Navigate back to the sequencing_project folder
	c) Choose: Set up a new conda environment for DNA and RNA-sequence analysis
	d) Exit the script and type "conda activate sequencing"

3. In the sequencing_project folder, run download_reference_genome.py to download and save GRCh38.p13.genome.fa.
	
	It will be saved in sequencing_project/reference_genome/GRCh38.p13.genome.fa



4. Copy your folder with WGS reads into seqquencing_project/dna_seq/reads/


5. In the sequencing_project folder, run dna_seq.py to start DNA analysis






3. Copy your folder with RNA-sequencing reads into the RNA folder.
   Run RNA/rna_se.py to start RNA analysis
