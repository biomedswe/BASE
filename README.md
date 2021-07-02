# sequencing_project
Allele-specific expression (ase) analysis pipeline

This script automates everything from setting up anaconda to downloading reference genome and performing DNA and RNA analysis to get a result about ase

## Requirements
- Platform: 
    - linux-64
    - python >= 3.6
    
    
## Installation and setup instructions
Open bash (Unix shell) and type the following:

### 1. Clone the git archive:


```
git clone https://github.com/biomedswe/sequencing_project.git $HOME/sequencing_project
```

### 2. Run setup_anaconda3.py and follow instructions in program
Please note that you must have python: ≥ 2.7 installed first.

Type the following in the shell:
```
python2 setup_anaconda3.py (after navigating to $HOME/sequencing_project)
```
or
```
python2 $HOME/sequencing_project/setup_anaconda3.py
```

Let Anaconda3 install at default location


### 3. Run main.py and follow instructions in program
Type the following in the shell:
```
python3 main.py -t <tumor clinical id> -n <normal clinical id> (after navigating to $HOME/sequencing_project)
```
or
```
python3 $HOME/sequencing_project/main.py -t <tumor clinical id> -n <normal clinical id>
```

### 4. Copy your DNA-seq/RNA-seq reads into the right folders

DNA-seq reads: (in fasta/fastq format) into $HOME/sequencing_project/dna_seq/reads 

RNA-seq reads: (in fasta.gz) into $HOME/sequencing_project/rna_seq/reads

### 4. Cite BASE

Please refer to this repository when using BASE in your project
    
 
