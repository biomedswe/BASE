# sequencing_project
Allele-specific expression (ASE) analysis pipeline

This script automates everything from setting up anaconda to downloading reference genome and performing DNA and RNA analysis to get a result about ASE

## Requirements
- Platform: 
    - linux-64
    - python ≥ 2.7
    
    
## Installation and setup instructions
Open bash (Unix shell) and type the following:

### 1. Clone the git archive:


```
git clone https://github.com/biomedswe/BASE.git $HOME/BASE
```

### 2. Run setup_anaconda3.py and follow instructions in program
Please note that you must have python: ≥ 2.7 installed first.

Type the following in the shell:
```
python2 $HOME/BASE/setup_anaconda3.py
```

Let Anaconda3 install at default location

### 3. Copy your DNA-seq/RNA-seq reads into the right folders

DNA-seq reads: (in fastq.gz format) into $HOME/BASE/dna_seq/reads 

RNA-seq reads: (in fastq.gz format) into $HOME/BASE/rna_seq/reads

or create soft links between these folders and the files with "ln -s"


### 4. Run main.py and follow instructions in program
Type the following in the shell:
```
python3 $HOME/BASE/main.py -t <tumor clinical id> -n <normal clinical id> -sg <sub group> -T <number of threads to use>

e.g. "python3 $HOME/BASE/main.py -t 2064-01 -n 987-02 -sg Heh -T 38"

```

### Optional. Copy Excel document with Copy number information into $HOME/BASE/rna_seq/star/[tumor-id]/[tumor-id]_CN.xlsx

If you have Copy number information, you can copy this document into $HOME/BASE/rna_seq/star/[tumor-id]/[tumor-id]_CN.xlsx 
in order to get calculated pValue and RNA/DNA ratio based on CNV

### 4. Cite BASE

Please refer to this repository when using BASE in your project
    
 
