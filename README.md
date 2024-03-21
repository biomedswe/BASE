                                                                                                                                                                                                            
# BASE

This toolkit provides a set of scripts to process WGS and RNA-seq data for genomic and transcriptomic analyses. The workflow includes downloading reference sequences, indexing, folder structure setup, DNA alignment for WGS data, and RNA alignment for ASE analysis.

## Setup

### Prerequisites

Before you can run BASE, you need to ensure that your environment is set up correctly. The following are required:

- **Python**: The core of this project is written in Python, so you'll need Python installed on your computer. We recommend using Python 3.7 or newer to ensure compatibility with all dependencies.

- **Virtual Environment (Optional, but recommended)**: Using a virtual environment for Python projects helps manage dependencies and avoid conflicts with other projects. You can create a virtual environment using Python's built-in `venv` module or with `virtualenv` for older versions of Python.

  To create a virtual environment using `venv`, run the following command in your terminal (replace `myenv` with your preferred environment name):
  ```bash
  python3 -m venv myenv
  ```
  To activate the virtual environment, on Windows, run:
  ```
  myenv\Scripts\activate.bat
  ```
  On Unix or MacOS, run:
  ```
  source myenv/bin/activate
  ```
  
- **Required Python Packages**:
  - `numpy`: version 1.23.4
  - `scipy`: version 1.9.1
  - `pandas`: version 1.5.3
  - `matplotlib`: version 3.3.2
  - `joblib`: version 1.2.0
  - `scikit-learn`: version 1.2.2
  - `pysam`: version 0.15.2
  - `cyvcf2`: version 0.30.18
  
   ```
   pip install -r requirements.txt
   ```
   This command will automatically install all the dependencies listed in the requirements.txt file, ensuring that your project environment is correctly set up and ready to 
   run the BASE.


- **R**: If you don't have R installed, you can download and install it from [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/). We recommend using the latest version of R to ensure compatibility with all packages.   

- **Required R Packages**:
  - `dplyr`: A grammar of data manipulation, providing a consistent set of verbs that help you solve the most common data manipulation challenges.
  - `outliers`: A collection of some tests commonly used for identifying outliers.
  - `future.apply`: Implementations of apply function variations that leverage the 'future' framework for asynchronous computing.
  - `facets`: An R package for allele-specific copy number analysis of tumors.
  - `DNAcopy`: A package for analyzing copy number data.
  - `GenomicRanges`: Representation and manipulation of genomic intervals and variables defined along a genome.
  - `Rsamtools`: Provides an interface to the 'samtools', 'bcftools', and 'tabix' utilities for manipulating SAM, BAM, and VCF files.
  
  To install all required packages, You will be using the R console for entering the commands provided below.
  Install CRAN Packages
  ```
  install.packages("dplyr")
  install.packages("outliers")
  install.packages("future.apply")
  ```
  To install packages from Bioconductor, you first need to ensure that you have the BiocManager package installed. If not, you can install it using the following command:
  Install Bioconductor Packages
  ```
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  ```
  Once BiocManager is installed, you can proceed to install the Bioconductor packages with the following commands:
  ```
  BiocManager::install("facets")
  BiocManager::install("DNAcopy")
  BiocManager::install("GenomicRanges")
  BiocManager::install("Rsamtools")
  ```

Please ensure all the required software and packages are installed before proceeding with the project setup and execution.


### Configuring Third-Party Programs

1. Before running the scripts, ensure all third-party programs (e.g., GATK, BWA, Samtools, etc.) are installed on your system.
2. Set the paths to these programs in the `config.ini` file. An example configuration is provided below:

    ```ini
    [third_party_programs]
    gatk_path = /path/to/gatk
    bwa_path = /path/to/bwa
    samtools_path = /path/to/samtools
    STAR_path = /path/to/STAR
    bedtools_path = /path/to/bedtools
    ...
    ```

    Replace `/path/to/...` with the actual paths to the installed third-party programs.


 detailed instructions on each script's functionality and options, refer to the script's inline comments or documentation sections.

### Initial Setup

  Run `reference_genome.py` to download reference sequences, index files, and set up the necessary folder structure.   
  ```bash
  python reference_genome.py
  ```
  This script will:   
  Download the reference genome and annotation files.   
  Index the genome using tools like Samtools, BWA, and STAR.   
  Create a folder structure for the analysis.   

Before running the analysis, ensure the reference files and annotation files required by the analysis are properly set up and accessible.

### Notes
The scripts assume a Unix-like environment with tools like wget, gunzip, and standard Bash commands available.
Ensure all scripts and the config.ini file are in the same directory or adjust paths accordingly.



## Usage 
### WGS Data Processing
To process WGS data, run DNA_alignment.py with the required arguments.
  ```
  python DNA_alignment.py --read1 path/to/read1.fastq --read2 path/to/read2.fastq --output_prefix sample_name [--sample_ploidy 2]
  ```
  
  --read1 and --read2: Paths to the paired-end FASTQ files.   
  --output_prefix: A prefix for output files (typically the sample name).   
  --sample_ploidy (optional): The ploidy of the sample, default is 2. It is recommended to specify this for your project.   
  
  This script will:   
  Run WGS reads alignment, copy number variations calling and exonic SNVs calling.

### RNA-seq Data Processing for ASE Analysis
After WGS data processing, run RNA_alignment.py to get ASE analysis results.

```
python RNA_alignment.py --read1 path/to/rna_read1.fastq --read2 path/to/rna_read2.fastq --output_prefix sample_name
```
--read1 and --read2: Paths to the RNA-seq paired-end FASTQ files.
--output_prefix: A prefix for output files (typically the sample name).
This script will:

Align RNA-seq reads using STAR.
Run RNA reads alignment, expressed heterozygous SNVs coverage counting and allele-specific expression analysis.

## Result
This comprehensive file delineates the results of the ASE (Allele-Specific Expression) analysis for individual genes, encapsulating aspects such as 
gene symbols, 
GENCODE identifiers, 
Chromosome numebr of the gene under examination, 
genomic coordinates of the gene under examination,
Haplotype reads count derived from Whole Genome Sequencing (WGS) data,
Haplotype reads count derived from RNA sequencing (RNAseq) data, 
Total and allele-specific copy numbers of the gene under examination, 
The number of informative (expressed) heterozygous SNVs of the gene under examination, 
P-values: those calculated from ASE analysis leveraging WGS/RNAseq read counts, those predicated on the assumption of equal expression across both alleles (model 1), those derived from the genomic allele ratio (model 2),
ASE odds ratios informed by models 1 and 2.

Review Results: Open and review the generated .tab file using a text editor or data analysis tool of your choice (e.g., Excel, R, Python pandas DataFrame). Each row corresponds to a gene analyzed, with columns providing the statistical analysis results and relevant genomic information.


# About
Biomedswe Allele-Specific Expression analyser (BASE) v.2.0. 2024

Developer: Minjun Yang, Jonas Andersson, Efe Aydın                                                                                                                                                                                              
Division of Clinical Genetics                                                                                           
Lund University, BMC C13                                                                                                
SE-221 84 Lund, Sweden 

### 4. Cite BASE

Please refer to this repository when using BASE in your project
    
 
