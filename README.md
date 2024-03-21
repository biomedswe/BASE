# About
Biomedswe Allele-Specific Expression analyser (BASE) v.2.0. 2024


Developer: Minjun Yang, Jonas Andersson                                                                                                                                                                                                
Division of Clinical Genetics                                                                                           
Lund University, BMC C13                                                                                                
SE-221 84 Lund, Sweden                                                                                                                                                                                                                          
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

Required Python Packages: This project depends on several external Python packages. These dependencies are listed in the requirements.txt file included in the project. The following is a brief overview of some key packages:

numpy: A fundamental package for scientific computing with Python.
pandas: An open-source data analysis and manipulation tool.
matplotlib: A comprehensive library for creating static, animated, and interactive visualizations in Python.
scikit-learn: Simple and efficient tools for predictive data analysis.
pysam: A module for reading, manipulating, and writing genomic data sets.
scipy: An open-source software for mathematics, science, and engineering.
To install all required packages, navigate to the root directory of this project and run the following command:

```
pip install -r requirements.txt
```
This command will automatically install all the dependencies listed in the requirements.txt file, ensuring that your project environment is correctly set up and ready to run the BASE.

- **R**: If you don't have R installed, you can download and install it from [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/). We recommend using the latest version of R to ensure compatibility with all packages.

- **Required R Packages**:
  - `dplyr`: A grammar of data manipulation, providing a consistent set of verbs that help you solve the most common data manipulation challenges.
  - `outliers`: A collection of some tests commonly used for identifying outliers.
  - `future.apply`: Implementations of apply function variations that leverage the 'future' framework for asynchronous computing.
  - `facets`: An R package for allele-specific copy number analysis of tumors.
  - `DNAcopy`: A package for analyzing copy number data.
  - `GenomicRanges`: Representation and manipulation of genomic intervals and variables defined along a genome.
  - `Rsamtools`: Provides an interface to the 'samtools', 'bcftools', and 'tabix' utilities for manipulating SAM, BAM, and VCF files.

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

### Initial Setup

Run `reference_genome.py` to download reference sequences, index files, and set up the necessary folder structure.

```bash
python reference_genome.py

This script will:

Download the reference genome and annotation files.
Index the genome using tools like Samtools, BWA, and STAR.
Create a folder structure for the analysis.

```


### WGS Data Processing
To process WGS data, run DNA_alignment.py with the required arguments.
```
python DNA_alignment.py --read1 path/to/read1.fastq --read2 path/to/read2.fastq --output_prefix sample_name [--sample_ploidy 2]
```
--read1 and --read2: Paths to the paired-end FASTQ files.
--output_prefix: A prefix for output files (typically the sample name).
--sample_ploidy (optional): The ploidy of the sample, default is 2. It is recommended to specify this for your project.
This script will:

Align reads to the reference genome.
Perform duplicate marking and indexing.
Run CNV and SNV callers to analyze genomic variations.

### RNA-seq Data Processing for ASE Analysis
After WGS data processing, run RNA_alignment.py to get ASE analysis results.

```
python RNA_alignment.py --read1 path/to/rna_read1.fastq --read2 path/to/rna_read2.fastq --output_prefix sample_name
```
--read1 and --read2: Paths to the RNA-seq paired-end FASTQ files.
--output_prefix: A prefix for output files (typically the sample name).
This script will:

Align RNA-seq reads using STAR.
Run the GATK ASEReadCounter for allele-specific expression counting.
Add WGS information and perform ASE analysis to obtain ASE analysis results.

### Notes
Ensure all scripts and the config.ini file are in the same directory or adjust paths accordingly.
The scripts assume a Unix-like environment with tools like wget, gunzip, and standard Bash commands available.
For detailed instructions on each script's functionality and options, refer to the script's inline comments or documentation sections.


### 4. Cite BASE

Please refer to this repository when using BASE in your project
    
 
