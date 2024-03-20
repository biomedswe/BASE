# About
Biomedswe Allele-Specific Expression analyser (BASE) v.1.0. 2021


Developer: Minjun Yang, Jonas Andersson                                                                                                         
Master's programme in biomedicine                                                                                       
Division of Clinical Genetics                                                                                           
Lund University, BMC C13                                                                                                
SE-221 84 Lund, Sweden                                                                                                                                                                                                                          
# BASE

This toolkit provides a set of scripts to process WGS and RNA-seq data for genomic and transcriptomic analyses. The workflow includes downloading reference sequences, indexing, folder structure setup, DNA alignment for WGS data, and RNA alignment for ASE analysis.

## Setup

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


# WGS Data Processing
To process WGS data, run DNA_alignment.py with the required arguments.

python DNA_alignment.py --read1 path/to/read1.fastq --read2 path/to/read2.fastq --output_prefix sample_name [--sample_ploidy 2]

--read1 and --read2: Paths to the paired-end FASTQ files.
--output_prefix: A prefix for output files (typically the sample name).
--sample_ploidy (optional): The ploidy of the sample, default is 2. It is recommended to specify this for your project.
This script will:

Align reads to the reference genome.
Perform duplicate marking and indexing.
Run CNV and SNV callers to analyze genomic variations.

# RNA-seq Data Processing for ASE Analysis
After WGS data processing, run RNA_alignment.py to get ASE analysis results.

python RNA_alignment.py --read1 path/to/rna_read1.fastq --read2 path/to/rna_read2.fastq --output_prefix sample_name

--read1 and --read2: Paths to the RNA-seq paired-end FASTQ files.
--output_prefix: A prefix for output files (typically the sample name).
This script will:

Align RNA-seq reads using STAR.
Run the GATK ASEReadCounter for allele-specific expression counting.
Add WGS information and perform ASE analysis to obtain ASE analysis results.

# Notes
Ensure all scripts and the config.ini file are in the same directory or adjust paths accordingly.
The scripts assume a Unix-like environment with tools like wget, gunzip, and standard Bash commands available.
For detailed instructions on each script's functionality and options, refer to the script's inline comments or documentation sections.


### 4. Cite BASE

Please refer to this repository when using BASE in your project
    
 
