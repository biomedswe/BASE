from os import listdir, getenv, sys
import subprocess
import multiprocessing
from Bio import SeqIO

class RnaSeqAnalysis():

    def __init__(self):
        pass

    # Shortcuts to files used in RNA sequencing analysis
    self.annotation_gtf_file = getenv("HOME")+"/sequencing_project/reference_genome/gencode.v37.primary_assembly.annotation.gtf"

    # Shortcuts to folders used in RNA sequencing analysis
    self.rna_reads_dir  = getenv("HOME")+"/sequencing_project/rna_seq/reads/"
    self.reference_genome_dir = getenv("HOME")+"/sequencing_project/reference_genome/"
    self.haplotypecaller_output_dir = getenv("HOME")+"/sequencing_project/dna_seq/GATK_haplotypecaller/"

    #---------------------------------------------------------------------------
    def create_directory(path):
        '''This function creates the output directories for the different analysis steps'''
        try:
            # create target directory
            mkdir(path)
        except FileExistsError:
            pass

    #---------------------------------------------------------------------------
    def run_command(command, step):
        '''This function executes a command and checks if it was executes without errors'''

        return_code = subprocess.run(command, shell=True)
        if return_code.returncode == 0:
            print(f"{step} without errors!, continuing with next step...\n")
        else:
            print('\nAn error has occured, see shell for more information. Exiting program...')
            sys.exit()

    #---------------------------------------------------------------------------
    def index_genome(self,choice):
        '''This function indexes either the whole genome or the chromosomes entered'''

        if choice == 1:
            threads = multiprocessing.cpu_count() - 2
            cmd_StarIndex = f'''
        STAR --runThreadN {threads} \\
        --runMode genomeGenerate \\
        --genomeDir {reference_dir}star_index/GRCh38_index \\
        --genomeFastaFiles {reference_file} \\
        --sjdbGTFfile {annotations_gtf}'''

        else:
            threads = multiprocessing.cpu_count() - 2
            print('Index parts of genome\n\n')
            print('''You can add chromosomes separated by a space.
            Use this syntax:

            chr1, chr2, chr3, ..., chr22
            chrX, chrY, chrM
            e.g. "chr11 chr12"
            ''')

            name = self.index_chromosomes()

            cmd_StarIndex = f'''
            STAR --runThreadN {threads} \\
            --runMode genomeGenerate \\
            --genomeDir {reference_dir}star_index/{name}_hg38_index \\
            --genomeFastaFiles {reference_dir}{name}.fa \\
            --sjdbGTFfile {reference_dir}{name}.gtf'''
            print(cmd_StarIndex)
            self.run_command(cmd_StarIndex, 'Index parts of genome completed')

            return name

    #---------------------------------------------------------------------------
    def index_chromosomes(self):
        '''takes one or more chromosome as input, check if syntax is valid and then prints the chromosom id and sequence to a new fasta file'''

        list = ['chr'+str(i) for i in range(1,23)]
        list.extend(['chrX', 'chrY', 'chrM'])

        while True:
            chromosomes = [chr for chr in input("Enter chromosomes to index: ").split()]
            if not chromosomes:
                print('You must enter at least on chromosome')
            if all(chr in list for chr in chromosomes):
                sequences = SeqIO.parse(refGenome_file, 'fasta')
                name = "".join(chromosomes)
                with open(f'reference_genome/{name}.fa', 'w+') as fa, open(f'reference_genome/{name}.bed', 'w') as bed:
                    for chr in chromosomes:
                        for line in sequences:
                            if line.id == chr:
                                fa.write(">" + str(line.id) + "\n")
                                fa.write(str(line.seq)[5:]+ "\n\n")
                                bed.write(str(line.id) + "\t")
                                bed.write("0\t")
                                bed.write(str(len(line.seq)))
                                break
                cmd_createGTF = f"bedtools intersect -a {annotations_gtf_file} -b {refGenome_dir}{name}.bed > {refGenome_dir}{name}.gtf"
                self.run_command(cmd_createGTF, 'New gtf file created')
                return name
            else:
                print('Invalid syntax, please check spelling!')

    #---------------------------------------------------------------------------
    def map_reads(self, options, name):
        '''This function map reads to the reference genome'''

        reads = []
        for read in listdir(rnaReads_dir):
            if options.tumor_id in read:
                reads.append(read)


        threads = multiprocessing.cpu_count() - 2
        print('Map reads to the genome\n\n')

        cmd_mapReads = f'''STAR --genomeDir {refGenome_dir}/star_index/{name}_hg38_index --readFilesIn {rnaReads_dir}{reads[0]} {rnaReads_dir}{reads[1]}
    --runThreadN {threads} --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15 --chimM
    ainSegmentMultNmax 1 --chimOutType Junctions SeparateSAMold WithinBAM SoftClip --chimSegmentMin 15 --genomeLoad NoSharedMemory --limitSjdbInsertNsj 1200000 --outFileNamePrefix star/2064-01 --outFilterIntronMotifs
     None --outFilterMatchNminOverLread 0.33 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterType BySJout --outSAMattributes N
    H HI AS nM NM MD XS ch vA vG vW --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat --waspOutputMode SAMtag --var
    VCFfile {outputFolder_gatk}{options.tumor_id}_filter_RD10.vcf --outSAMattrRGline ID:2064-01_AGTCA_L004 SM:2064-01 LB:2064-01_AGTCA_L004 PL:"ILLUMINA" PU:2064-01_AGTCA_L004 --twopassMode Basic'''
