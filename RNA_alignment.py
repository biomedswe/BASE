import os
import argparse
from shortcuts import Shortcuts
from Misc2 import Misc
import shutil
from add_WGSinfo import WGSInfoProcessor
from ASE_analysis import ASEAnalysis


class RNA_alignment:
    def __init__(self):
        self.shortcuts = Shortcuts()
        self.misc = Misc(dry_run=False)
        
    def align_reads_with_star(self, read1, read2, output_prefix):
        try:
            output_dir = os.path.join(self.shortcuts.rna_seq_dir)
            os.makedirs(output_dir, exist_ok=True)
            
            vcf_file = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.filtered.hc.vcf")

            # Construct the STAR command with specified options
            star_cmd = [
                self.shortcuts.STAR_path,
                "--genomeDir", self.shortcuts.star_index_dir,
                "--runThreadN", "32",
                "--alignIntronMax", "1000000",
                "--alignIntronMin", "20",
                "--alignMatesGapMax", "1000000",
                "--alignSJDBoverhangMin", "1",
                "--alignSJoverhangMin", "8",
                "--chimJunctionOverhangMin", "15",
                "--chimMainSegmentMultNmax", "1",
                "--chimOutType", "Junctions SeparateSAMold WithinBAM SoftClip",
                "--chimSegmentMin", "15",
                "--limitSjdbInsertNsj", "1200000",                
                "--outSAMtype", "BAM", "Unsorted",
                "--readFilesCommand", "zcat",
                "--outFilterMatchNminOverLread", "0.33",
                "--outFilterMismatchNmax", "999",
                "--outFilterMismatchNoverLmax", "0.1",
                "--outFilterMultimapNmax", "20",
                "--outFilterScoreMinOverLread", "0.33",
                "--outFilterType", "BySJout",
                "--outSAMattributes", "NH HI AS nM NM MD XS ch vA vG vW",
                "--outSAMstrandField", "intronMotif",
                "--outSAMunmapped", "Within",
                "--quantMode", "TranscriptomeSAM GeneCounts",
                "--waspOutputMode", "SAMtag",
                "--readFilesIn", read1, read2,
                "--varVCFfile", vcf_file,
                "--outFileNamePrefix", os.path.join(output_dir, output_prefix + "."),
                "--outSAMattrRGline", f"ID:{output_prefix} SM:{output_prefix} LB:{output_prefix} PL:ILLUMINA PU:{output_prefix}",
                "--twopassMode", "Basic"
            ]

            self.misc.log_to_file("INFO", "Aligning RNA-seq reads with STAR...")
            self.misc.run_command(" ".join(star_cmd))

            # After STAR alignment, sort the BAM file with samtools
            unsorted_bam = os.path.join(output_dir, f"{output_prefix}.Aligned.out.bam")
            sorted_bam = os.path.join(output_dir, f"{output_prefix}.STAR_WASP.bam")
            sort_cmd = [
                self.shortcuts.samtools_path, "sort", "-@ 4", "-m 5G", "-o", sorted_bam, unsorted_bam
            ]
            self.misc.run_command(" ".join(sort_cmd))
            
            # Filter the sorted BAM file for WASP-passing reads
            pass_bam = os.path.join(output_dir, f"{output_prefix}.STAR_WASP.pass.bam")
            filter_cmd = f"samtools view -h {sorted_bam} | LC_ALL=C egrep '^@|vW:i:1' | {self.shortcuts.samtools_path} view -bo {pass_bam} -"
            self.misc.run_command(filter_cmd)

            index_cmd = [self.shortcuts.samtools_path, "index", pass_bam]
            self.misc.run_command(" ".join(index_cmd))

            self.misc.log_to_file("INFO", "RNA-seq alignment with STAR and post-processing with samtools completed successfully.")
            
            self.misc.log_to_file("INFO", "Cleaning up STAR output files...")
            files_to_remove = [
                f"{output_prefix}._STARgenome",
                f"{output_prefix}._STARpass1",
                f"{output_prefix}.Aligned.out.bam",
                f"{output_prefix}.Chimeric.out.sam",
                f"{output_prefix}.Aligned.toTranscriptome.out.bam",
                f"{output_prefix}.SJ.out.tab",
                f"{output_prefix}.Chimeric.out.junction"
            ]
            for file in files_to_remove:
                file_path = os.path.join(output_dir, file)
                if os.path.exists(file_path):
                    if os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                    else:
                        os.remove(file_path)

            self.misc.log_to_file("INFO", "RNA-seq alignment with STAR and cleanup completed successfully.")

        except Exception as e:
            self.misc.log_exception(self, e, "Error in align_reads_with_star")


    def run_ase_read_counter(self, output_prefix):
        try:
            self.misc.log_to_file("INFO", "Running GATK ASEReadCounter...")

            input_bam = os.path.join(self.shortcuts.rna_seq_dir, f"{output_prefix}.STAR_WASP.pass.bam")
            filtered_vcf = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.filtered.hc.vcf")
            output_csv = os.path.join(self.shortcuts.rna_seq_dir, f"{output_prefix}.STAR_WASP.csv")

            # Constructing the GATK ASEReadCounter command
            ase_read_counter_cmd = [
                self.shortcuts.gatk_path,
                "ASEReadCounter",
                "--reference", self.shortcuts.reference_genome_file,
                "--min-mapping-quality", "10",
                "--min-depth-of-non-filtered-base", "10",
                "--min-base-quality", "2",
                "--disable-read-filter", "NotDuplicateReadFilter",
                "--disable-read-filter", "NotSecondaryAlignmentReadFilter",
                "--variant", filtered_vcf,
                "-I", input_bam,
                "--output", output_csv
            ]
            self.misc.run_command(" ".join(ase_read_counter_cmd))
            self.misc.log_to_file("INFO", "GATK ASEReadCounter completed successfully.")

        except Exception as e:
            self.misc.log_exception(self, e, "Error in run_ase_read_counter")
            

    def add_wgs_info(self, output_prefix):
        try:
            self.misc.log_to_file("INFO", "Adding WGS information...")

            input_path = os.path.join(self.shortcuts.rna_seq_dir, f"{output_prefix}.STAR_WASP.csv")
            output_path = os.path.join(self.shortcuts.rna_seq_dir, f"{output_prefix}.ASEinput.tab")
            
            vcf_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.filtered.hc.vcf.gz")            
            cnv_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}_cnv_prediction.bed.gz")
            sample_name = output_prefix

            addwgs_processor = WGSInfoProcessor(vcf_path, cnv_path, sample_name)
            addwgs_processor.process_sites(input_path, output_path)

            self.misc.log_to_file("INFO", "WGS information added successfully.")

        except Exception as e:
            self.misc.log_exception(self, e, "Error in add_wgs_info")


    def perform_ase_analysis(self, output_prefix):
        try:
            self.misc.log_to_file("INFO", "Starting ASE analysis...")

            self.add_wgs_info(output_prefix)

            ase_input_file = os.path.join(self.shortcuts.rna_seq_dir, f"{output_prefix}.ASEinput.tab")
            bed_file = self.shortcuts.annotation_exon_file + ".gz"
            ase_output_file = os.path.join(self.shortcuts.rna_seq_dir, f"{output_prefix}.ASEAnalysisResults.tab")

            ase_analyzer = ASEAnalysis(ase_input_file, bed_file, ase_output_file)
            ase_analyzer.run_analysis()

            self.misc.log_to_file("INFO", "ASE analysis completed successfully.")
        
        except Exception as e:
            self.misc.log_exception(self, e, "Error in perform_ase_analysis")
 
            
            
def main():
    parser = argparse.ArgumentParser(description='Process RNAseq input files and parameters.')
    parser.add_argument('--read1', required=True, help='Path to the first read file')
    parser.add_argument('--read2', required=True, help='Path to the second read file')
    parser.add_argument('--output_prefix', required=True, help='Prefix for the output files')
    args = parser.parse_args()

    print(f"Read1: {args.read1}")
    print(f"Read2: {args.read2}")
    print(f"Output Prefix: {args.output_prefix}")
    
    rna_alignment = RNA_alignment()
    
    # Align reads with STAR
    rna_alignment.align_reads_with_star(args.read1, args.read2, args.output_prefix)

    # Run GATK ASEReadCounter
    rna_alignment.run_ase_read_counter(args.output_prefix)

    # Perform ASE analysis including adding WGS information
    rna_alignment.perform_ase_analysis(args.output_prefix)

if __name__ == "__main__":
    main()