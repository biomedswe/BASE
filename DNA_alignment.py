import os
import argparse
from shortcuts import Shortcuts
from Misc2 import Misc
from gvcf2reads import gvcf2cnv
from cnv_calling import cnv_calling
from ImbaCalculator import ImbaCalculator
from adj_logr import LogRAdjuster
from cnv_train import CNVTrainerPredictor


class DNA_alignment:
    def __init__(self):
        self.shortcuts = Shortcuts()
        self.misc = Misc(dry_run=False)
        
    def align_and_process_reads(self, read1, read2, output_prefix):
        try:
            bam_file = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.bam")
            read_group = f"@RG\\tID:{output_prefix}\\tSM:{output_prefix}\\tLB:lib1\\tPL:ILLUMINA\\tPU:{output_prefix}"
            self.misc.log_to_file("INFO", "Aligning reads with bwa mem and converting to BAM...")
            align_and_convert_cmd = f"{self.shortcuts.bwa_path} mem -M -R \"{read_group}\" -t 8 {self.shortcuts.reference_genome_file} {read1} {read2} | {self.shortcuts.samtools_path} view -bS - > {bam_file}"
            self.misc.run_command(align_and_convert_cmd)

            sorted_bam = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}_sorted.bam")
            self.misc.log_to_file("INFO", "Sorting alignments with sambamba sort...")
            sort_cmd = f"{self.shortcuts.sambamba_path} sort -o {sorted_bam} -t 8 {bam_file}"
            self.misc.run_command(sort_cmd)

            dedup_bam = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}_dedup.bam")
            self.misc.log_to_file("INFO", "Removing duplicates with sambamba markdup...")
            dedup_cmd = f"{self.shortcuts.sambamba_path} markdup -t 8 --tmpdir={self.shortcuts.tmp} {sorted_bam} {dedup_bam}"
            self.misc.run_command(dedup_cmd)

            self.misc.log_to_file("INFO", "Indexing deduplicated BAM with sambamba index...")
            index_cmd = f"{self.shortcuts.sambamba_path} index -t 8 {dedup_bam}"
            self.misc.run_command(index_cmd)

            if os.path.exists(dedup_bam + ".bai"):
                self.misc.log_to_file("INFO", "Cleaning up intermediate files...")
                os.remove(bam_file)
                os.remove(sorted_bam)

            self.misc.log_to_file("INFO", "Alignment, sorting, deduplication, and indexing completed successfully.")
        except Exception as e:
            self.misc.log_exception("align_and_process_reads", e, "Error in align_and_process_reads")



    def run_cnv_haplotypecaller(self, bam_file, output_prefix):
        try:
            self.misc.log_to_file("INFO", "Running GATK HaplotypeCaller for CNV analysis...")

            cnv_intervals_file = self.shortcuts.reference_genome_cnv_region
            cnv_output_vcf = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.cnv.gvcf")

            cnv_haplotypecaller_cmd = [
                self.shortcuts.gatk_path,
                "HaplotypeCaller",
                "-R", self.shortcuts.reference_genome_file,
                "-I", bam_file,
                "-O", cnv_output_vcf,
                "--emit-ref-confidence GVCF",
                "-L", cnv_intervals_file
            ]
            self.misc.log_to_file("INFO", f"Executing: {cnv_haplotypecaller_cmd}")
            self.misc.run_command(" ".join(cnv_haplotypecaller_cmd))
            self.misc.log_to_file("INFO", "GATK HaplotypeCaller execution for CNV analysis completed successfully.")
            
            self.misc.log_to_file("INFO", "Compressing VCF with bgzip...")
            bgzip_cmd = f"{self.shortcuts.bgzip_path} -c {cnv_output_vcf} > {cnv_output_vcf}.gz"
            self.misc.run_command(bgzip_cmd)

            self.misc.log_to_file("INFO", "Indexing VCF with tabix...")
            tabix_cmd = f"{self.shortcuts.tabix_path} -p vcf {cnv_output_vcf}.gz"
            self.misc.run_command(tabix_cmd)

            self.misc.log_to_file("INFO", "Filtering VCF with bcftools and excluding lines with './.'...")
            filtered_vcf_gz = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.cnv.use.gvcf.gz")
            bcftools_cmd = f"{self.shortcuts.bcftools_path} filter -g 20 -e 'FORMAT/DP < 10 || FORMAT/GQ < 35' {cnv_output_vcf}.gz | grep -v '\\./\\.' | {self.shortcuts.bgzip_path} -c > {filtered_vcf_gz}"
            self.misc.run_command(bcftools_cmd)
            
            self.misc.log_to_file("INFO", f"Indexing filtered VCF {filtered_vcf_gz} with tabix...")
            tabix_index_cmd = f"{self.shortcuts.tabix_path} -p vcf {filtered_vcf_gz}"
            self.misc.run_command(tabix_index_cmd)
            self.misc.log_to_file("INFO", "Tabix indexing of filtered VCF completed successfully.")
            
        except Exception as e:
            self.misc.log_exception(self, e, "Error in run_cnv_haplotypecaller")


    def run_cnv_analysis(self, output_prefix):
        
        input_vcf_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.cnv.use.gvcf.gz")
        gvcf2cnv_output_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}")

        cnv_processor = gvcf2cnv(input_vcf_path, output_prefix, gvcf2cnv_output_path)
        cnv_processor.run()
 
        cnv_caller = cnv_calling()     
        cnv_output_dir = self.shortcuts.aligned_output_dir
        cnv_caller.run_r_script(gvcf2cnv_output_path, cnv_output_dir)

        probes_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.probes")
        snps_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.snps")
    
        def remove_first_line(file_path):
            with open(file_path, 'r') as f:
                lines = f.readlines()[1:]
            with open(file_path, 'w') as f:
                f.writelines(lines)
    
        remove_first_line(probes_path)
        remove_first_line(snps_path)
    

        def compress_and_index(file_path):
            compressed_path = file_path + '.gz'
            self.misc.run_command(f"{self.shortcuts.bgzip_path} -c {file_path} > {compressed_path}")
            self.misc.run_command(f"{self.shortcuts.tabix_path} -p vcf {compressed_path}")
    
        compress_and_index(probes_path)
        compress_and_index(snps_path)


    def correct_cnv_ratio(self, output_prefix, ploidy=2):
        snps_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.snps.gz")
        seg_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.sd1.sdundo.seg")        
        imba_output_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.imba")

        # Run the IMBA calculation
        imba_calculator = ImbaCalculator(baf_path=snps_path, interval_path=seg_path, output_path=imba_output_path)
        imba_calculator.run()
        self.misc.log_to_file("INFO", "IMBA calculation completed successfully.")
        
        # Run the LogRAdjuster
        probes_path_gz = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.probes.gz")
        adj_logr_output_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.adj_logr")        
        logr_adjuster = LogRAdjuster(logr_path=probes_path_gz, segment_path=imba_output_path, output_path=adj_logr_output_path, ploidy=ploidy)
        logr_adjuster.adjust_logR()
        self.misc.log_to_file("INFO", "LogR adjustment completed successfully.")       

        # Call the cnv_predict method of CNVTrainerPredictor
        cnv_prediction_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}_cnv_prediction.bed")
        cnv_trainer = CNVTrainerPredictor()
        cnv_trainer.cnv_predict(segment_file=adj_logr_output_path, output_file=cnv_prediction_path)
        
        cnv_prediction_compressed_path = f"{cnv_prediction_path}.gz"
        bgzip_command = f"{self.shortcuts.bgzip_path} -c {cnv_prediction_path} > {cnv_prediction_compressed_path}"
        self.misc.run_command(bgzip_command)
        tabix_command = f"{self.shortcuts.tabix_path} -p bed {cnv_prediction_compressed_path}"
        self.misc.run_command(tabix_command)

        os.remove(cnv_prediction_path)


    def run_snv_haplotypecaller_for_target(self, dedup_bam, output_prefix):
        try:
            self.misc.log_to_file("INFO", "Running GATK HaplotypeCaller for target...")
            intervals_file = os.path.join(self.shortcuts.reference_genome_dir, "parsed_gtf_merged.bed")
            output_vcf = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.hc.vcf")

            # Run GATK HaplotypeCaller
            haplotypecaller_cmd = [
                self.shortcuts.gatk_path,
                "HaplotypeCaller",
                "-R", self.shortcuts.reference_genome_file,
                "-I", dedup_bam,
                "-O", output_vcf,
                "-L", intervals_file
            ]
            self.misc.run_command(" ".join(haplotypecaller_cmd))
            
            hc_running_vcf = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.hc.running.vcf")
            filter_cmd = [
                self.shortcuts.gatk_path,
                "VariantFiltration",
                "-R", self.shortcuts.reference_genome_file,
                "-V", output_vcf,
                "-O", hc_running_vcf,
                "--filter-expression", "\"DP > 15 && QD > 2.0 && MQ > 35.0 && MQRankSum > -12.5 && ReadPosRankSum > -8.0 && FS < 60.0\"",
                "--filter-name", "\"SNV_FILTER\""
            ]
            self.misc.run_command(" ".join(filter_cmd))
            self.misc.log_to_file("INFO", "GATK HaplotypeCaller execution and filtering completed successfully.")


            filtered_vcf = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.filtered.hc.vcf")
            selectvariants_cmd = [
                self.shortcuts.gatk_path,
                "SelectVariants",
                "-R", self.shortcuts.reference_genome_file,
                "-V", hc_running_vcf,
                "-O", filtered_vcf,
                "--select-type-to-include", "SNP"
            ]
            self.misc.run_command(" ".join(selectvariants_cmd))
            self.misc.log_to_file("INFO", "Selecting only SNPs from the VCF completed successfully.")


            index_cmd = [
                self.shortcuts.gatk_path,
                "IndexFeatureFile",
                "--input", filtered_vcf
            ]
            self.misc.run_command(" ".join(index_cmd))
            self.misc.log_to_file("INFO", "GATK HaplotypeCaller execution and VCF indexing completed successfully.")

            compressed_vcf_path = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.filtered.hc.vcf.gz")
            bgzip_command = f"{self.shortcuts.bgzip_path} -c {filtered_vcf} > {compressed_vcf_path}"
            self.misc.run_command(bgzip_command)
            tabix_command = f"{self.shortcuts.tabix_path} -p vcf {compressed_vcf_path}"
            self.misc.run_command(tabix_command)            

            filtered_het_vcf = os.path.join(self.shortcuts.aligned_output_dir, f"{output_prefix}.filtered_het.hc.vcf")
            bcftools_command = f"{self.shortcuts.bcftools_path} view -i 'GT=\"0/1\"' -g het {compressed_vcf_path} > {filtered_het_vcf}"
            self.misc.run_command(bcftools_command)
            
            index_het_vcf_command = f"{self.shortcuts.gatk_path} IndexFeatureFile -I {filtered_het_vcf}"
            self.misc.run_command(index_het_vcf_command)
            self.misc.log_to_file("INFO", "Indexing of heterozygous filtered VCF completed successfully.")
            

        except Exception as e:
            self.misc.log_exception(self, e, "Error in run_snv_haplotypecaller_for_target")
            

def main():
    parser = argparse.ArgumentParser(description='Process input files and parameters.')
    parser.add_argument('--read1', required=True, help='Path to the first read file')
    parser.add_argument('--read2', required=True, help='Path to the second read file')
    parser.add_argument('--sample_ploidy', type=int, default=2, help='Sample ploidy (default: 2)')
    parser.add_argument('--output_prefix', required=True, help='Prefix for the output files')
    args = parser.parse_args()

    print(f"Read1: {args.read1}")
    print(f"Read2: {args.read2}")
    print(f"Sample Ploidy: {args.sample_ploidy}")
    print(f"Output Prefix: {args.output_prefix}")
    
    dna_alignment = DNA_alignment()
    dedup_bam = os.path.join(dna_alignment.shortcuts.aligned_output_dir, f"{args.output_prefix}_dedup.bam")

    dna_alignment.align_and_process_reads(args.read1, args.read2, args.output_prefix)

    # CNV calling 
    dna_alignment.run_cnv_haplotypecaller(dedup_bam, args.output_prefix)
    dna_alignment.run_cnv_analysis(args.output_prefix)
    dna_alignment.correct_cnv_ratio(args.output_prefix, ploidy=args.sample_ploidy)
    
    # SNV calling 
    dna_alignment.run_snv_haplotypecaller_for_target(dedup_bam, args.output_prefix)


if __name__ == "__main__":
    main()
