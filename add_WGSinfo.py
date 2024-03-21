from __future__ import division
import re
from cyvcf2 import VCF
import pysam
import argparse

class WGSInfoProcessor:
    def __init__(self, vcf_path, cnv_path, sample_name):
        self.vcf_path = vcf_path
        self.cnv_path = cnv_path
        self.sample_name = sample_name

    def get_dp(self, format_field, record):
        format_tokens = format_field.split(":")
        record_tokens = record.split(":")
        dp_index = format_tokens.index('DP') if 'DP' in format_tokens else -1
        ad_index = format_tokens.index('AD') if 'AD' in format_tokens else -1

        if ad_index != -1:
            return [int(i) for i in record_tokens[ad_index].split(",")[:2]]
        elif dp_index != -1:
            return int(record_tokens[dp_index]), 0

    def get_gt_dp(self, site):
        vcf_reader = VCF(self.vcf_path)
        vcf_reader.set_samples(self.sample_name)
        chrom, pos = site.split('\t')[:2]
        region = f"{chrom}:{int(pos)-1}-{pos}"

        for record in vcf_reader(region):
            hap_a_dp, hap_b_dp = self.get_dp(str(record).split("\t")[8], str(record).split("\t")[9])
            if record.gt_phases[0]:  # phased result
                return [hap_a_dp, hap_b_dp, int(record.genotypes[0][0] != 0), 1]
            else:  # unphased result
                return sorted([hap_a_dp, hap_b_dp], reverse=True) + [hap_a_dp < hap_b_dp, 0]

        return [0, 0, -1, 0]

    def get_cn(self, site):
        tbx = pysam.TabixFile(self.cnv_path)
        chrom, pos = site.split('\t')[:2]
        records = tbx.fetch(chrom, int(pos)-1, int(pos), parser=pysam.asBed())

        for record in records:
            return [int(record[8]), int(record[9])]
        return [0, -1]

    def adjust_read_depth(self, site, swap_flag):
        if swap_flag:
            tokens = site.split('\t')
            tokens[5], tokens[6] = tokens[6], tokens[5]
            return '\t'.join(tokens)
        return site

    def process_sites(self, input_path, output_path):
        with open(input_path, "r") as ase, open(output_path, "w") as result:
            header = ase.readline().strip()
            result.write(f"#{header}\tWGS_Hap1_DP\tWGS_Hap2_DP\ttotal_cn\thom1_cn\tRS_flag\tPhase_stat\n")

            for line in ase:
                line = line.strip()
                hap1, hap2, rs_flag, phase_stat = self.get_gt_dp(line)
                t_cn, h1_cn = self.get_cn(line) if self.cnv_path else [0, -1]
                h2_cn = t_cn - h1_cn

                if h1_cn == h2_cn and phase_stat == 0:
                    rs_flag = 2

                if t_cn != -1 and rs_flag != -1:
                    rna_info = self.adjust_read_depth(line, rs_flag)
                    result.write(f"{rna_info}\t{hap1}\t{hap2}\t{t_cn}\t{h1_cn}\t{int(rs_flag)}\t{phase_stat}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AI_CNV')
    parser.add_argument('-i', '--input', required=True, help='GATK readscount')
    parser.add_argument('-v', '--vcf', required=True, help='WXS_coding.het.vcf.gz')
    parser.add_argument('-c', '--cnv', required=True, help='segment.cnv.gz')
    parser.add_argument('-s', '--sample', required=True, help='sample_name')
    parser.add_argument('-o', '--output', required=True, help='segment.imba.bed')
    options = parser.parse_args()

    processor = WGSInfoProcessor(options.vcf, options.cnv, options.sample)
    processor.process_sites(options.input, options.output)
