from __future__ import division
import re
import random
import numpy as np
from cyvcf2 import VCF
from shortcuts import Shortcuts


class gvcf2cnv:
    def __init__(self, vcf_path, sample_name, output_path):
        self.shortcuts = Shortcuts()
        self.vcf_path = vcf_path
        self.sample_name = sample_name
        self.reference_panel_path = self.shortcuts.reference_cnv_panel
        self.output_path = output_path
        self.BASE_COVERAGE = 40

    def read_panel_npz(self):
        loaded_data = np.load(self.reference_panel_path, allow_pickle=True)
        if 'data' in loaded_data:
            loaded_structured_array = loaded_data['data']
            data_dict = {}
            for row in loaded_structured_array:
                chromosome = row['Chromosome']
                position = row['Position']
                value = row['Value']
                data_dict.setdefault(chromosome, {})[position] = value
            return data_dict
        else:
            print("Key 'data' not found in the .npz file. Available keys:", loaded_data.files)
            return {}

    def getDP(self, info1, format1, record):
        array_format = format1.split(":")
        array_record = record.split(":")
        array_info = info1.split(";")
        
        raw_mq_and_dp_field = [field for field in array_info if "RAW_MQandDP" in field]
        second_number = 30
        if raw_mq_and_dp_field:
            second_number = int(re.search(r'\d+', raw_mq_and_dp_field[0].split(",")[1]).group())
        
        GT_idx = next((i for i, item in enumerate(array_format) if item == 'GT'), None)
        DP_idx = next((i for i, item in enumerate(array_format) if item == 'DP'), None)
        AD_idx = next((i for i, item in enumerate(array_format) if item == 'AD'), None)
        
        if GT_idx is not None:
            gts = array_record[GT_idx].replace("|", "/").split("/")
            gt_value = [int(gt) for gt in gts if gt.isdigit()]
            if len(gt_value) == 2 and (gt_value[0] == 1 and gt_value[1] == 0):
                return None, None, None, None

        if AD_idx is not None and len(gt_value) == 2:
            ad_values = [int(i) for i in (array_record[AD_idx].split(",")[:2])]
            return second_number, sum(ad_values), ad_values[0], 0
        elif DP_idx is not None and len(gt_value) == 2:
            t_depth = int(array_record[DP_idx])
            if gt_value[0] == 0 and gt_value[0] == gt_value[1]:
                t_depth1 = t_depth - int(t_depth * (random.random() / 50) + 0.5)
                return second_number, t_depth, t_depth1, 1
            elif gt_value[0] == 1 and gt_value[0] == gt_value[1]:
                t_depth1 = int(t_depth * (random.random() / 50) + 0.5)
                return second_number, t_depth, t_depth1, 0
        return None, None, None, None

    def get_gt_DP(self):
        ref_panel = self.read_panel_npz()
        vcf_reader = VCF(self.vcf_path)
        vcf_reader.set_samples(self.sample_name)

        results = []
        sex_depth = []
        autosome_depth = []
        record_idx = 1

        for v in vcf_reader:
            vx = str(v).split("\t")
            split_var = vx[4].split(",")
            if len(split_var) > 2:
                continue

            second_number, t_depth, ref_depth, ref_flag = self.getDP(vx[7], vx[8], vx[9])
            if second_number is None or second_number < 25:
                continue
            if vx[0] in ["chrX", "chrY"]:
                sex_depth.append(t_depth)
                ref_flag = 0
            else:
                autosome_depth.append(t_depth)
            if vx[0] in ref_panel and int(vx[1]) in ref_panel[vx[0]]:
                value = ref_panel[vx[0]][int(vx[1])]
                if ref_flag:
                    record_idx += 1
                if not ref_flag or record_idx % 2 !=0:
                    value1 = int(value*self.BASE_COVERAGE + 0.5)
                    results.append((vx[0], int(vx[1]), value1, value1,  t_depth, ref_depth))

        min_sexsome_depth = np.median(sex_depth) * 0.25 if sex_depth else 0
        min_autosome_depth = np.median(autosome_depth) * 0.25 if autosome_depth else 0

        new_dtype = [('chromosome', 'U10'), ('position', int), ('value', int), ('value1', int), ('depth', int), ('refDepth', int)]
        res_array = np.array(results, dtype=new_dtype)
        return res_array, min_autosome_depth, min_sexsome_depth

    def get_result(self, raw_res, min_auto, min_sex):
        sex_array = raw_res[np.isin(raw_res['chromosome'], ['chrX', 'chrY'])]
        autosome_array = raw_res[~np.isin(raw_res['chromosome'], ['chrX', 'chrY'])]

        filtered_sex_array = sex_array[sex_array['depth'] >= min_sex]
        filtered_autosome_array = autosome_array[autosome_array['depth'] >= min_auto]
        combined_filtered_array = np.concatenate((filtered_autosome_array, filtered_sex_array))
        
        with open(self.output_path, "w") as result:
            np.savetxt(result, combined_filtered_array, delimiter='\t', fmt='%s')

    def run(self):
        raw_res, min_auto, min_sex = self.get_gt_DP()
        self.get_result(raw_res, min_auto, min_sex)


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Process VCF files')
    parser.add_argument('-v', '--vcf', help='VCF file path', required=True)
    parser.add_argument('-s', '--sample', help='Sample name', required=True)
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    options = parser.parse_args()

    gvcf_parser = gvcf2cnv(options.vcf, options.sample, options.output)
    gvcf_parser.run()


if __name__ == "__main__":
    main()
