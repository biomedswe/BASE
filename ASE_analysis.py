import numpy as np
from scipy.stats import binomtest
import pysam
import copy
from argparse import ArgumentParser

class ASEAnalysis:
    class GeneR:
        def __init__(self, ghap1, ghap2, rhap1, rhap2):
            self.ghap1 = ghap1
            self.ghap2 = ghap2
            self.rhap1 = rhap1
            self.rhap2 = rhap2
            self.tcn = []
            self.h1cn = []

        def __iadd__(self, other):
            self.ghap1 += other.ghap1
            self.ghap2 += other.ghap2
            self.rhap1 += other.rhap1
            self.rhap2 += other.rhap2
            self.tcn.extend(other.tcn)
            self.h1cn.extend(other.h1cn)
            return self

        def add_h1cn_item(self, item):
            self.h1cn.append(item)

        def add_tcn_item(self, item):
            self.tcn.append(item)

        def get_median(self):
            return np.median(self.tcn), np.median(self.h1cn), len(self.tcn)

        def deepcopy(self):
            new_instance = ASEAnalysis.GeneR(self.ghap1, self.ghap2, self.rhap1, self.rhap2)
            new_instance.tcn = copy.deepcopy(self.tcn)
            new_instance.h1cn = copy.deepcopy(self.h1cn)
            return new_instance


    def ase_test(self, geneR_item, t_cn, h1_cn):
        r_total = geneR_item.rhap1 + geneR_item.rhap2
        g_total = geneR_item.ghap1 + geneR_item.ghap2
        g_frac = geneR_item.ghap1 / g_total
        p_wxs = binomtest(geneR_item.rhap1, r_total, g_frac)
        p_balanced = binomtest(geneR_item.rhap1, r_total, 0.5)
    
        rhap2_adjusted = max(geneR_item.rhap2, 1)
        ase_bias = geneR_item.rhap1 / rhap2_adjusted
    
        g_bias = 1.0
        p_ai_adj = p_balanced
        h2_cn = t_cn - h1_cn
    
        if h1_cn == 0:
            g_bias = 0.7 if geneR_item.ghap1 >= geneR_item.ghap2 else 0.3
            p_ai_adj = binomtest(geneR_item.rhap1, r_total, 0.95 if g_bias == 0.7 else 0.05)
        elif h1_cn != h2_cn:
            ratio = h2_cn / h1_cn if geneR_item.ghap1 >= geneR_item.ghap2 else h1_cn / h2_cn
            g_bias = ratio
            p_ai_adj = binomtest(geneR_item.rhap1, r_total, 1 - h1_cn / t_cn if geneR_item.ghap1 >= geneR_item.ghap2 else h1_cn / t_cn)
    
        ase_cnv_adj = ase_bias / g_bias
        return p_wxs.pvalue, p_balanced.pvalue, p_ai_adj.pvalue, ase_bias, ase_cnv_adj


    def __init__(self, input_file, bed_file, output_file):
        self.input_file = input_file
        self.bed_file = bed_file
        self.output_file = output_file

    def run_analysis(self):
        gene_dict = {}
        gene_name = {}
        gene_site = {}
        gene_chr = {}

        with open(self.input_file, "r") as ase, open(self.output_file, "w") as result:
            result.write("Gene_Symbol\tENS_ID\tChr\tSites\tgHap1\tgHap2\trHap1\trHap2\ttotal_cn\thom1_cn\tNsites\tP_by_WXS\tP_balanced\tP_imba\tASE_ratio\tCNadj_ASE_ratio\n")
            ase.readline()  # Skip header
            tbx = pysam.TabixFile(self.bed_file)
            for line in ase:
                line = line.strip()
                tseg = line.split('\t')
                for record in tbx.fetch(tseg[0], (int(tseg[1])-1), int(tseg[1]), parser=pysam.asBed()):
                    gene_name[record[3]] = record[4]
                    gene_chr[record[3]] = tseg[0]
                    gene_site.setdefault(record[3], []).append(int(tseg[1]) - 1)

                    tmp = ASEAnalysis.GeneR(int(tseg[13]), int(tseg[14]), max(int(tseg[5]), int(tseg[6])), min(int(tseg[5]), int(tseg[6])))
                    tmp.add_tcn_item(float(tseg[15]))
                    tmp.add_h1cn_item(float(tseg[16]))
                    if record[3] in gene_dict:
                        gene_dict[record[3]] += tmp
                    else:
                        gene_dict[record[3]] = tmp.deepcopy()

            for key, value in gene_dict.items():
                total_cn, hom1_cn, nsites = value.get_median()
                if value.ghap2 == 0 or total_cn < 1 or hom1_cn < 0:
                    continue
                P_wxs, P_balanced, P_ai, ASE_bias, ASE_cnv_adj = self.ase_test(value, total_cn, hom1_cn)
                gene_sites = '|'.join(map(str, gene_site[key]))
                result.write(f"{gene_name[key]}\t{key}\t{gene_chr[key]}\t{gene_sites}\t{value.ghap1}\t{value.ghap2}\t{value.rhap1}\t{value.rhap2}\t{int(total_cn)}\t{int(hom1_cn)}\t{nsites}\t{P_wxs:.16f}\t{P_balanced:.16f}\t{P_ai:.16f}\t{ASE_bias:.3f}\t{ASE_cnv_adj:.3f}\n")


if __name__ == "__main__":
    parser = ArgumentParser(description="ASE Analysis")
    parser.add_argument("-i", "--input", help="RNA_ase_cn file path", required=True)
    parser.add_argument("-b", "--bed", help="Bed file path", required=True)
    parser.add_argument("-o", "--output", help="Output file path", required=True)
    args = parser.parse_args()

    analyzer = ASEAnalysis(args.input, args.bed, args.output)
    analyzer.run_analysis()

