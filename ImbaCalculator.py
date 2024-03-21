import numpy as np
import pysam
from sklearn.cluster import KMeans
import random
from argparse import ArgumentParser

class ImbaCalculator:
    def __init__(self, baf_path, interval_path, output_path):
        self.baf_path = baf_path
        self.interval_path = interval_path
        self.output_path = output_path
        random.seed(10)

    def kmeans_cluster_core(self, af):
        cluster_number = 2
        clusters_res = KMeans(n_clusters=cluster_number, init='random', algorithm='elkan', n_init=20).fit(af.reshape(-1, 1))
        value_mean = np.zeros((cluster_number,), dtype=float)
        value_size = np.zeros((cluster_number,), dtype=int)
        for i in range(cluster_number):
            value_mean[i] = np.median(af[clusters_res.labels_ == i])
            value_size[i] = af[clusters_res.labels_ == i].size
        max_size_idx = np.argsort(value_size)[::-1][:2]
        value_mean2 = value_mean[max_size_idx]

        if value_size[max_size_idx][0] / af.size > 0.93:
            return np.max(value_mean2), np.max(value_mean2) - 0.02
        else:
            return np.max(value_mean2), np.min(value_mean2)

    def kmeans_cluster(self, af, n_site):
        if n_site < 20:
            return -1

        s_af_value = np.sort(af)
        if len(np.unique(s_af_value)) <= 3:
            return np.max(s_af_value)

        af_max, af_min = self.kmeans_cluster_core(s_af_value)
        return af_min / af_max

    def cala_imba_tab(self, seg):
        tb = pysam.TabixFile(self.baf_path)
        tseg = seg.split('\t')
        n_site = 0
        init_num = max(20, int(tseg[3]) * 5)
        af_value = np.zeros(init_num, dtype=float)
        
        for record in tb.fetch(tseg[0], int(tseg[1]), int(tseg[2]), parser=pysam.asBed()):
            af_value[n_site] = abs(float(record[2]) - 0.5)
            n_site += 1
        
        if n_site == 0:
            return -1, tseg[3]
        
        return self.kmeans_cluster(af_value[:n_site], n_site), str(n_site)

    def run(self):
        with open(self.interval_path, "rt") as seg_file, open(self.output_path, "w") as result_file:
            result_file.write(f"#{seg_file.readline().strip()}\tai_value\n")
            for line in seg_file:
                imba, n_sites = self.cala_imba_tab(line.strip())
                if imba != -1:
                    line_list = line.strip().split("\t")
                    line_list[3] = n_sites
                    line_list[5] = str(int(line_list[2]) - int(line_list[1]))
                    line_content = '\t'.join(line_list)
                    result_file.write(f"{line_content}\t{imba:.3f}\n")

def main():
    parser = ArgumentParser(description='Calculate IMBA from segmented BAF data')
    parser.add_argument('-b', '--baf', help='BAF file path (baf.tab.gz)', required=True)
    parser.add_argument('-l', '--interval', help='Segment interval file path (DNAcopy.segment)', required=True)
    parser.add_argument('-o', '--output', help='Output file path (segment.imba.bed)', required=True)
    args = parser.parse_args()
    
    imba_calculator = ImbaCalculator(args.baf, args.interval, args.output)
    imba_calculator.run()

if __name__ == '__main__':
    main()
