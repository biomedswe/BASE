import numpy as np
import pandas as pd
import pysam
import random
from argparse import ArgumentParser

class LogRAdjuster:
    def __init__(self, logr_path, segment_path, output_path, ploidy=2):
        self.logr_path = logr_path
        self.segment_path = segment_path
        self.output_path = output_path
        self.ploidy = ploidy
        random.seed(10)

    def read_segments(self):
        column_names = ['seqnames', 'starts', 'ends', 'Markers', 'Value', 'Size', 'ai_value']
        segments = pd.read_csv(self.segment_path, sep='\s+', comment='#', names=column_names, 
                               dtype={'seqnames': str, 'starts': int, 'ends': int, 'Markers': int, 
                                      'Value': float, 'Size': int, 'ai_value': float})
        return segments

    def disomy_logr(self, segments, logr_value_list):
        tb = pysam.TabixFile(self.logr_path)
        idx = 0
        for row in segments.itertuples(index=False):
            for record in tb.fetch(row.seqnames, row.starts, row.ends, parser=pysam.asBed()):    
                logr_value_list[idx] = float(record[2])
                idx += 1

    def adjust_logR(self):
        segments = self.read_segments()
        logr_values = np.zeros((segments['Markers'].sum(), ), dtype=float)
        
        disomy_segments = segments[(segments['ai_value'] <= 0.2) & (segments['Value'] < 0.2)]
        use_segments = disomy_segments if len(disomy_segments) > 3 else segments

        self.disomy_logr(use_segments, logr_values)

        median_value = np.median(logr_values) + 0.5 * (2 - int(self.ploidy))  # Directly use np.median here
        segments['adj_logr'] = segments['Value'] - median_value

        title_line = "seqnames\tstarts\tends\tMarkers\tValue\tSize\tai_value\tadj_logr"
        with open(self.output_path, 'w') as f:
            f.write(f"#{title_line}\n")  # "#" and title line are now combined
            segments.to_csv(f, mode='a', sep='\t', index=False, float_format='%.3f', header=False)


def main():
    parser = ArgumentParser(description='Adjust logR values based on median logR and ploidy.')
    parser.add_argument('-r', '--logr', required=True, help='Path to logr file (logr.tab.gz).')
    parser.add_argument('-l', '--interval', required=True, help='Path to interval file (DNAcopy.segment).')
    parser.add_argument('-p', '--ploidy', default=2, type=int, help='Ploidy level (default: 2).')
    parser.add_argument('-o', '--output', required=True, help='Path to output file (segment.imba.bed).')

    args = parser.parse_args()

    adjuster = LogRAdjuster(logr_path=args.logr, segment_path=args.interval, output_path=args.output, ploidy=args.ploidy)
    adjuster.adjust_logR()


if __name__ == '__main__':
    main()