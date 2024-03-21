import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift
import joblib
from argparse import ArgumentParser
from shortcuts import Shortcuts

class CNVTrainerPredictor:
    def __init__(self):
        self.shortcuts = Shortcuts()
        self.training_set_path = self.shortcuts.reference_cnv_trainingset_raw
        self.database_path = self.shortcuts.reference_cnv_trainingset

    def read_training_set_idx(self):
        training_set = pd.read_csv(self.training_set_path, dtype={0: int, 2: float, 3: float, 4: float}, compression='gzip', sep='\t', header=0)
        return training_set[["adj_logr", "ai_value"]].to_numpy()

    def MeanShiftClustering(self):
        training_input = self.read_training_set_idx()
        ms = MeanShift(bandwidth=0.1111, bin_seeding=True, n_jobs=-1)
        ms.fit(training_input)

        labels = ms.labels_
        n_clusters = len(np.unique(labels))
        colors = plt.cm.rainbow(np.linspace(0, 1, n_clusters))

        plt.figure(figsize=(11, 9))
        for k, col in zip(range(n_clusters), colors):
            cluster_members = labels == k
            plt.scatter(training_input[cluster_members, 0], training_input[cluster_members, 1], color=col, s=10, alpha=0.5)
        plt.title("MeanShift Clustering")
        plt.savefig("MeanShift_Clustering.pdf", format='pdf')

        joblib.dump(ms, self.shortcuts.reference_cnv_trainingset)  # Save the model using shortcuts

    def cnv_predict(self, segment_file, output_file):
        calling_model = joblib.load(self.database_path)
        min_segment = 100000

        with open(segment_file, "r") as dnacopy_ai, open(output_file, "w") as result:
            result.write(f"{dnacopy_ai.readline().strip()}\ttotal_cn\thom1_cn\n")

            for line in dnacopy_ai:
                line = line.strip().split("\t")
                if int(line[5]) < min_segment:
                    continue

                seg_data = np.array([[float(line[7]) - 0.005, float(line[6])]], dtype=float)
                predicted_idx = calling_model.predict(seg_data)[0]
                _, cnv_stat = self.map_cnv_predict_res(predicted_idx)

                joined_line = '\t'.join(line)
                result.write(f"{joined_line}\t{cnv_stat[0]}\t{cnv_stat[1]}\n")


    def map_cnv_predict_res(self, idx):
        cnv_state = {
            0: [2, 1], 1: [3, 1], 2: [4, 2], 3: [2, 0], 4: [1, 0], 5: [1, 0], 6: [2, 0],
            7: [2, 0], 8: [4, 1], 9: [1, 0], 10: [5, 2], 11: [2, 0], 12: [3, 0],
            13: [4, 1], 14: [4, 0], 15: [6, 2], 16: [5, 1],
        }
        return idx, cnv_state.get(idx, [0, 0])



def main():
    parser = ArgumentParser(description='AI_CNV Training and Prediction')
    parser.add_argument('-d', '--database', action='store_true', help='Flag to indicate training using the default training set')
    parser.add_argument('-s', '--segment', required=False, help='Segment file path for prediction')
    parser.add_argument('-o', '--output', required=False, help='Output file path for prediction results')
    args = parser.parse_args()

    cnv_trainer_predictor = CNVTrainerPredictor()

    # Training mode using the default training set
    if args.database:
        print("Training using the default training set...")
        cnv_trainer_predictor.MeanShiftClustering()

    # Prediction mode using custom segment and output paths
    elif args.segment and args.output:
        print("Predicting using custom segment and output files...")
        cnv_trainer_predictor.cnv_predict(args.segment, args.output)
    
    else:
        print("Insufficient arguments provided. Use -d for training with the default training set, or -s and -o for custom prediction.")

if __name__ == '__main__':
    main()
