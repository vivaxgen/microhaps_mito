from sklearn.naive_bayes import CategoricalNB
from sklearn.neighbors import NearestCentroid
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from Bio.motifs import Motif
from Bio.Align import PairwiseAligner, MultipleSeqAlignment

class MitoCNB(CategoricalNB):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.roi_start = 0
        self.roi_end = -1
        self.remap_base = {"N":0, "A":1, "C":2, "G":3, "T":4, "-":5}
        self.spec_min_f_pair_diff = 0
        
    def transform_base(self, sequences):
        return [[self.remap_base[base] for base in str(sequence.seq).upper()] for sequence in sequences]
    
    def get_full_result(self, sequences):
        sequences = [seq[self.roi_start:self.roi_end] for seq in sequences]
        result = {"read_id": [seq.id for seq in sequences]}
        X = self.transform_base(sequences)
        y_pred_proba = self.predict_joint_log_proba(X)
        sorted_y_proba = np.sort(y_pred_proba, axis=1)
        f_pair_diff = sorted_y_proba[:,-1] - sorted_y_proba[:,-2]
        prob_result = {class_name: y_pred_proba[:,class_id] for class_id, class_name in enumerate(self.classes_)}
        result["species"] = self.classes_[np.argmax(y_pred_proba, axis=1)]
        result.update(prob_result)
        result["f_pair_diff"] = f_pair_diff
        result_df = pd.DataFrame.from_dict(result, orient = "columns")
        return result_df
    
    def fit2(self, X, y=None):
        filtered_spec = X[:, self.roi_start:self.roi_end]
        if y is None:
            y = pd.Series([a.id for a in filtered_spec])
        filtered_spec_remap = self.transform_base(filtered_spec)
        self.fit(filtered_spec_remap, y)
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()



class MitoNC(NearestCentroid):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.roi_start = 0
        self.roi_end = -1
        self.kmer = []
        self.k = 6
        self.spec_min_f_pair_diff = 0
    
    def transform_kmer_count(self, sequences, k = 3):
        all_kmers = {q for p in [["".join(sequence[i:i+k]) for i in range(len(sequence) - k + 1)] for sequence in sequences] for q in p}
        results = []
        kmer = sorted(all_kmers)

        for sequence in sequences:
            kmer_count = [0]*len(kmer)
            for i in range(len(sequence) - k + 1):
                kmer_count[kmer.index("".join(sequence[i:i+k]))] += 1
            results.append(kmer_count)
        return results, kmer
    
    def transform_kmer_to_list(self, sequences, kmers=list()):
        results = []
        for sequence in sequences:
            kmer_count = [sequence.count(kmer) for kmer in kmers]
            results.append(kmer_count)
        return results
    
    def get_full_result(self, sequences):
        sequences = [seq[self.roi_start:self.roi_end] for seq in sequences]
        result = {"read_id": [seq.id for seq in sequences]}
        X = self.transform_kmer_to_list([seq.upper() for seq in sequences], kmers = self.kmer)
        y_distance = pairwise_distances(X, self.centroids_)
        sorted_y_distance = np.sort(y_distance, axis=1)
        f = sorted_y_distance[:,0]
        t = sorted_y_distance[:,1]
        f_pair_diff = f - t
        dist_result = {class_name: y_distance[:,class_id] for class_id, class_name in enumerate(self.classes_)}
        result["species"] = self.classes_[np.argmin(y_distance, axis=1)]
        result.update(dist_result)
        result["f_pair_diff"] = f_pair_diff
        result_df = pd.DataFrame.from_dict(result, orient = "columns")
        return result_df
    
    def fit2(self, X, y = None):
        filtered_spec = X[:, self.roi_start:self.roi_end]
        filtered_spec_kmerised, kmer = self.transform_kmer_count([a.seq.upper() for a in filtered_spec], k = self.k)
        self.kmer = kmer
        if y is None:
            y = pd.Series([a.id for a in filtered_spec])
        self.fit(filtered_spec_kmerised, y)
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()

class MitoPA:
    def __init__(self):
        self.consensus = []
        self.species = []
        self.classes_ = []
        self.roi_start = 0
        self.roi_end = -1
        self.spec_min_f_pair_diff = 0
        self.pa = PairwiseAligner(scoring = "blastn")

    
    def get_full_result(self, sequences):
        sequences = [seq[self.roi_start:self.roi_end] for seq in sequences]
        result = {"read_id": [seq.id for seq in sequences]}
        X = [seq.seq.replace("-", "").upper() for seq in sequences]

        pa_score = np.array([[self.pa.score(sp.replace("-", "").upper(), seq) for sp in self.consensus] for seq in X])
        sorted_pa_score = np.sort(pa_score, axis=1)
        f = sorted_pa_score[:,-1]
        t = sorted_pa_score[:,-2]
        f_pair_diff = np.abs(f - t)
        dist_result = {class_name: pa_score[:,class_id] for class_id, class_name in enumerate(self.classes_)}
        result["species"] = self.classes_[np.argmax(pa_score, axis=1)]
        result.update(dist_result)
        result["f_pair_diff"] = f_pair_diff
        result_df = pd.DataFrame.from_dict(result, orient = "columns")
        return result_df
    
    def fit2(self, X, y = None):
        if y is None:
            y = pd.Series([a.id for a in X])
        filtered_spec = X[:, self.roi_start:self.roi_end]
        for uy in y.unique():
            subspec = MultipleSeqAlignment([a for a in filtered_spec if a.id == uy])
            mot = Motif('ACGTN-', subspec.alignment)
            cons = mot.degenerate_consensus
            self.consensus.append(cons)
            self.classes_.append(uy)
            self.species.append(uy)
        self.classes_ = np.array(self.classes_)
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()

        