from ngs_pipeline import arg_parser
from Bio import AlignIO
import pandas as pd
from sklearn.naive_bayes import CategoricalNB
from sklearn.neighbors import NearestCentroid
from pickle import dump
import numpy as np
from sklearn.metrics import pairwise_distances

def init_argparser():
    p = arg_parser('Generate model for Plasmodium species classification based on mitochondrial DNA')

    p.add_argument('-o', '--outfile', default='output.pickle', help='output file name (default: output.pickle)')

    p.add_argument('-k', default=10, type=int,
                   help='Size of kmer for Nearest Centroid model (default: 10)')

    p.add_argument('-d', '--distance', default="manhattan", type=str, choices=["manhattan", "euclidean"],
                     help='Distance metric for Nearest Centroid model (default: manhattan)')
    
    p.add_argument('--fit_prior', default=False, action="store_true",
                     help='Whether to learn class prior probabilities (default: False)')
    
    p.add_argument('--alpha', default=1, type=float,
                        help='Additive (Laplace/Lidstone) smoothing parameter (default: 1)')

    p.add_argument('--type', default="cnb", type=str, choices=["cnb", "nc"],
                   help='Type of model to be generated (default: cnb) - cnb: CategoricalNB, nc: NearestCentroid')

    p.add_argument('-b', '--bedfile', default=None, type=str, required = False,
                   help='bedfile to specifying region of interest')

    p.add_argument('infiles', nargs='*')
    return(p)

def build_nc_model(filtered_spec, y, k, distance, outfile, roi_start=0, roi_end=-1):
    def transform_kmer_count(sequences, k = 3):
        all_kmers = {q for p in [["".join(sequence[i:i+k]) for i in range(len(sequence) - k + 1)] for sequence in sequences] for q in p}
        results = []
        kmer = sorted(all_kmers)

        for sequence in sequences:
            kmer_count = [0]*len(kmer)
            for i in range(len(sequence) - k + 1):
                kmer_count[kmer.index("".join(sequence[i:i+k]))] += 1
            results.append(kmer_count)
        return results, kmer

    filtered_spec_kmerised, kmer = transform_kmer_count([a.seq.upper() for a in filtered_spec], k = k)
    X = filtered_spec_kmerised
    nrc = NearestCentroid(metric = distance)
    nrc.fit(X, y)
    nrc.kmer = kmer
    def transform_kmer_to_list(self, sequences, kmers=list()):
        results = []
        for sequence in sequences:
            kmer_count = [sequence.count(kmer) for kmer in kmers]
            results.append(kmer_count)
        return results
    NearestCentroid.transform_kmer_to_list = transform_kmer_to_list
    nrc.roi_start = roi_start
    nrc.roi_end = roi_end if not roi_end == -1 else len(filtered_spec[0])

    def get_full_result(self, sequences):
        sequences = [seq[self.roi_start:self.roi_end] for seq in sequences]
        result = {"read_id": [seq.id for seq in sequences]}
        X = self.transform_kmer_to_list([seq.upper() for seq in sequences], kmers = self.kmer)
        y_distance = pairwise_distances(X, self.centroids_)
        sorted_y_distance = np.sort(y_distance, axis=1)
        f = np.divide(1, sorted_y_distance[:,0], out=np.ones_like(sorted_y_distance[:,0]), where=sorted_y_distance[:,0]!=0)
        t = np.divide(1, sorted_y_distance[:,1], out=np.ones_like(sorted_y_distance[:,1]), where=sorted_y_distance[:,1]!=0)
        f_pair_diff = (2*f)/(f+t) - 1
        dist_result = {class_name: y_distance[:,class_id] for class_id, class_name in enumerate(self.classes_)}
        result["species"] = self.classes_[np.argmin(y_distance, axis=1)]
        result.update(dist_result)
        result["f_pair_diff"] = f_pair_diff
        result_df = pd.DataFrame.from_dict(result, orient = "columns")
        return result_df
    NearestCentroid.get_full_result = get_full_result

    result_model = {"model_type": "Nearest-centroid", "model":nrc, "additional_data": [kmer]}
    with open(outfile, "wb") as f:
        dump(result_model, f)


def build_cnb_model(filtered_spec, y, outfile, alpha, fit_prior, roi_start=0, roi_end=-1):
        remap_base = {"N":0, "A":1, "C":2, "G":3, "T":4, "-":5}
        filtered_spec_remap = [[remap_base[base] for base in str(a.seq).upper()] for a in filtered_spec]
        X = filtered_spec_remap
        cnb = CategoricalNB(alpha = alpha, fit_prior = fit_prior, min_categories = 6)
        cnb.fit(X, y)
        cnb.remap_base = remap_base
        def transform_base(self, sequences):
            return [[self.remap_base[base] for base in str(sequence.seq).upper()] for sequence in sequences]
        CategoricalNB.transform_base = transform_base
        cnb.roi_start = roi_start
        cnb.roi_end = roi_end if not roi_end == -1 else len(filtered_spec[0])
        def get_full_result(self, sequences):
            sequences = [seq[self.roi_start:self.roi_end] for seq in sequences]
            result = {"read_id": [seq.id for seq in sequences]}
            X = self.transform_base(sequences)
            y_pred_proba = self.predict_joint_log_proba(X)
            sorted_y_proba = np.sort(y_pred_proba, axis=1)
            f_pair_diff = sorted_y_proba[:,-1] - sorted_y_proba[:,-2]
            prob_result = {class_name: y_pred_proba[:,class_id] for class_id, class_name in enumerate(self.classes_)}
            result["species"] = self.classes_[self.argmax(y_pred_proba, axis=1)]
            result.update(prob_result)
            result["f_pair_diff"] = f_pair_diff
            result_df = pd.DataFrame.from_dict(result, orient = "columns")
            return result_df
        CategoricalNB.get_full_result = get_full_result

        result_model = {"model_type": "CategoricalNB", "model":cnb, "additional_data": [remap_base]}
        with open(outfile, "wb") as f:
            dump(result_model, f)

def main(args):
    if args.bedfile:
        roi = pd.read_csv(args.bedfile, sep="\t", header=None)
        roi_start = roi.iloc[0,1]# 0-based
        roi_end = roi.iloc[0,2] # 1-based, end is exclusive
    else:
        roi_start = 0
        roi_end = -1

    spec = AlignIO.read(args.infiles[0], "fasta")
    filtered_spec = spec[:, roi_start:roi_end]

    y = pd.Series([a.id for a in filtered_spec])
    if args.type == "nc":
        build_nc_model(filtered_spec, y, args.k, args.distance, args.outfile, roi_start, roi_end)
    elif args.type == "cnb":
        build_cnb_model(filtered_spec, y, args.outfile, args.alpha, args.fit_prior, roi_start, roi_end)