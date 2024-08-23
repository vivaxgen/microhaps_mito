from sklearn.naive_bayes import CategoricalNB
from sklearn.neighbors import NearestCentroid
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from Bio.motifs import Motif
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from itertools import zip_longest
import copy

class MitoCNB(CategoricalNB):
    def __init__(self, alpha = 1, fit_prior = False, roi_start=0, roi_end=-1, remap_base = {"N":0, "A":1, "C":2, "G":3, "T":4, "-":5}):
        super().__init__(alpha = alpha, fit_prior = fit_prior, min_categories = len(remap_base))
        self.roi_start = roi_start
        self.roi_end = roi_end
        self.remap_base = remap_base
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
        if self.roi_end == -1:
            self.roi_end = len(X[0])
        filtered_spec = X[:, self.roi_start:self.roi_end]
        if y is None:
            y = pd.Series([a.id for a in filtered_spec])
        filtered_spec_remap = self.transform_base(filtered_spec)
        self.fit(filtered_spec_remap, y)
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()



class MitoNC(NearestCentroid):
    def __init__(self, metric = "euclidean", roi_start=0, roi_end=-1, kmer = [], k=6):
        super().__init__(metric = metric)
        self.roi_start = roi_start
        self.roi_end = roi_end
        self.kmer = kmer
        self.k = k
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
        if self.roi_end == -1:
            self.roi_end = len(X[0])
        filtered_spec = X[:, self.roi_start:self.roi_end]
        filtered_spec_kmerised, kmer = self.transform_kmer_count([a.seq.upper() for a in filtered_spec], k = self.k)
        self.kmer = kmer
        if y is None:
            y = pd.Series([a.id for a in filtered_spec])
        self.fit(filtered_spec_kmerised, y)
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()

class MitoPA:
    def __init__(self, roi_start=0, roi_end=-1):
        self.consensus = []
        self.species = []
        self.classes_ = []
        self.roi_start = roi_start
        self.roi_end = roi_end
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
        if self.roi_end == -1:
            self.roi_end = len(X[0])
        filtered_spec = X[:, self.roi_start:self.roi_end]
        for uy in np.unique(y):
            subspec = [filtered_spec[i,:] for i, v in enumerate(y) if v == uy]
            subspec = MultipleSeqAlignment(subspec)
            mot = Motif('ACGTN-', subspec.alignment)
            cons = mot.degenerate_consensus
            self.consensus.append(cons)
            self.classes_.append(uy)
            self.species.append(uy)

        self.classes_ = np.array(self.classes_)
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()


class MitoEnsemble:
    def __init__(self, models=[], models_n_args = None, roi_start = 0, roi_end = -1):
        # models_n_args = [(model_class, args), (model_class, args), ...]
        self.models = models if models_n_args is None else [model(**args) for model, args in models_n_args]
        self.roi_start = roi_start
        self.roi_end = roi_end
        self.spec_min_f_pair_diff = 0
        self.classes_ = models[0].classes_ if len(models) > 0 else []
    
    def fit2(self, X, y = None):
        if y is None:
            y = pd.Series([a.id for a in X])
        if self.roi_end == -1:
            self.roi_end = len(X[0])
        for model in self.models:
            model.roi_start = self.roi_start
            model.roi_end = self.roi_end
            model.fit2(X, y)
        self.classes_ = self.models[0].classes_
        self.spec_f_pair_diff = self.get_full_result(X)["f_pair_diff"]
        self.spec_min_f_pair_diff = self.spec_f_pair_diff.min()

    def get_full_result(self, sequences):
        result = {"read_id": [seq.id for seq in sequences]}
        scores = None
        for model in self.models:
            result_df = model.get_full_result(sequences)
            if (model.classes_ != self.classes_).all():
                raise ValueError("Species classes do not match")
            if model.__class__.__name__ == "MitoNC":
                weighted = result_df[self.classes_]
                weighted = 1 - (weighted - weighted.min())/(weighted.max() - weighted.min())
            else:
                weighted = result_df[self.classes_]
                weighted = (weighted - weighted.min())/(weighted.max() - weighted.min())
            if scores is None:
                scores = weighted
            else:
                scores += weighted
        scores = scores/len(self.models)
        result["species"] = self.classes_[np.argmax(scores.values, axis=1)]
        result.update({class_name: scores.iloc[:,class_id] for class_id, class_name in enumerate(self.classes_)})
        sorted_scores = np.sort(scores.values, axis=1)
        result["f_pair_diff"] = sorted_scores[:,-1] - sorted_scores[:,-2]
        return pd.DataFrame.from_dict(result, orient = "columns")
    
class cascadingSpeciationModel:
    def __init__(self, models_n_args = None, roi_start = 0, roi_end = -1):
        assert len(models_n_args) == 2, "A model is needed"
        assert models_n_args[0].__name__ in ["MitoCNB", "MitoNC", "MitoPA", "MitoEnsemble"], "Invalid model"
        self.cascade_models = {}
        self.model = models_n_args[0](**models_n_args[1])
        self.roi_start = roi_start
        self.roi_end = roi_end
        self.spec_min_f_pair_diff = {}
        self.classes_ = {}

        # [{"Species": "Plasmodium", "model": XXX, "submodels": [{"Species": "Plasmodium Vivax", "model":}, {}, {}]}]
    def fit2(self, X, y = None):
        if y is None:
            y_name = pd.Series([a.id for a in X])
            y = [a.split(" ") for a in y_name]
        
        # y = [["Plasmodium", "Vivax"], ["Plasmodium", "Falciparum"], ["Plasmodium", "Malariae"], ["Plasmodium", "Ovale", "Curtisi"], ["Plasmodium", "Ovale", "Wallikeri"]]
        # Determine how many submodels are needed
        y = [a.split(" ") for a in y]
        y = np.array([list(tpl) for tpl in zip(*zip_longest(*y, fillvalue=""))])

        # build_search_map
        # Ensure all the samples have common "genus"
        assert len(set(y[:,0])) == 1, "All samples must have the same genus"
        
        for spp_i in range(1, y.shape[1]):
            to_diff = np.unique([yi[0:spp_i] for yi in y if yi[spp_i] != "" ], axis = 0)
            if to_diff.shape[1] < 1:
                continue
            else:
                for sub_diff in to_diff:
                    filtered_y = [n[spp_i] for n in y[:,0:spp_i+1] if n[spp_i]!= "" and (n[0:spp_i] == sub_diff).all()]
                    filtered_y_index = [i for i, v in enumerate(y[:,spp_i]) if v in filtered_y]
                    filtered_X = MultipleSeqAlignment([X[i] for i in filtered_y_index])
                    model = copy.deepcopy(self.model)
                    filtered_y_fullname = [" ".join(n) for n in y[:,0:spp_i+1] if n[spp_i]!= "" and (n[0:spp_i] == sub_diff).all()]
                    model.fit2(filtered_X, filtered_y_fullname)
                    if len(self.cascade_models) == 0:
                        self.cascade_models[0] = model
                        self.classes_[0] = model.classes_
                        self.spec_min_f_pair_diff[0] = model.spec_min_f_pair_diff
                    else:
                        submodel_name = " ".join(sub_diff)
                        self.cascade_models[submodel_name] = model
                        self.classes_[submodel_name] = model.classes_
                        self.spec_min_f_pair_diff[submodel_name] = model.spec_min_f_pair_diff

    def get_full_result(self, sequences):
        # get result from the primary model and expand if needed
        result = self.cascade_models[0].get_full_result(sequences)

        for _class in self.classes_.keys():
            if _class == 0:
                continue
            
            subset_indexes = np.where(result["species"] == _class)[0]
            if len(subset_indexes) > 0:
                model = self.cascade_models[_class]
                subset_sequences = MultipleSeqAlignment([sequences[int(i)] for i in subset_indexes])
                subset_result = model.get_full_result(subset_sequences)
                subset_result.rename(columns = {"species": f"{_class} subspecies", "f_pair_diff": f"{_class}_f_pair_diff"}, inplace = True)
                result = pd.merge(result, subset_result, how = "left", on = "read_id")
        return result

        # subset the result to those that need submodels
        
        # get result with submodels


    # {"Plasmodium": model1, "Plasmodium Ovale": model2} --> while Result in dict.keys()


        # # Train the primary model
        # current_model = self.model
        # current_model.fit2(X, y[:,0])
        # self.cascade_models = {0: self.model}
        # to_diff = np.unique([yi[0:y.shape[1]-1] for yi in y if yi[y.shape[1]-1] != "" ], axis = 0)
        # spp_i = y.shape[1]-1
        # if to_diff.shape[1] < 1:
        #     raise ValueError("No submodels needed")
        # else:
        #     for sub_diff in to_diff:
        #         filtered_y = [n[spp_i] for n in y[:,0:spp_i+1] if n[spp_i]!= "" and (n[0:spp_i] == sub_diff).all()]
        #         filtered_y_index = [i for i, v in enumerate(y[:,spp_i]) if v in filtered_y]
        #         filtered_X = X[filtered_y_index]
        #         model = MitoEnsemble(models = self.models)
        #         model.fit2(filtered_X, filtered_y)
        #         self.cascade_models.append({"species": sub_diff, "model": model})


# "pk", "pv", "pv" "pmo","pmv", "pmm", "pva","pvb"
# result = result = {"read_id": [1,2,3,4,5],
#                    "species": ["pv", "pmo", "pmv", "pk", "pf"],}
# pd.DataFrame.from_dict(result, orient = "columns")