from ngs_pipeline import arg_parser
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
from microhaps_mito.mito_models import MitoCNB, MitoNC, MitoPA, MitoEnsemble
from pickle import dump
import numpy as np
from itertools import zip_longest

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

    p.add_argument('--type', default="cnb", type=str, choices=["cnb", "nc", "pa", "ensemble", "cons"],
                   help='Type of model to be generated (default: cnb) - cnb: CategoricalNB, nc: NearestCentroid, pa: PairwiseAligner, ensemble: combination of all previous three, cons: get consensus for mapping only')

    p.add_argument('-b', '--bedfile', default=None, type=str, required = False,
                   help='bedfile to specifying region of interest')

    p.add_argument('-m', '--meta', default=None, type=str, required = False,
                   help='metadata file for species label, two columns with `fasta_name` and `species` as header (default: None), same number as defined in fasta file')

    p.add_argument('-c', '--cascade', default = False, action="store_true", 
                   help=("Whether to cascade the model building process (default: False), submodels will be generated, e.g."
                         "`Plasmodium Vivax, Plasmodium Falciparum, Plasmodium Ovale Curtisi, Plasmodium Ovale Wallikeri, Plasmodium Malariae` will generate:"
                         "model1 to determine --> Vivax, Falciparum Ovale Malariae, model2 to determine --> Ovale Curtisi, Ovale Wallikeri")
    )

    p.add_argument('infiles', nargs='*')
    return(p)

def build_nc_model(spec, y, k, distance, roi_start=0, roi_end=-1):
    nrc = MitoNC(metric = distance, k = k, roi_start = roi_start, roi_end = roi_end)
    nrc.fit2(spec, y)
    return nrc


def build_cnb_model(spec, y, alpha, fit_prior, roi_start=0, roi_end=-1):
    cnb = MitoCNB(alpha = alpha, fit_prior = fit_prior, roi_start = roi_start, roi_end = roi_end)
    cnb.fit2(spec, y)
    return cnb

def build_pa_model(spec, y, roi_start=0, roi_end=-1):
    pa = MitoPA(roi_start = roi_start, roi_end = roi_end)
    pa.fit2(spec, y)
    return pa

def build_consensus_seq(spec, roi_start=0, roi_end=-1):
    from Bio.motifs import Motif
    roi_end = roi_end if not roi_end == -1 else len(spec[0])
    spec = spec[:, roi_start:roi_end]
    mot = Motif('ACGTN-', spec.alignment)
    mot.counts["-"] = [0 for _ in range(mot.length)]
    cons = mot.consensus
    outstr = f">consensus\n{str(cons)}\n"
    return outstr


def main(args):
    
    def write_to_pickle(obj):
        with open(args.outfile, "wb+") as f:
            dump(obj, f)
    
    if args.bedfile:
        roi = pd.read_csv(args.bedfile, sep="\t", header=None)
        roi_start = roi.iloc[0,1]# 0-based
        roi_end = roi.iloc[0,2] # 1-based, end is exclusive
    else:
        roi_start = 0
        roi_end = -1
    
    def cascadeSpecMeta(y, spec):
        result = []
        y_name = [a.split(" ") for a in y]
        ny = np.array([list(tpl) for tpl in zip(*zip_longest(*y_name, fillvalue=""))])
        for spp_i in range(ny.shape[1]):
            to_diff = np.unique([yi[0:spp_i] for yi in ny if yi[spp_i] != "" ], axis = 0)
            if to_diff.shape[1] < 1:
                continue
            else:
                for sub_diff in to_diff:
                    filtered_y = [n[spp_i] for n in ny[:,0:spp_i+1] if n[spp_i]!= "" and (n[0:spp_i] == sub_diff).all()]
                    filtered_y_index = [i for i, v in enumerate(ny[:,spp_i]) if v in filtered_y]
                    filtered_spec = [spec[i,:] for i in filtered_y_index]
                    filtered_spec = MultipleSeqAlignment(filtered_spec)
                    result.append((" ".join(sub_diff), np.array(filtered_y), filtered_spec))
        return result

    spec = AlignIO.read(args.infiles[0], "fasta")
    if args.meta:
        y = pd.Series([a for a in open(args.meta).read().strip().split("\n")])
    else:
        y = pd.Series([a.id for a in spec])
    
    assert len(spec) == len(y), "Number of sequences in fasta and meta file do not match"

    def run_build_model(subspec, suby):
        if args.type == "nc":
            res = build_nc_model(subspec, suby, args.k, args.distance, roi_start, roi_end)
        elif args.type == "cnb":
            res = build_cnb_model(subspec, suby, args.alpha, args.fit_prior, roi_start, roi_end)
        elif args.type == "pa":
            res = build_pa_model(subspec, suby, roi_start, roi_end)
        elif args.type == "ensemble":
            cnb = build_cnb_model(subspec, suby, args.alpha, args.fit_prior, roi_start, roi_end)
            nc = build_nc_model(subspec, suby, args.k, args.distance, roi_start, roi_end)
            pa = build_pa_model(subspec, suby, roi_start, roi_end)
            res = MitoEnsemble([cnb, nc, pa])
        return res
    
    if args.type == "cons":
        res = build_consensus_seq(spec, roi_start, roi_end)
        with open(args.outfile, "w+") as f:
            f.write(res)
        return

    res = {}
    if args.cascade:
        cascade = cascadeSpecMeta(y, spec)
        for name, suby, subspec in cascade:
            m = run_build_model(subspec, suby)
            if len(res.keys()) == 0:
                res[0] = m
            else:
                res[name] = m
        res["is_cascade"] = True
        write_to_pickle(res)
    else:
        #import IPython; IPython.embed()
        res[0] = run_build_model(spec, y)
        res["is_cascade"] = False
        write_to_pickle(res)