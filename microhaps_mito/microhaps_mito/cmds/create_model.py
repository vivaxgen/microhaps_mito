from ngs_pipeline import arg_parser
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
from microhaps_mito.mito_models import MitoCNB, MitoNC, MitoPA, MitoEnsemble, cascadingSpeciationModel
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
                   help='metadata file for species label, single column with no header (default: None), same number as defined in fasta file')

    p.add_argument('-c', '--cascade', default = False, action="store_true", 
                   help=("Whether to cascade the model building process (default: False), submodels will be generated, e.g."
                         "`Plasmodium Vivax, Plasmodium Falciparum, Plasmodium Ovale Curtisi, Plasmodium Ovale Wallikeri, Plasmodium Malariae` will generate:"
                         "model1 to determine --> Vivax, Falciparum Ovale Malariae, model2 to determine --> Ovale Curtisi, Ovale Wallikeri")
    )

    p.add_argument('infiles', nargs='*')
    return(p)

def build_nc_model(spec, y, k, distance, roi_start=0, roi_end=-1, cascade = False):
    if cascade:
        nrc = cascadingSpeciationModel(
            (MitoNC, {"metric": distance, "k": k, "roi_start": roi_start, "roi_end": roi_end}))
        nrc.fit2(spec, y)
    else:
        nrc = MitoNC(metric = distance, k = k, roi_start = roi_start, roi_end = roi_end)
        nrc.fit2(spec, y)
    return nrc


def build_cnb_model(spec, y, alpha, fit_prior, roi_start=0, roi_end=-1, cascade = False):
    if cascade:
        cnb = cascadingSpeciationModel(
            (MitoCNB, {"alpha": alpha, "fit_prior": fit_prior, "roi_start": roi_start, "roi_end": roi_end}))
        cnb.fit2(spec, y)
    else:
        cnb = MitoCNB(alpha = alpha, fit_prior = fit_prior, roi_start = roi_start, roi_end = roi_end)
        cnb.fit2(spec, y)
    return cnb

def build_pa_model(spec, y, roi_start=0, roi_end=-1, cascade = False):
    if cascade:
        pa = cascadingSpeciationModel(
            (MitoPA, {"roi_start": roi_start, "roi_end": roi_end}))
        pa.fit2(spec, y)
    else:
        pa = MitoPA(roi_start = roi_start, roi_end = roi_end)
        pa.fit2(spec, y)
    return pa

def build_ensemble_model(spec, y, k, distance, alpha, fit_prior, roi_start=0, roi_end=-1, cascade = False):
    if cascade:
        ensemble = cascadingSpeciationModel(
            (MitoEnsemble, {
                "models_n_args": [
                    (MitoCNB, {"alpha": alpha, "fit_prior": fit_prior, "roi_start": roi_start, "roi_end": roi_end}),
                    (MitoNC, {"metric": distance, "k": k, "roi_start": roi_start, "roi_end": roi_end}),
                    (MitoPA, {"roi_start": roi_start, "roi_end": roi_end})
                    ]}),
            roi_start = roi_start, roi_end = roi_end
        )
        ensemble.fit2(spec, y)
    else:
        ensemble = MitoEnsemble(models_n_args = [(MitoCNB, {"alpha": alpha, "fit_prior": fit_prior, "roi_start": roi_start, "roi_end": roi_end}),
                              (MitoNC, {"metric": distance, "k": k, "roi_start": roi_start, "roi_end": roi_end}),
                              (MitoPA, {"roi_start": roi_start, "roi_end": roi_end})],
                              roi_start = roi_start, roi_end = roi_end)
        ensemble.fit2(spec, y)
    return ensemble

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

    spec = AlignIO.read(args.infiles[0], "fasta")
    if args.meta:
        y = pd.Series([a for a in open(args.meta).read().strip().split("\n")])
    else:
        y = pd.Series([a.id for a in spec])
    
    assert len(spec) == len(y), "Number of sequences in fasta and meta file do not match"

    def run_build_model(subspec, suby):
        if args.type == "nc":
            res = build_nc_model(subspec, suby, args.k, args.distance, roi_start, roi_end, cascade = args.cascade)
        elif args.type == "cnb":
            res = build_cnb_model(subspec, suby, args.alpha, args.fit_prior, roi_start, roi_end, cascade = args.cascade)
        elif args.type == "pa":
            res = build_pa_model(subspec, suby, roi_start, roi_end, cascade = args.cascade)
        elif args.type == "ensemble":
            res = build_ensemble_model(subspec, suby, args.k, args.distance, args.alpha, args.fit_prior, roi_start, roi_end, cascade = args.cascade)
        return res
    
    if args.type == "cons":
        res = build_consensus_seq(spec, roi_start, roi_end)
        with open(args.outfile, "w+") as f:
            f.write(res)
        return
    else:
        res = run_build_model(spec, y)
        write_to_pickle(res)