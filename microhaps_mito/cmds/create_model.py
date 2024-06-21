from ngs_pipeline import arg_parser
from Bio import AlignIO
import pandas as pd
from mito_models import MitoCNB, MitoNC, MitoPA
from pickle import dump


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

    p.add_argument('--type', default="cnb", type=str, choices=["cnb", "nc", "pa", "cons"],
                   help='Type of model to be generated (default: cnb) - cnb: CategoricalNB, nc: NearestCentroid, pa: PairwiseAligner, cons: get consensus for mapping only')

    p.add_argument('-b', '--bedfile', default=None, type=str, required = False,
                   help='bedfile to specifying region of interest')

    p.add_argument('-m', '--meta', default=None, type=str, required = False,
                   help='metadata file for species label, single column with no header (default: None), same number as defined in fasta file')

    p.add_argument('infiles', nargs='*')
    return(p)

def build_nc_model(spec, y, k, distance, outfile, roi_start=0, roi_end=-1):
    nrc = MitoNC(metric = distance)
    nrc.roi_start = roi_start
    nrc.roi_end = roi_end if not roi_end == -1 else len(spec[0])
    nrc.k = k
    nrc.fit2(spec, y)

    with open(outfile, "wb") as f:
        dump(nrc, f)


def build_cnb_model(spec, y, outfile, alpha, fit_prior, roi_start=0, roi_end=-1):
    cnb = MitoCNB(alpha = alpha, fit_prior = fit_prior)
    cnb.roi_start = roi_start
    cnb.roi_end = roi_end if not roi_end == -1 else len(spec[0])
    cnb.fit2(spec, y)
    with open(outfile, "wb") as f:
        dump(cnb, f)

def build_pa_model(spec, y, outfile, roi_start=0, roi_end=-1):
    pa = MitoPA()
    pa.roi_start = roi_start
    pa.roi_end = roi_end if not roi_end == -1 else len(spec[0])
    pa.fit2(spec, y)
    with open(outfile, "wb") as f:
        dump(pa, f)

def build_consensus_seq(spec, outfile, roi_start=0, roi_end=-1):
    from Bio.motifs import Motif
    roi_end = roi_end if not roi_end == -1 else len(spec[0])
    spec = spec[:, roi_start:roi_end]
    mot = Motif('ACGTN-', spec.alignment)
    mot.counts["-"] = [0 for _ in range(mot.length)]
    cons = mot.consensus
    with open(outfile, "w") as f:
        f.write(">consensus\n")
        f.write(str(cons))
        f.write("\n")


def main(args):
    if args.bedfile:
        roi = pd.read_csv(args.bedfile, sep="\t", header=None)
        roi_start = roi.iloc[0,1]# 0-based
        roi_end = roi.iloc[0,2] # 1-based, end is exclusive
    else:
        roi_start = 0
        roi_end = -1

    spec = AlignIO.read(args.infiles[0], "fasta")
    if args.meta:
        y = pd.Series([a for a in open(args.meta).read().split("\n")])
    else:
        y = pd.Series([a.id for a in spec])
    
    assert len(spec) == len(y), "Number of sequences in fasta and meta file do not match"
    
    if args.type == "nc":
        build_nc_model(spec, y, args.k, args.distance, args.outfile, roi_start, roi_end)
    elif args.type == "cnb":
        build_cnb_model(spec, y, args.outfile, args.alpha, args.fit_prior, roi_start, roi_end)
    elif args.type == "pa":
        build_pa_model(spec, y, args.outfile, roi_start, roi_end)
    elif args.type == "cons":
        build_consensus_seq(spec, args.outfile, roi_start, roi_end)