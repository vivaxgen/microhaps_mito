from ngs_pipeline import arg_parser

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

    p.add_argument('-b', '--bedfile', default=None, type=str,
                   help='bedfile to specifying region of interest')

    p.add_argument('infiles', nargs='*')
    return(p)

def main(args):
    from Bio import AlignIO
    import pandas as pd
    from sklearn.naive_bayes import CategoricalNB
    from sklearn.neighbors import NearestCentroid
    from pickle import dump
    import numpy as np
    roi = pd.read_csv(args.bedfile, sep="\t", header=None)
    roi_start = roi.iloc[0,1]# 0-based
    roi_end = roi.iloc[0,2] # 1-based, end is exclusive

    spec = AlignIO.read(args.infiles[0], "fasta")
    filtered_spec = spec[:, roi_start:roi_end]

    y = pd.Series([a.id for a in filtered_spec])

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

    if args.type == "nc":
        filtered_spec_kmerised, kmer = transform_kmer_count([a.seq.upper() for a in filtered_spec], k = args.k)

        X = filtered_spec_kmerised

        nrc = NearestCentroid(metric = args.d)
        nrc.fit(X, y)
        result_model = {"model_type": "Nearest-centroid", "model":nrc, "additional_data": [kmer]}
        with open(args.outfile, "wb") as f:
            dump(result_model, f)

    elif args.type == "cnb":
        remap_base = {"N":0, "A":1, "C":2, "G":3, "T":4, "-":5}
        filtered_spec_remap = [[remap_base[base] for base in str(a.seq).upper()] for a in filtered_spec]
        X = filtered_spec_remap
        cnb = CategoricalNB(alpha = args.alpha, fit_prior = args.fit_prior, min_categories = 6)
        cnb.fit(X, y)
        result_model = {"model_type": "CategoricalNB", "model":cnb, "additional_data": [remap_base]}
        with open(args.outfile, "wb") as f:
            dump(result_model, f)

