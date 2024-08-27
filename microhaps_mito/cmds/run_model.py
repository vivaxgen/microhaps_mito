from ngs_pipeline import cerr, cexit, get_snakefile_path, check_NGSENV_BASEDIR, snakeutils
from ngs_pipeline.cmds import run_snakefile
from pickle import load
from Bio import AlignIO

def init_argparser():
    p = snakeutils.init_argparser('Infer species based on mitochondrial DNA - per read')
    p.add_argument('--model', help='model to use for inference')
    p.add_argument('--output', help='output filename')
    p.add_argument('--sample', help='sample name')
    p.add_argument('infile', nargs=1, help='read pair file')
    return p


def main(args):
    print(args.infile[0])
    print(args.model)
    model = load(open(args.model, 'rb'))
    species_classes = sum([(a.tolist()) for a in model.classes_.values()], [])

    sequences = None
    # import ipdb; ipdb.set_trace()
    try:
        sequences = AlignIO.read(open(args.infile[0]), "fasta")
    # No read
    except Exception as e:
        print(e)
        classes = "\t".join(species_classes)
        output_string = f"sample\tspecies\tread_id\tf_pair_diff\t{classes}\n"
        output_string += f"{args.sample}\tNA\tNA\tNA\t{"\t".join(["NA" for _ in species_classes])}\n"
        open(args.output, "w+").write(output_string)

    if not sequences is None:
        result_df = model.get_full_result(sequences)
        result_df["sample"] = args.sample
        result_df.to_csv(args.output, sep="\t", index = False)