from ngs_pipeline import cerr, cexit, get_snakefile_path, check_NGSENV_BASEDIR, snakeutils
from ngs_pipeline.cmds import run_snakefile
from pickle import load
from Bio import AlignIO

def init_argparser():
    p = snakeutils.init_argparser('Infer species based on mitochondrial DNA - per read')
    p.add_argument('--model', help='model to use for inference')
    p.add_argument('infile', nargs=1, help='read file')
    return p


def main(args):
    result = []
    model = load(open(args.model, 'rb'))
    sequences = AlignIO.read(open(args.infile), "fasta")
    is_cascade = model.get("is_cascade", False)

    result.append(model[0].get_full_result(sequences))

    if is_cascade:
        pass