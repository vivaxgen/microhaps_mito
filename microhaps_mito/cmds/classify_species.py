import os
import pathlib
from ngs_pipeline import cerr, cexit, get_snakefile_path, check_NGSENV_BASEDIR, snakeutils
from ngs_pipeline.cmds import run_snakefile


def init_argparser():
    p = snakeutils.init_argparser('Microhaps Mito pipeline task runner')
    p.add_argument('-o', '--outdir', required=True,
                   help='output directory')
    p.add_argument('-u', '--underscore', default=4, type=int,
                   help='number of undercore character to be stripped, '
                   'counted in reverse')
    p.add_argument('--model', action='append', default = ['*'], choices = ['*', 'cnb', 'nc', 'pa', 'ensemble'],
                    help='model to run (default: *)')
    p.add_argument('--strictness', action='append', default = ['conservative'], choices = ['conservative', 'sensitive'],
                    help='strictness of the filtering (default: conservative)')
    p.add_argument('infiles', nargs="*", help='read files')
    return p


def main(args):

    import microhaps_mito

    # NGSENV_BASEDIR is the base directory of the current pipeline (G6PD)
    # NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    # smk_basepath = NGSENV_BASEDIR / 'pipeline' / 'rules'

    if "*" in args.model:
        args.model = ['cnb', 'nc', 'pa', 'ensemble']

    # args.snakefile = smk_basepath / 'index_reference.smk'
    config=dict(
        outdir=args.outdir,
        underscore=args.underscore,
        infiles=args.infiles,
        model=args.model,
        strictness=args.strictness,
    )

    args.snakefile = get_snakefile_path('classify_species.smk',
                                        from_module=microhaps_mito)

    args.no_config_cascade = True
    args.force = True

    status, elapsed_time = snakeutils.run_snakefile(args, config=config)

    if not status:
        cerr('[WARNING: targeted variant calling did not successfully complete]')
    cerr(f'[Finish targeted variant calling (time: {elapsed_time})]')

# EOF