import os
import pathlib
from ngs_pipeline import cerr, cexit, get_snakefile_path, check_NGSENV_BASEDIR, snakeutils
from ngs_pipeline.cmds import run_snakefile


def init_argparser():
    p = snakeutils.init_argparser('Microhaps Mito pipeline models builder')
    return p


def main(args):

    import microhaps_mito

    # NGSENV_BASEDIR is the base directory of the current pipeline (G6PD)
    # NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    # smk_basepath = NGSENV_BASEDIR / 'pipeline' / 'rules'

    # args.snakefile = smk_basepath / 'index_reference.smk'

    args.snakefile = get_snakefile_path('create_model.smk',
                                        from_module=microhaps_mito)

    args.no_config_cascade = True
    args.force = True

    status, elapsed_time = snakeutils.run_snakefile(args, config=dict())

    if not status:
        cerr('[WARNING: models building did not successfully complete]')
    cerr(f'[Finish models building (time: {elapsed_time})]')

# EOF