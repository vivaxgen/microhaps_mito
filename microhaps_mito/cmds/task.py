import os
import pathlib
from ngs_pipeline import cerr, cexit, get_snakefile_path, check_NGSENV_BASEDIR
from ngs_pipeline.cmds import run_snakefile


def init_argparser():
    p = run_snakefile.init_argparser()
    return p


def main(args):

    import pipeline

    # NGSENV_BASEDIR is the base directory of the current pipeline (G6PD)
    # NGSENV_BASEDIR = pathlib.Path(check_NGSENV_BASEDIR())
    # smk_basepath = NGSENV_BASEDIR / 'pipeline' / 'rules'

    # args.snakefile = smk_basepath / 'index_reference.smk'
    args.snakefile = get_snakefile_path('sometask.smk',
                                        from_module=pipeline)
    args.no_config_cascade = True
    args.force = True

    run_snakefile.main(args)


# EOF