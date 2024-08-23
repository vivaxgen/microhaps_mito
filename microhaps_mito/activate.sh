# activation script
# this will run the pipeline using batch job submission system

# edit the following variables as necessary
export OPT=${VVG_BASEDIR}/opt # **UPDATE** the pipeline name (if changed)

PROMPT=uhaps_mito # **UPDATE** the pipeline name (if changed)
#NGS_PIPELINE_BASH_SOURCE=${VVG_BASEDIR}/env/ngs-pipeline/bin/activate.sh

# only edit the following lines if you do know what you are doing

_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

export NGSENV_BASEDIR=${_mydir}


# snakemake can handle sbatch
# export JOBCMD="sbatch --cpus-per-task={threads} --job-name=smk-{rule} --output=${NGSENV_BASEDIR}/slurm-logs/j-%j.out"

#. ${NGS_PIPELINE_BASH_SOURCE}
PATH=${NGSENV_BASEDIR}/bin:${PATH}
PYTHONPATH=${NGSENV_BASEDIR}:${PYTHONPATH}
export NGS_PIPELINE_CMD_MODS=microhaps_mito.cmds:${NGS_PIPELINE_CMD_MODS}

PS1="(${PROMPT}) [\u@\h \W]\$ "

export HISTFILE=${HOME}/.bash_history.${PROMPT}
history -c; history -r


# EOF
