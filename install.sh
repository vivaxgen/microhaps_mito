#!/usr/bin/bash

# installation script for NGS-PL extension

# optional variable:
# - BASEDIR
# - OMIT

set -eu

# run the base.sh
# Detect the shell from which the script was called
parent=$(ps -o comm $PPID |tail -1)
parent=${parent#-}  # remove the leading dash that login shells have
case "$parent" in
  # shells supported by `micromamba shell init`
  bash|fish|xonsh|zsh)
    shell=$parent
    ;;
  *)
    # use the login shell (basename of $SHELL) as a fallback
    shell=${SHELL##*/}
    ;;
esac

# Parsing arguments
if [ -t 0 ] && [ -z "${BASEDIR:-}" ]; then
  printf "Pipeline base directory? [./microhaps_mito] " # **UPDATE** the pipeline name (if changed)
  read BASEDIR
fi

# default value
BASEDIR="${BASEDIR:-./microhaps_mito}" # **UPDATE** the pipeline name (if changed)

uMAMBA_ENVNAME='microhaps_mito' # **UPDATE** the pipeline name (if changed)
OMIT='GATK4'
source <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

# **UPDATE** the pipeline dependencies (if changed)
#echo Installing apptainer
#micromamba -y install apptainer -c conda-forge -c bioconda
#micromamba -y install squashfuse -c conda-forge
pip3 install scikit-learn 


echo "Cloning Microhaps-mito pipeline"
git clone https://github.com/vivaxgen/microhaps_mito.git ${ENVS_DIR}/microhaps_mito # **UPDATE** the pipeline name (if changed)

ln -sr ${ENVS_DIR}/microhaps_mito/bin/update-pipeline.sh ${BASEDIR}/bin/update-pipeline.sh # **UPDATE** the pipeline name (if changed)
chmod +x ${BASEDIR}/bin/update-pipeline.sh

echo "Creating model"
ngs-pl create-model --type nc -o ${ENVS_DIR}/microhaps_mito/refs/models/nc.pickle -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.fasta
ngs-pl create-model --type cnb-o ${ENVS_DIR}/microhaps_mito/refs/models/cnb.pickle -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.fasta


#echo "source \${VVG_BASEDIR}/env/G6PD-pipeline/activate.sh" >> ${BASEDIR}/bin/activate.sh
ln -sr ${ENVS_DIR}/microhaps_mito/etc/bashrc.d/50-pipeline ${BASHRC_DIR}/  # **UPDATE** the pipeline name (if changed)

echo "Reloading source files"
reload_vvg_profiles


echo
echo "microhaps_mito has been successfully installed." # **UPDATE** the pipeline name (if changed)
echo "To activate the microhaps_mito environment, either run the activation script" # **UPDATE** the pipeline name (if changed)
echo "to get a new shell:"
echo
echo "    `realpath ${BINDIR}/activate`"
echo
echo "or source the activation script (eg. inside another script):"
echo
echo "    source `realpath ${BINDIR}/activate`"
echo

# EOF