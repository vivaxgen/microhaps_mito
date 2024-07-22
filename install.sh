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

#echo "source \${VVG_BASEDIR}/env/G6PD-pipeline/activate.sh" >> ${BASEDIR}/bin/activate.sh
ln -sr ${ENVS_DIR}/microhaps_mito/etc/bashrc.d/50-microhaps-mito-pipeline ${BASHRC_DIR}/  # **UPDATE** the pipeline name (if changed)

git clone https://github.com/trmznt/seqpy.git ${ENVS_DIR}/seqpy
ln -sr ${ENVS_DIR}/seqpy/bin/spcli ${BASEDIR}/bin/spcli

mkdir -p ${ENVS_DIR}/seqpy/etc/bashrc.d

cat > ${ENVS_DIR}/seqpy/etc/bashrc.d/20-seqpy << 'EOF'
_script="$(readlink -f ${BASH_SOURCE[0]})"

# Delete last component from $_script
_mydir="$(dirname $_script)"

export PYTHONPATH=${PYTHONPATH}:"$(dirname $(dirname $_mydir))"
EOF

ln -sr ${ENVS_DIR}/seqpy/etc/bashrc.d/20-seqpy ${BASHRC_DIR}/ 

## Wait for update
curl https://raw.githubusercontent.com/trmznt/pys/master/seq/primersim.py > ${BASEDIR}/opt/primersim.py

pip3 install numba
micromamba install -y bioconda:muscle==5.1

echo "Reloading source files"
reload_vvg_profiles

echo "Creating model"
mkdir -p ${ENVS_DIR}/microhaps_mito/refs/models
#ngs-pl create-model --type nc -o ${ENVS_DIR}/microhaps_mito/refs/models/nc.pickle -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/ORIGINAL_Spec2-1638-1893.fasta
#ngs-pl create-model --type cnb -o ${ENVS_DIR}/microhaps_mito/refs/models/cnb.pickle -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/ORIGINAL_Spec2-1638-1893.fasta
#ngs-pl create-model --type pa -o ${ENVS_DIR}/microhaps_mito/refs/models/pa.pickle -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/ORIGINAL_Spec2-1638-1893.fasta
#ngs-pl create-model --type cons -o ${ENVS_DIR}/microhaps_mito/refs/consensus.fasta -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/ORIGINAL_Spec2-1638-1893.fasta
#ngs-pl create-model --type ensemble -o ${ENVS_DIR}/microhaps_mito/refs/models/ensemble.fasta -b ${ENVS_DIR}/microhaps_mito/refs/Mito_sync_cons05.bed ${ENVS_DIR}/microhaps_mito/refs/ORIGINAL_Spec2-1638-1893.fasta

ngs-pl create-model --type nc -o ${ENVS_DIR}/microhaps_mito/refs/models/nc.pickle --cascade -m ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product_muscle_spec.txt ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product.muscle.fasta
ngs-pl create-model --type cnb -o ${ENVS_DIR}/microhaps_mito/refs/models/cnb.pickle --cascade -m ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product_muscle_spec.txt ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product.muscle.fasta
ngs-pl create-model --type pa -o ${ENVS_DIR}/microhaps_mito/refs/models/pa.pickle --cascade -m ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product_muscle_spec.txt ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product.muscle.fasta
ngs-pl create-model --type cons -o ${ENVS_DIR}/microhaps_mito/refs/consensus.fasta --cascade ${ENVS_DIR}/microhaps_mito/refs/cons.fasta ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product.muscle.fasta
ngs-pl create-model --type ensemble -o ${ENVS_DIR}/microhaps_mito/refs/models/ensemble.pickle --cascade -m ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product_muscle_spec.txt ${ENVS_DIR}/microhaps_mito/refs/mit_filtered_product.muscle.fasta

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