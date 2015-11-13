#!/bin/bash 

module load mpi

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/euryale4/scratch/momm/NEMO-FABM-AMM7-3.6/AMM7-MYSRC
export XIOS_HOME=$HOME/git/XIOS-1.0-svn
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=GCC_PMPC
#ARCH=GCC_PMPC_DEBUG

cd $NEMO_BUILD_DIR
echo "Building NEMO-FABM..."

./makenemo -m $ARCH -n AMM7_TEST clean_config
./makenemo -m $ARCH -r AMM7 -n AMM7_TEST && rsync -a $NEMO_BUILD_DIR/AMM7_TEST/BLD/bin/nemo.exe $RUNDIR/
