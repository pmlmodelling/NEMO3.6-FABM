#!/bin/bash 

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/ORCA2
export XIOS_HOME=$HOME/local/xios-gnu
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=GCC_UBUNTU

cd $NEMO_BUILD_DIR
echo "Building NEMO-FABM..."

./makenemo -m $ARCH -n ORCA2_LIM_FABM_TEST clean_config
./makenemo -m $ARCH -r ORCA2_LIM_FABM_ERSEM -n ORCA2_LIM_FABM_TEST && rsync -a $NEMO_BUILD_DIR/ORCA2_LIM_FABM_TEST/BLD/bin/nemo.exe $RUNDIR/
echo "Done."
