#!/bin/bash

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/C1D_PAPA
export XIOS_HOME=$HOME/local/xios-gnu
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=GCC_UBUNTU

cd $NEMO_BUILD_DIR
echo "Building NEMO..."

./makenemo -m $ARCH -n C1D_PAPA_TEST clean_config
./makenemo -m $ARCH -r C1D_PAPA_FABM_ERSEM -n C1D_PAPA_TEST && rsync -a $NEMO_BUILD_DIR/C1D_PAPA_TEST/BLD/bin/nemo.exe $RUNDIR/
echo "Done."
