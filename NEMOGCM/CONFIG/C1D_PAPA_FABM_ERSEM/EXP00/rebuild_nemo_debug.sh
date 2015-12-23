#!/bin/bash

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/C1D_PAPA

export XIOS_HOME=$HOME/local/xios-gnu
cd $NEMO_BUILD_DIR
echo "Building NEMO..."
./makenemo -m "GCC_UBUNTU_DEBUG" -n C1D_PAPA_TEST && rsync -a $NEMO_BUILD_DIR/C1D_PAPA_TEST/BLD/bin/nemo.exe $RUNDIR/nemo-debug.exe
echo "Done."

