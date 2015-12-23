#!/bin/bash

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/C1D_PAPA

export XIOS_HOME=$HOME/local/xios-gnu-debug
cd $NEMO_BUILD_DIR
echo "Building NEMO..."
./makenemo -m "GCC_UBUNTU_DEBUG" -n C1D_PAPA_TEST clean_config
#./makenemo -m "GCC_UBUNTU_DEBUG" -r C1D_PAPA_FABM -n C1D_PAPA_TEST && rsync -a $NEMO_BUILD_DIR/ORCA2_LIM_FABM/BLD/bin/nemo.exe $RUNDIR/
./makenemo -m "GCC_UBUNTU_DEBUG" -r C1D_PAPA_FABM_ERSEM -n C1D_PAPA_TEST && rsync -a $NEMO_BUILD_DIR/C1D_PAPA_TEST/BLD/bin/nemo.exe $RUNDIR/nemo-debug.exe
#./makenemo -m "GCC_UBUNTU" -n C1D_PAPA_TEST && rsync -a $NEMO_BUILD_DIR/C1D_PAPA_TEST/BLD/bin/nemo.exe $RUNDIR/
echo "Done."

