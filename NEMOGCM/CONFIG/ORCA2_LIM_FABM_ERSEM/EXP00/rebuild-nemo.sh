#!/bin/bash 

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/ORCA2
export XIOS_HOME=$HOME/local/xios-gnu
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=GCC_UBUNTU

cd $NEMO_BUILD_DIR
echo "Building NEMO-FABM..."
rm -f $RUNDIR/nemo.exe
./makenemo -m $ARCH -n ORCA2_LIM_FABM_BLD_SCRATCH > $RUNDIR/compile.log && rsync -a $NEMO_BUILD_DIR/ORCA2_LIM_FABM_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/ $$ echo "Done."
