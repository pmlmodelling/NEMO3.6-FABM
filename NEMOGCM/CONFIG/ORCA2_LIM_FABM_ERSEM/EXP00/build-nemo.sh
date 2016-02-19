#!/bin/bash 

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/data/ORCA2
export XIOS_HOME=$HOME/local/xios-gnu
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=GCC_UBUNTU

cd $NEMO_BUILD_DIR
echo "Cleaning old build..."
rm -f $RUNDIR/nemo.exe
./makenemo -m $ARCH -n ORCA2_LIM_FABM_BLD_SCRATCH clean_config

echo "Building NEMO-FABM..."
./makenemo -m $ARCH -r ORCA2_LIM_FABM_ERSEM -n ORCA2_LIM_FABM_BLD_SCRATCH > $RUNDIR/compile.log && rsync -a $NEMO_BUILD_DIR/ORCA2_LIM_FABM_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/ && echo "Done."
