#!/bin/bash 

module load intel intel-mpi netcdf-intelmpi hdf5-intelmpi

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=~/build/NEMO-shelf
export XIOS_HOME=$HOME/local/xios-intel
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=CETO_INTEL_NOSIGNEDZERO

cd $NEMO_BUILD_DIR
echo "Cleaning old build..."
rm -f $RUNDIR/nemo.exe
./makenemo -m $ARCH -n AMM7_BLD_SCRATCH clean_config

echo "Building NEMO-FABM..."
./makenemo -m $ARCH -r AMM7_FABM_ERSEM -n AMM7_BLD_SCRATCH > $RUNDIR/compile.log && rsync -a $NEMO_BUILD_DIR/AMM7_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/ && echo "Done."
