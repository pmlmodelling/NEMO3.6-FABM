#!/bin/bash 

module load intel intel-mpi netcdf-intelmpi hdf5-intelmpi

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/work/n01/n01/momme/AMM7
export XIOS_HOME=$HOME/local/xios-intel-debug
export FABM_HOME=$HOME/local/fabm/nemo-debug

ARCH=CETO_INTEL_NOSIGNEDZERO_DEBUG

cd $NEMO_BUILD_DIR
echo "Cleaning old build ..."
rm -f $RUNDIR/nemo-debug.exe $NEMO_BUILD_DIR/AMM7_DBG_SCRATCH/BLD/bin/nemo.exe
./makenemo -m $ARCH -n AMM7_DBG_SCRATCH clean_config

echo "Building NEMO-FABM..."
./makenemo -m $ARCH -r AMM7_FABM_ERSEM -n AMM7_DBG_SCRATCH | tee $RUNDIR/compile-debug.log && mv $NEMO_BUILD_DIR/AMM7_DBG_SCRATCH/BLD/bin/nemo.exe $RUNDIR/nemo-debug.exe && echo "Done."
