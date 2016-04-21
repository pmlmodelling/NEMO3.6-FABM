#!/bin/bash 

module unload PrgEnv-cray PrgEnv-gnu
module load PrgEnv-intel
module unload cray-netcdf
module load cray-netcdf-hdf5parallel

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/work/n01/n01/momme/AMM7
export XIOS_HOME=$HOME/local/xios-intel
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=XC_ARCHER_INTEL_NOSIGNEDZERO

cd $NEMO_BUILD_DIR
echo "Cleaning old build..."
rm -f $RUNDIR/nemo.exe
./makenemo -m $ARCH -n AMM7_BLD_SCRATCH clean_config

echo "Building NEMO-FABM..."
./makenemo -m $ARCH -r AMM7_FABM_ERSEM -n AMM7_BLD_SCRATCH > $RUNDIR/compile.log && rsync -a $NEMO_BUILD_DIR/AMM7_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/ && echo "Done."
