#!/bin/bash 

module unload PrgEnv-cray PrgEnv-gnu
module load PrgEnv-intel
module load cray-netcdf-hdf5parallel

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/work/n01/n01/momme/AMM7
export XIOS_HOME=$HOME/local/xios-intel
export FABM_HOME=$HOME/local/fabm/nemo

ARCH=XC_ARCHER_INTEL_NOSIGNEDZERO
#ARCH=XC_ARCHER_INTEL_NOSIGNEDZERO_DEBUG

cd $NEMO_BUILD_DIR
echo "Building NEMO-FABM..."

./makenemo -m $ARCH -n AMM7_TEST clean_config
./makenemo -m $ARCH -r AMM7 -n AMM7_TEST && rsync -a $NEMO_BUILD_DIR/AMM7_TEST/BLD/bin/nemo.exe $RUNDIR/
