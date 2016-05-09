#!/bin/bash 

module unload PrgEnv-cray PrgEnv-gnu
module load PrgEnv-intel
module unload cray-netcdf
module load cray-netcdf-hdf5parallel

NEMO_BUILD_DIR=$HOME/git/NEMO-shelf/NEMOGCM/CONFIG
RUNDIR=/work/n01/n01/yuti/AMM7-v1
export XIOS_HOME=$HOME/local/xios-intel-debug
export FABM_HOME=$HOME/local/fabm/nemo-debug

ARCH=XC_ARCHER_INTEL_NOSIGNEDZERO_DEBUG

cd $NEMO_BUILD_DIR
echo "Cleaning old build ..."
rm -f $RUNDIR/nemo-debug.exe
./makenemo -m $ARCH -n AMM7_BLD_SCRATCH clean_config

echo "Building NEMO-FABM..."
./makenemo -m $ARCH -r AMM7_FABM_ERSEM -n AMM7_BLD_SCRATCH | tee $RUNDIR/compile-debug.log && rsync -a $NEMO_BUILD_DIR/AMM7_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/nemo-debug.exe && echo "Done."
