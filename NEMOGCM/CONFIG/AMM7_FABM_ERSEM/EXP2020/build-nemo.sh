#!/bin/bash 

module unload PrgEnv-cray PrgEnv-gnu
module load PrgEnv-intel
module unload cray-netcdf
module load cray-netcdf-hdf5parallel

NEMO_BUILD_DIR=<put here the path to the CONFIG folder, e.g. .../NEMO-shelf/NEMOGCM/CONFIG>
RUNDIR=<put here the path to the running folder>
export XIOS_HOME=<puth here the path to the compiled XIOS, e.g. $HOME/local/xios-intel>
export FABM_HOME=<puth here the path to the compiled FABM, e.g. $HOME/local/fabm/nemo>

#ARCH is the name of the Architecture file to be used in compilation
ARCH=XC_ARCHER_INTEL_NOSIGNEDZERO

cd $NEMO_BUILD_DIR
echo "Cleaning old build..."
rm -f $RUNDIR/nemo.exe $NEMO_BUILD_DIR/AMM7_BLD_SCRATCH/BLD/bin/nemo.exe
./makenemo -m $ARCH -n AMM7_BLD_SCRATCH clean_config

echo "Building NEMO-FABM..."
./makenemo -m $ARCH -r AMM7_FABM_ERSEM -n AMM7_BLD_SCRATCH > $RUNDIR/compile.log && mv $NEMO_BUILD_DIR/AMM7_BLD_SCRATCH/BLD/bin/nemo.exe $RUNDIR/ && echo "Done."
