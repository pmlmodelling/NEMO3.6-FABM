#!/bin/bash 

module -s restore /work/n01/shared/acc/n01_modules/ucx_env

work=/work/n01/n01/gle

NEMO_BUILD_DIR=$work/NEMO-ERSEM-shelf/NEMOGCM/CONFIG

MY_CONFIG=MY_AMM7_FABM_ERSEM

RUNDIR=`pwd`

export XIOS_HOME=$work/XIOS1
export FABM_HOME=$work/local/fabm/nemo-fabm-ersem

ARCH=X86_ARCHER2-Cray

cd $NEMO_BUILD_DIR

echo "Cleaning old build..."

rm -f $RUNDIR/nemo.exe $NEMO_BUILD_DIR/$MY_CONFIG/BLD/bin/nemo.exe
./makenemo -m $ARCH -n $MY_CONFIG clean_config

echo "Building NEMO-FABM..."

./makenemo -n $MY_CONFIG -r AMM7_FABM_ERSEM -m $ARCH -j 16 > $RUNDIR/compile.log && mv $NEMO_BUILD_DIR/$MY_CONFIG/BLD/bin/nemo.exe $RUNDIR/ && echo "Done."
