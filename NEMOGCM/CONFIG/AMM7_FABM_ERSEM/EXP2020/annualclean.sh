#!/bin/bash

# the run dir is now passed as argument from the pbs script
#RUNDIR=/work/n01/n01/yuti/AMM7-v1-hindcast
RUNDIR=$1

echo "Cleaning year " $RUNDIR "..."
rm -f $RUNDIR/bdy/amm7_bdy?_*m??d??.nc $RUNDIR/bdy/amm7_bt_bdy?_*m??d??.nc  $RUNDIR/bdy/amm7skag_bt_bdy?_*m??d??.nc  $RUNDIR/bdy/amm7skag_bdy?_*m??d??.nc $RUNDIR/fluxes/CUT_ERAI_INCLUDE_MSLP_y????m??d??.nc
rm -f $RUNDIR/AMM7-EMEP-NDeposition_y????.nc
echo "Done."

