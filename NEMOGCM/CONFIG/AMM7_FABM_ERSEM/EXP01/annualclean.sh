#!/bin/bash

RUNDIR=/work/n01/n01/momme/AMM7-v0

echo "Cleaning" $RUNDIR "..."
rm -f $RUNDIR/bdy/amm7_bdy?_*m??d??.nc $RUNDIR/bdy/amm7_bt_bdy?_*m??d??.nc  $RUNDIR/bdy/amm7skag_bt_bdy?_*m??d??.nc  $RUNDIR/bdy/amm7skag_bdy?_*m??d??.nc $RUNDIR/fluxes/CUT_ERAI_INCLUDE_MSLP_y????m??d??.nc
rm -f $RUNDIR/AMM7-EMEP-NDeposition_y????.nc
echo "Done."
