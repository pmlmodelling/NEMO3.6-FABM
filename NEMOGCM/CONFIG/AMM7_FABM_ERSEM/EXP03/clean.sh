#!/bin/bash

RUNDIR=/work/momm/AMM7

echo "Cleaning" $RUNDIR "..."
rm -f $RUNDIR/bdy/amm7_bdy?_*m??d??.nc $RUNDIR/bdy/amm7_bt_bdy?_*m??d??.nc  $RUNDIR/bdy/amm7skag_bt_bdy?_*m??d??.nc  $RUNDIR/bdy/amm7skag_bdy?_*m??d??.nc $RUNDIR/fluxes/CUT_ERAI_INCLUDE_MSLP_y????m??d??.nc
rm -f $RUNDIR/AMM7-EMEP-NDeposition_y????.nc
rm -f $RUNDIR/amm7_?d_*_grid_*.nc $RUNDIR/amm7_1m_*_grid*.nc $RUNDIR/amm7_1y_*_grid*.nc
rm -f $RUNDIR/amm7_?d_*_ptrc_*.nc $RUNDIR/amm7_1m_*_ptrc*.nc $RUNDIR/amm7_1y_*_ptrc*.nc
echo "Done."
