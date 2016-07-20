#!/bin/bash

RUNDIR=/work/n01/n01/momme/AMM7-v0

echo "Cleaning" $RUNDIR "..."
rm -f $RUNDIR/amm7_?d_*_grid_*.nc $RUNDIR/amm7_1m_*_grid*.nc $RUNDIR/amm7_1y_*_grid*.nc
rm -f $RUNDIR/amm7_?d_*_ptrc_*.nc $RUNDIR/amm7_1m_*_ptrc*.nc $RUNDIR/amm7_1y_*_ptrc*.nc
rm -rf $RUNDIR/restart.nc $RUNDIR/restart_trc.nc
echo "Done."
