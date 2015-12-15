#!/bin/bash
# Command line arguments:
#   1: simulation year
#

INPUTPATH="/work/n01/n01/momme/AMM7-inputs"
RUNPATH="/work/n01/n01/momme/AMM7"
y0=1995

yn=$1
yb=$(( $yn-1 ))
ya=$(( $yn+1 ))

#restarts:
rm -rf "$RUNPATH"/restart.nc
if [ $yn -eq $y0 ]
then
    ln -sf "$INPUTPATH"/RESTART/temp_ic_noriver.nc "$RUNPATH"/restart.nc
fi

#prepare folders for annual atmospheric and lateral forcings:
mkdir -p "$RUNPATH"/fluxes
rm -rf "$RUNPATH"/fluxes/flx_y????.nc
rm -rf "$RUNPATH"/fluxes/met_y????.nc
rm -rf "$RUNPATH"/fluxes/strd_y????.nc
rm -rf "$RUNPATH"/fluxes/ssrd24_y????.nc
mkdir -p "$RUNPATH"/bdy
rm -rf "$RUNPATH"/bdy/amm7*_y????m??d??.nc

for y in $yb $yn $ya
do
  #atmospheric forcings:
  ln -sf "$INPUTPATH"/FLUXES/ERA_INT_LSM/flx_y"$y".nc "$RUNPATH"/fluxes/
  ln -sf "$INPUTPATH"/FLUXES/ERA_INT_LSM/met_y"$y".nc "$RUNPATH"/fluxes/
  ln -sf "$INPUTPATH"/FLUXES/ERA_INT_LSM/strd_y"$y".nc "$RUNPATH"/fluxes/
  ln -sf "$INPUTPATH"/FLUXES/ERA_INT_LSM/ssrd24_y"$y".nc "$RUNPATH"/fluxes/

  #lateral boundary conditions:
  ln -sf "$INPUTPATH"/BDY/NWS_INPUT/amm7_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/NWS_INPUT/amm7_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/NWS_INPUT/amm7_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/NWS_INPUT/amm7_bt_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/NWS_INPUT/amm7_bt_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/NWS_INPUT/amm7_bt_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/amm7_skag/amm7skag_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/amm7_skag/amm7skag_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/amm7_skag/amm7skag_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/amm7_skag/amm7skag_bt_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/amm7_skag/amm7skag_bt_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/amm7_skag/amm7skag_bt_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy/
done
