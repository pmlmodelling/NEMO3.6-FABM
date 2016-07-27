#!/bin/bash
# Command line arguments:
#   1: simulation year
#

INPUTPATH="/work/n01/n01/momme/AMM7-INPUTS"
RUNPATH="/work/n01/n01/momme/AMM7"

yn=$1
yb=$(( $yn-1 ))
ya=$(( $yn+1 ))

#rivers:
if [ $(( yn % 4)) -ne 0 -o $(( yn % 100)) -eq 0 -a $(( yn % 400)) -ne 0 ]
then
    ln -sf "$INPUTPATH"/RIVERS/rivers.ersem.nc "$RUNPATH"/rivers.nc
else
    ln -sf "$INPUTPATH"/RIVERS/rivers.ersem.leap.nc "$RUNPATH"/rivers.nc
fi

#prepare folders for annual atmospheric and lateral forcings:
mkdir -p "$RUNPATH"/fluxes
rm -rf "$RUNPATH"/fluxes/CUT_ERAI_INCLUDE_MSLP_y????m??d??.nc
mkdir -p "$RUNPATH"/bdy
rm -rf "$RUNPATH"/bdy/amm7*_y????m??d??.nc

for y in $yb $yn $ya
do
  #Atmospheric deposition of nutrients:
  ln -sf "$INPUTPATH"/AtmosphericDeposition/AMM7-EMEP-NDeposition."$y".nc "$RUNPATH"/AMM7-EMEP-NDeposition_y"$y".nc
  #ln -sf "$INPUTPATH"/AtmosphericDeposition/AMM7-EMEP-NDeposition."$y".nc "$RUNPATH"/AMM7-EMEP-NDeposition.nc

  #atmospheric forcings:
  ln -sf "$INPUTPATH"/FLUXES/CUT_ERAI_INCLUDE_MSLP_y"$y"m??d??.nc "$RUNPATH"/fluxes/

  #lateral boundary conditions:
  ln -sf "$INPUTPATH"/BDY/amm7_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy
  #ln -sf "$INPUTPATH"/BDY/amm7_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy #only needed id full baroclinic velocity fields are required
  #ln -sf "$INPUTPATH"/BDY/amm7_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy
  ln -sf "$INPUTPATH"/BDY/amm7_bt_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/amm7_bt_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy/
  ln -sf "$INPUTPATH"/BDY/amm7_bt_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy/
  if [ $y -lt 1990 -o $y -gt 2009 ]
  then
     ln -sf "$INPUTPATH"/BDY/amm7skag_bdyT_m??d??.nc "$RUNPATH"/bdy
     #ln -sf "$INPUTPATH"/BDY/amm7skag_bdyU_m??d??.nc "$RUNPATH"/bdy #only needed if full baroclinic velocity fields are required
     #ln -sf "$INPUTPATH"/BDY/amm7skag_bdyV_m??d??.nc "$RUNPATH"/bdy
     ln -sf "$INPUTPATH"/BDY/amm7skag_bt_bdyT_m??d??.nc "$RUNPATH"/bdy/
     ln -sf "$INPUTPATH"/BDY/amm7skag_bt_bdyU_m??d??.nc "$RUNPATH"/bdy/
     ln -sf "$INPUTPATH"/BDY/amm7skag_bt_bdyV_m??d??.nc "$RUNPATH"/bdy/
  else
     ln -sf "$INPUTPATH"/BDY/amm7skag_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy
     #ln -sf "$INPUTPATH"/BDY/amm7skag_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy #only needed if full baroclinic velocity fields are required
     #ln -sf "$INPUTPATH"/BDY/amm7skag_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy
     ln -sf "$INPUTPATH"/BDY/amm7skag_bt_bdyT_y"$y"m??d??.nc "$RUNPATH"/bdy/
     ln -sf "$INPUTPATH"/BDY/amm7skag_bt_bdyU_y"$y"m??d??.nc "$RUNPATH"/bdy/
     ln -sf "$INPUTPATH"/BDY/amm7skag_bt_bdyV_y"$y"m??d??.nc "$RUNPATH"/bdy/
  fi
done
