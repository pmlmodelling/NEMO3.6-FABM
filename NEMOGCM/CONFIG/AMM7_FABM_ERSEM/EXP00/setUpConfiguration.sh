#!/bin/bash

INPUTPATH="/work/n01/n01/momme/AMM7-INPUTS"
RUNPATH="/work/n01/n01/momme/AMM7"
CONFIGPATH="$HOME"/git/NEMO-shelf/NEMOGCM/CONFIG/AMM7_FABM_ERSEM

#prepare folders for annual atmospheric and lateral forcings:

mkdir -p "$RUNPATH"/fluxes
rm -rf "$RUNPATH"/fluxes/flx_y????.nc
rm -rf "$RUNPATH"/fluxes/met_y????.nc
rm -rf "$RUNPATH"/fluxes/strd_y????.nc
rm -rf "$RUNPATH"/fluxes/ssrd24_y????.nc

mkdir -p "$RUNPATH"/bdy
rm -rf "$RUNPATH"/bdy/amm7*_y????m??d??.nc

#namelists:
#rsync -a "$CONFIGPATH"/EXP00/namelist*_cfg.template $RUNPATH/
#rsync -a "$CONFIGPATH"/EXP00/namelist*_ref $RUNPATH/
#rsync -a "$CONFIGPATH"/EXP00/namelist*_ref $RUNPATH/
#rsync -a "$CONFIGPATH"/EXP00/iodef.xml $RUNPATH/
#rsync -a "$CONFIGPATH"/EXP00/domain_def.xml $RUNPATH/
#rsync -a "$CONFIGPATH"/EXP00/field_def.xml $RUNPATH/

#atmospheric deposition:

ln -sf "$INPUTPATH"/AtmosphericDeposition/AMM7_atn_weights.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/AtmosphericDeposition/EMEP_atmosnuts_1980_2009.nc "$RUNPATH"/
#
#grid:
#ln -sf "$INPUTPATH"/GRID/coordinates.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/GRID/AMM7_coordinates.nc "$RUNPATH"/coordinates.nc
ln -sf "$INPUTPATH"/GRID/coordinates.bdy.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/GRID/coordinates.skagbdy.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/GRID/bathy_meter.nc "$RUNPATH"/

#Light:
ln -sf "$INPUTPATH"/Light/AMM7-ADY-broadband.nc "$RUNPATH/"
ln -sf "$INPUTPATH"/Light/kd490.nc "$RUNPATH/"

#interpolation weigths for surface flux forcing:
ln -sf "$INPUTPATH"/FLUXES/weights_erai_amm7_bicubic.nc "$RUNPATH"/fluxes/
ln -sf "$INPUTPATH"/FLUXES/weights_erai_amm7_bilin.nc "$RUNPATH"/fluxes/
ln -sf "$INPUTPATH"/FLUXES/CUT_ERAI_LSM.nc "$RUNPATH"/fluxes/

#rivers:
ln -sf "$INPUTPATH"/RIVERS/rivers.365.nc "$RUNPATH"/rivers.nc

#lateral boundaries from climatology:
ln -sf "$INPUTPATH"/BDY/amm7skagbdy_trc.nc "$RUNPATH"/
#ln -sf "$INPUTPATH"/BDY/amm7bdy_trc.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY/amm7bdy_trc.fixed.nc "$RUNPATH"/amm7bdy_trc.nc
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_K1_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_K1_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_K1_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_K2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_K2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_K2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_L2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_L2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_L2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_M2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_M2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_M2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_M4_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_M4_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_M4_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_MU2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_MU2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_MU2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_2N2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_2N2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_2N2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_N2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_N2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_N2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_NU2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_NU2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_NU2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_O1_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_O1_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_O1_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_P1_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_P1_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_P1_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_Q1_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_Q1_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_Q1_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_S1_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_S1_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_S1_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_S2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_S2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_S2_grid_V.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_T2_grid_T.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_T2_grid_U.nc "$RUNPATH"/bdy/
ln -sf "$INPUTPATH"/BDY/amm7_bdytide_T2_grid_V.nc "$RUNPATH"/bdy/

#restart files:
ln -sf "$INPUTPATH"/RESTARTS/restart.nc restart.nc
