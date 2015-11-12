#!/bin/bash

INPUTPATH="/data/euryale4/to_archive/momm-AMM7-inputs"
RUNPATH="."
CONFIGPATH="$HOME"/git/NEMO-shelf/NEMOGCM/CONFIG/AMM7

#prepare folders for annual atmospheric and lateral forcings:

mkdir -p "$RUNPATH"/fluxes
rm -rf "$RUNPATH"/fluxes/flx_y????.nc
rm -rf "$RUNPATH"/fluxes/met_y????.nc
rm -rf "$RUNPATH"/fluxes/strd_y????.nc
rm -rf "$RUNPATH"/fluxes/ssrd24_y????.nc

mkdir -p "$RUNPATH"/bdy
rm -rf "$RUNPATH"/bdy/amm7*_y????m??d??.nc

#namelists:
rsync -a "$CONFIGPATH"/EXP00/namelist*_cfg.template $RUNPATH/
rsync -a "$CONFIGPATH"/EXP00/namelist*_ref $RUNPATH/
rsync -a "$CONFIGPATH"/EXP00/namelist*_ref $RUNPATH/
rsync -a "$CONFIGPATH"/EXP00/iodef.xml $RUNPATH/
rsync -a "$CONFIGPATH"/EXP00/domain_def.xml $RUNPATH/
rsync -a "$CONFIGPATH"/EXP00/field_def.xml $RUNPATH/

#atmospheric deposition:

ln -sf "$INPUTPATH"/AtmosphericDeposition/AMM7_atn_weights.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/AtmosphericDeposition/EMEP_atmosnuts_1980_2009.nc "$RUNPATH"/

#grid:
ln -sf "$INPUTPATH"/GRID/coordinates.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/GRID/coordinates.bdy.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/GRID/coordinates.skagbdy.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/GRID/bathy_meter.nc "$RUNPATH"/

#IOPs:
ln -sf "$INPUTPATH"/ADY/ady_data_netcdf3_AMM7_interp.nc ady_data_netcdf3.nc

#interpolation weigths for surface flux forcing:
ln -sf "$INPUTPATH"/FLUXES/ERA_INT_LSM/ERA_INT_AMM7_bilin.nc "$RUNPATH"/fluxes/

#rivers:
ln -sf "$INPUTPATH"/RIVERS/AMM_v36_rivers_trc_365.nc "$RUNPATH"/rivers.nc

#lateral boundaries from climatology:
ln -sf "$INPUTPATH"/BDY-clim/amm7skagbdy_TR1.nc "$RUNPATH"/amm7skagbdy_trc.nc
ln -sf "$INPUTPATH"/BDY-clim/amm7_bdyT_bio_WOA09_GLODAP_51.nc "$RUNPATH"/amm7bdy_trc.nc
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_K1_grid_T.nc "$RUNPATH"/

ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_K1_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_K1_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_K2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_K2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_K2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_L2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_L2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_L2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_M2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_M2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_M2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_M4_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_M4_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_M4_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_MU2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_MU2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_MU2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_2N2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_2N2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_2N2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_N2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_N2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_N2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_NU2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_NU2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_NU2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_O1_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_O1_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_O1_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_P1_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_P1_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_P1_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_Q1_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_Q1_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_Q1_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_S1_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_S1_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_S1_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_S2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_S2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_S2_grid_V.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_T2_grid_T.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_T2_grid_U.nc "$RUNPATH"/
ln -sf "$INPUTPATH"/BDY-clim/VN34TIDE/amm7_bdytide_T2_grid_V.nc "$RUNPATH"/

