#!/bin/bash

NPMAX=600 # Maximum of Processes to possibly run
DATADIR=/work/n01/n01/slwa/NEMO/data/
TIDEDIR=/work/n01/n01/slwa/NEMO/data/Tides/VN34TIDE/
GRIDDIR=$DATADIR/GRID_v3.4/
RIVERDIR=$DATADIR/Rivers/ 

NUTATMODIR=/work/n01/n01/slwa/NEMO/data/AtmosphericNutrients
RSTDIR=/work/n01/n01/slwa/NEMO/data/ICs/
ADYDIR=/work/n01/n01/slwa/NEMO/data/ADY/
BDYDIR=$DATADIR/BDY/NWS_INPUT/
BDYDIR2=$DATADIR/BDY/amm7_skag/
WDIRBDY=bdy/
ATMOSDIR=$DATADIR/ERA_INT_LSM/
INITDIR=/work/n01/n01/slwa/NEMO/data/ICs/INIT_TS

#Boundary and TIDAL FORCING
if [ ! -d $WDIRBDY ];then
mkdir $WDIRBDY
fi

ln -sf $BDYDIR/* $WDIRBDY/
ln -sf $BDYDIR2/* $WDIRBDY/


ln -sf $TIDEDIR/amm7_bdytide_* $WDIRBDY/



#GRID
ln -sf $GRIDDIR/bathy_meter.nc bathy_meter.nc
ln -sf $GRIDDIR/coordinates.nc coordinates.nc

ln -sf $GRIDDIR/coordinates.bdy.nc coordinates.bdy.nc
ln -sf $GRIDDIR/coordinates.skagbdy.nc coordinates.skagbdy.nc


#RIVER FORCING 
#ln -sf $RIVERDIR/amm7_rivers_biogeo_corrected_fillvals.nc AMM_rivers.nc
#ln -sf $RIVERDIR/AMM_v34_rivers_trc_365.nc AMM_rivers.nc
#ln -sf $RIVERDIR/AMM_v34_rivers_trc_365.nc river.nc

#ERSEM boundary
#ln -sf $DATADIR/ERSEM/amm7_bdyT_bio_WOA09_GLODAP_51.nc  amm7bdy_trc.nc

#ATMOSPHERIC FORCING
if [ ! -d fluxes ];then
mkdir fluxes
fi
ln -sf $ATMOSDIR/*nc fluxes/

#INITIAL CONDITIONS
#ln -sf $INITDIR/ORCA1_EXP550_AMM7_1980_inits.nc inits.nc

#GENERIC RESTART
#ln -sf $RSTDIR/restart-phys-temp.nc restart.nc
#ln -sf $RSTDIR/restart_202612.nc restart.nc
ln -sf /work/n01/n01/slwa/NEMO/data/RESTART/qwes00.amm7.restarts.20100101/restart_20100101.nc restart.nc


#ERSEM INPUTS
#  ADY
#ln -sf $ADYDIR/ady_data_netcdf3_Filled.nc ady_data_netcdf3.nc
#ln -sf $ADYDIR/ady_data_netcdf3_Filled.nc ady_data_netcdf3
#ln -sf $ADYDIR/AMM7_ady_weights.nc AMM7_ady_weights.nc
#ln -sf $ADYDIR/ady_data_netcdf3_AMM7.nc ady_data_netcdf3.nc
#ln -sf $ADYDIR/ady_data_netcdf3_AMM7.nc ady_data_netcdf3
ln -sf $ADYDIR/ady_data_netcdf3_AMM7_interp.nc ady_data_netcdf3.nc
ln -sf $ADYDIR/ady_data_netcdf3_AMM7_interp.nc ady_data_netcdf3

#  ATMOSPHERIC NUTRIENTS
ln -sf $NUTATMODIR/AMM7_atn_weights.nc AMM7_atn_weights.nc
ln -sf $NUTATMODIR/EMEP_atmosnuts_1980_2009.nc EMEP_atmosnuts_1980_2009.nc
#  RIVER LOADS
#ln -sf $RIVERDIR/amm7_rivers_biogeo_corrected_fillvals.nc AMM_rivers_trc.nc
ln -sf $RIVERDIR/AMM_v34_rivers_trc_365.nc AMM_rivers_trc.nc
#  OPEN BC CLIMATOLOGY
ln -sf $DATADIR/ERSEM/amm7_bdyT_bio_WOA09_GLODAP_51.nc  amm7bdy_trc.nc
ln -sf $DATADIR/BDY/amm7skagbdy_trc_sample.nc amm7skagbdy_trc.nc
#  ERSEM GENERIC RESTART
#ln -sf $RSTDIR/nemoersem_amm_7k_restart_19691231_ady_bioalk-clean_l2_woa.nc restart_trc.nc
ln -sf $RSTDIR/nemoersem_amm_7k_restart_SSB.nc  restart_trc.nc


ln -sf /work/n01/n01/slwa/NEMO/data/Rivers/AMM_v36_rivers_trc_365.nc rivers.nc
ln -sf /work/n01/n01/slwa/NEMO/data/RESTART/temp_ic_noriver.nc restart.nc
