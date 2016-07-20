#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=19
#SBATCH --ntasks-per-node=20
#SBATCH --threads-per-core=1
#SBATCH --job-name=AMM7-hind
#SBATCH --partition=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=momm
#SBATCH --exclusive


# NOTE:
# the simulation year and month must be passed to this pbs script at runtime
# by the command line, e.g.:
# "sbatch 1980 1 MonthlyChainHindcast.sh"
#

module load intel
module load intel-mpi
module load hdf5-intelmpi
module load netcdf-intelmpi

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_PIN_PROCS=0-19
export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_EXTRA_FILESYSTEM_LIST=gpfs
export OMP_NUM_THREADS=1
ulimit -s 10240
#ulimit -s unlimited
set -u #break on unset variables

ystart=1981
yend=2010
y0=$1
m0=$2

#Compute previous and next month:
mp=$(( $m0 + 1 ))
if [ $mp -eq 13 ]
then
   mp=01
   yp=$(( y0 + 1 ))
else
   yp=$y0
fi
mm=$(( $m0 - 1 ))
if [ $mm -eq 0 ]
then
   mm=12
   ym=$(( y0 - 1 ))
else
   ym=$y0
fi

m0str=$(printf %02d $m0)
mpstr=$(printf %02d $mp)
mmstr=$(printf %02d $mm)

echo "Previous month" $ym $mm "..."
echo "Launching" $y0 $m0 "..."
echo "Next month" $yp $mp "..."

RUNDIR=$HOME/run/AMM7
INPUTS=$HOME/AMM7-INPUTS
ARCHIVEDIR=$RUNDIR/$y0/$m0str

cd $RUNDIR

#cleanup:
./clean.sh

#link annual forcing files:
./setAnnualLinks.sh $y0

#restarts:
rm -rf restart.nc restart_trc.nc restart_[0-9]???.nc restart_trc_[0-9]???.nc
if [ $y0 -eq $ystart ] && [ $m0 -eq 1 ]
then
  ln -sf $INPUTS/RESTARTS/restart.nc restart.nc
  ln -sf $INPUTS/RESTARTS/restart_trc.full.nc restart_trc.nc
  rst=0
  euler=0
else
  if [ -s $RUNDIR/$ym/$mmstr/restart_0000.nc ]
  then 
     ln -sf $RUNDIR/$ym/$mmstr/restart_????.nc .
  elif [ -s $RUNDIR/$ym/$mmstr/restart.nc ]
  then
     ln -sf $RUNDIR/$ym/$mmstr/restart.nc .
  fi
  if [ -s $RUNDIR/$ym/$mmstr/restart_trc_0000.nc ]
  then 
     ln -sf $RUNDIR/$ym/$mmstr/restart_trc_????.nc .
  elif [ -s $RUNDIR/$ym/$mmstr/restart_trc.nc ]
  then
     ln -sf $RUNDIR/$ym/$mmstr/restart_trc.nc .
  fi
  rst=2
  euler=0
fi

#compute run-time:
case $m0 in
   4|6|9|11) nit=8640 ;;
   2) if [ $(( y0 % 4 )) -ne 0 -o $(( y0 % 100)) -eq 0 -a $(( $y0 % 400 )) -ne 0 ]; then nit=8064; else nit=8352; fi ;;
   *) nit=8928 ;;
esac

#compute start iteration and end iteration:
dt=300 #time step
nsstart=$(date -d $ystart-01-01 +%s) #seconds since EPOCH for total simulation start
ns0=$(date -d $y0-${m0str}-01 +%s) #seconds since EPOCH for this chunk
n0=$(( ns0 - nsstart ))
n0=$(( n0 / dt + 1 ))
nend=$(( n0 + nit -1 ))
d0=$y0${m0str}01
if [ $y0 -lt 1990 ]
then
   lclim=.true.
else
   lclim=.false.
fi
cat namelist.template \
                | sed "s,__DATE0__,$d0,g" \
                | sed "s,__RST__,$rst,g" \
                | sed "s,__N0__,$n0,g" \
                | sed "s,__NEND__,$nend,g" \
                | sed "s,__LCLIM__,$lclim,g" \
                | sed "s,__EULER__,$euler,g" \
                > namelist_cfg
cat namelist_top.template \
                | sed "s,__RST__,$rst,g" \
                > namelist_top_cfg

XIOSCORES=1
XIOSBLOCKEDCORES=19
CORES=$(( SLURM_NTASKS ))
COMPUTECORES=$(( CORES - XIOSCORES - XIOSBLOCKEDCORES ))
COMPUTETAG=0-$(( COMPUTECORES - 1 ))
if [ $XIOSCORES = 1 ]
then
   XIOSTAG=$COMPUTECORES
else
   XIOSTAG=$COMPUTECORES-$(( COMPUTECORES + XIOSCORES - 1 ))
fi

cat > ./xios.conf <<EOL
$COMPUTETAG  ./nemo.exe
$XIOSTAG ./xios_server.exe
EOL

echo "Total number of nodes / cores used:" $SLURM_NNODES "/" $CORES
echo "Compute cores used: $COMPUTECORES"
echo "XIOS cores used:" $XIOSCORES
echo "Buffer cores blocked for XIOS:" $XIOSBLOCKEDCORES

#cat > xios.conf <<EOL
#0-359  ./nemo.exe
#360    ./xios_server.exe
#EOL

#srun -n $COMPUTECORES ./nemo.exe &
#srun -n 1 ./xios_server.exe &
#wait

echo "Launching NEMO at $(date +%s) seconds since 1970-01-01 00:00:00"

srun -n $(( COMPUTECORES + XIOSCORES )) --multi-prog ./xios.conf
#srun -n $CORES ./nemo.exe

echo "Finished NEMO at $(date +%s) seconds since 1970-01-01 00:00:00"

#prepare archive directory
echo "Archiving to" $ARCHIVEDIR "..."
mkdir -p $ARCHIVEDIR
#prepare restart files:
for file in amm7_*${nend}_restart_????.nc
do
  fn=${file: -7:4} #file number
  mv $file $ARCHIVEDIR/restart_$fn.nc
done
for file in amm7_*${nend}_restart_trc_????.nc
do
  fn=${file: -7:4} #file number
  mv $file $ARCHIVEDIR/restart_trc_$fn.nc
done

#move outputs:
mv amm7_1d_${d0}_[1-2]???????_*.nc $ARCHIVEDIR
mv amm7_1m_${d0}_[1-2]???????_ptrc_T.nc $ARCHIVEDIR
mv amm7_1m_${d0}_[1-2]???????_grid_?.nc $ARCHIVEDIR
bzip2 ocean.output && mv -f ocean.output.bz2 $ARCHIVEDIR

#resubmit:
if [ $yp -le $yend ]
then
   echo "Submitting" $mp $yp '...'
   sbatch MonthlyChainHindcast.sh $yp $mp
   echo "Done."
else
   echo "All done."
fi

#archive:
./archiveFolder.sh $y0/$m0str >& ./archive.$y0.$m0str.log &
