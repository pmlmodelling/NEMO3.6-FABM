#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --threads-per-core=1
#SBATCH --job-name=ORCA2
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

#NOTE:
# start and end cycle of perpetual cycle need to be passed by command line argument as first and second argument respectively.

CORES=$(( SLURM_NTASKS ))

ncycle=$1
nend=$2
nstart=1
nm=$(( ncycle - 1 ))
np=$(( ncycle + 1 ))
set -u #break on unset variables

RUNDIR=/work/momm/ORCA2
ARCHIVEDIR=$RUNDIR/$ncycle

cd $RUNDIR

stopflag=0

#restarts:
rm -rf restart.nc restart_trc.nc restart_[0-9]???.nc restart_trc_[0-9]???.nc
if [ $ncycle -eq $nstart ]
then
  ln -sf $RUNDIR/0/restart_trc.nc restart_trc.nc
  lrst=.FALSE.
  rst=0
  euler=0
else
  lrst=.TRUE.
  if [ -s $RUNDIR/$nm/restart_0000.nc ]
  then
     ln -sf $RUNDIR/$nm/restart_????.nc .
  elif [ -s $RUNDIR/$nm/restart.nc ]
  then
     ln -sf $RUNDIR/$nm/restart.nc .
  else
     stopflag=1
  fi
  if [ -s $RUNDIR/$nm/restart_trc_0000.nc ]
  then
     ln -sf $RUNDIR/$nm/restart_trc_????.nc .
  elif [ -s $RUNDIR/$nm/restart_trc.nc ]
  then
     ln -sf $RUNDIR/$nm/restart_trc.nc .
  else
     stopflag=1
  fi
  rst=2
  euler=0
fi

#compute run-time:
dt=5760 #time step in seconds
nit=$(( 86400*365/dt ))

#compute start iteration and end iteration:
ns0=$(( 86400*365 * (ncycle - 1) ))
n0=$(( ns0 / dt + 1 ))
nend=$(( n0 + nit -1 ))
d0=$(printf %04d $ncycle)0101

#restart
cat namelist.template \
                | sed "s,__DATE0__,$d0,g" \
                | sed "s,__LRST__,$lrst,g" \
                | sed "s,__RST__,$rst,g" \
                | sed "s,__N0__,$n0,g" \
                | sed "s,__NEND__,$nend,g" \
                | sed "s,__EULER__,$euler,g" \
                > namelist_cfg
cat namelist_top.template \
                | sed "s,__RST__,$rst,g" \
                > namelist_top_cfg

#cat > ./xios.conf <<EOL
#$COMPUTETAG  ./nemo.exe
#$XIOSTAG ./xios_server.exe
#EOL

echo "Total number of nodes / cores used:" $SLURM_NNODES "/" $CORES
#echo "Compute cores used: $COMPUTECORES"
#echo "XIOS cores used:" $XIOSCORES
#echo "Buffer cores blocked for XIOS:" $XIOSBLOCKEDCORES

#cat > xios.conf <<EOL
#0-359  ./nemo.exe
#360    ./xios_server.exe
#EOL

if [ $stopflag -ne 0 ]
then
   echo "Not ready to launch. Forced exit."
   exit
fi

echo "Launching NEMO at $(date +%s) seconds since 1970-01-01 00:00:00"

#srun -K1 -n $(( COMPUTECORES + XIOSCORES )) -m plane=20 --multi-prog ./xios.conf
srun -K1 -n $CORES ./nemo.exe

echo "Finished NEMO at $(date +%s) seconds since 1970-01-01 00:00:00"

#prepare archive directory
echo "Archiving to" $ARCHIVEDIR "..."
mkdir -p $ARCHIVEDIR
#prepare restart files:
for file in ORCA2_*${nend}_restart_????.nc
do
  fn=${file: -7:4} #file number
  mv $file $ARCHIVEDIR/restart_$fn.nc
done
for file in ORCA2_*${nend}_restart_trc_????.nc
do
  fn=${file: -7:4} #file number
  mv $file $ARCHIVEDIR/restart_trc_$fn.nc
done
for file in ORCA2_*${nend}_restart_ice_????.nc
do
  fn=${file: -7:4} #file number
  mv $file $ARCHIVEDIR/restart_ice_$fn.nc
done

mv ORCA2_5d_*.nc $ARCHIVEDIR
mv ORCA2_1m_ptrc_T.nc $ARCHIVEDIR
mv ORCA2_1m_grid_?.nc $ARCHIVEDIR
bzip2 ocean.output && mv -f ocean.output.bz2 $ARCHIVEDIR

#resubmit:
if [ $np -le $nendcycle ]
then
   echo "Submitting cycle" $np 'of' $nendcycle '...'
   sbatch Perpetual.sh $np $nendcycle
   echo "Done."
else
   echo "All done."
fi

#archive:
./archiveFolder.sh $ncycle >& ./archive.$ncycle.log &

