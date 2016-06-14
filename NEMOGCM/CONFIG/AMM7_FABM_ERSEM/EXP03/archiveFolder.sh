#!/bin/bash

set -u

ARCHIVEDIR=/pmldata/euryale6/to_archive/ceto-AMM7-HINDCAST-v0
MERGETOOL=$HOME/local/bin/nocscombine4
TMPDIR=$HOME/tmp

folder=$1

module load intel
module load intel-mpi
module load hdf5-intelmpi
module load netcdf-intelmpi

if [ -d $folder ]
then
  cd $folder
  echo "Archiving $folder:"
  mkdir -p $TMPDIR/$folder
  rm -rf $TMPDIR/$folder/*
  for file in amm7*.nc
  do
     nccopy -k 4 -d 9 $file $TMPDIR/$folder/$file && rm $file && echo "   " $file || echo "   " $file 'failed!!!'
     #nccopy -k 4 -d 9 $file $TMPDIR/$folder/$file && echo "   " $file || echo "   " $file 'failed!!!'
  done
  $MERGETOOL -f restart_0000.nc && nccopy -k 4 -d 9 restart.nc $TMPDIR/$folder/restart.nc && rm restart_????.nc && echo "   restart.nc" || echo '   restart.nc failed!!!'
  $MERGETOOL -f restart_trc_0000.nc && nccopy -k 4 -d 9 restart_trc.nc $TMPDIR/$folder/restart_trc.nc && rm restart_trc_????.nc && echo "   restart_trc.nc" || echo '   restart_trc.nc failed!!!' 
  #$MERGETOOL -f restart_0000.nc && nccopy -k 4 -d 9 restart.nc $TMPDIR/$folder/restart.nc && echo "   restart.nc" || echo '   restart.nc failed!!!'
  #$MERGETOOL -f restart_trc_0000.nc && nccopy -k 4 -d 9 restart_trc.nc $TMPDIR/$folder/restart_trc.nc && echo "   restart_trc.nc" || echo '   restart_trc.nc failed!!!' 
  mv -f ocean.output.bz2 $TMPDIR/$folder
fi
echo "Syncing to $ARCHIVEDIR"
ssh ceto7 "mkdir -p $ARCHIVEDIR/"
# sync with -a without preserving group and ownser:
ssh ceto7 "rsync -rlptD --exclude=*.tmp --exclude=*.nc.* $TMPDIR/* $ARCHIVEDIR/" && rm -rf $TMPDIR/$folder
echo "Done."
