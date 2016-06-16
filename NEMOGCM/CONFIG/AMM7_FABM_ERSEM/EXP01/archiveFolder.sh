#!/bin/bash

set -u

ARCHIVEDIR=/nerc/n01/n01/momme/AMM7-HINDCAST-15.06
MERGETOOL=$HOME/local/bin/nocscombine-5

folder=$1

module load cray-netcdf

if [ -d $folder ]
then
  cd $folder
  echo "Archiving $folder:"
  mkdir -p $ARCHIVEDIR/$folder
  for file in amm7*.nc
  do
     nccopy -k 4 -d 9 $file $ARCHIVEDIR/$folder/$file && rm $file && echo "   " $file || echo "   " $file 'failed!!!'
  done
  $MERGETOOL -f restart_0000.nc && nccopy -k 4 -d 9 restart.nc $ARCHIVEDIR/$folder/restart.nc && rm restart_????.nc && echo "   restart.nc" || echo '   restart.nc failed!!!'
  $MERGETOOL -f restart_trc_0000.nc && nccopy -k 4 -d 9 restart_trc.nc $ARCHIVEDIR/$folder/restart_trc.nc && rm restart_trc_????.nc && echo "   restart_trc.nc" || echo '   restart_trc.nc failed!!!' 
fi
echo "Done."
