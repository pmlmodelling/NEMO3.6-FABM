#!/bin/bash

set -u

# since /nerc is not available anymore, this script will compress the outputs file and "archive" them in a folder  $RUNDIR/ARCHIVE/yyyy/mm/
# tagging has been removed until a revised script is generated

MERGETOOL=$HOME/local/bin/nocscombine-5

folder=$1
ARCHIVEDIR=./ARCHIVE/

module load cray-netcdf

if [ -d $folder ]
then
  cd $folder
  echo "Archiving $folder:"
  mkdir -p $ARCHIVEDIR/$folder
  $MERGETOOL -f restart_0000.nc # && nccopy -k 4 -d 9 restart.nc $ARCHIVEDIR/$folder/restart.nc && echo "   restart.nc" || echo '   restart.nc failed!!!'
  $MERGETOOL -f restart_trc_0000.nc # && nccopy -k 4 -d 9 restart_trc.nc $ARCHIVEDIR/$folder/restart_trc.nc && echo "   restart_trc.nc" || echo '   restart_trc.nc failed!!!' 
  rsync -a ocean.output.bz2 $ARCHIVEDIR/$folder  
  for file in amm7*.nc
  do
     nccopy -k 4 -d 9 $file $ARCHIVEDIR/$folder/$file && rm $file && echo "   " $file || echo "   " $file 'failed!!!'
  done
  chmod -R g+rX $ARCHIVEDIR/$folder/*
  module load anaconda
  cd $ARCHIVEDIR
  #for file in $folder/*.nc
  #do
  #  echo "Tagging" $file "..." 
  #  python /work/n01/n01/momme/AMM7-v1/ncAddPMLLicenseMomme.py $file
  #  chmod g+rX $folder
  #  chmod -R g+rX $folder/..
  #done
fi
echo "Done."
