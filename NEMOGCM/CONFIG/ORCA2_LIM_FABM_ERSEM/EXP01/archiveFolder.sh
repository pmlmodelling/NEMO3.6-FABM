set -u

folder=$1
ARCHIVEDIR=/pmldata/euryale1/momm/ORCA2
MERGETOOL=$HOME/local/bin/nocscombine4
TMPDIR=/work/momm/tmp

if [ -d $folder ]
then
  cd $folder
  echo "Archiving $folder:"
  mkdir -p $TMPDIR/$folder
  rm -rf $TMPDIR/$folder/*
  $MERGETOOL -f restart_0000.nc && nccopy -k 4 -d 9 restart.nc $TMPDIR/$folder/restart.nc && echo "   restart.nc" || echo '   restart.nc failed!!!'
  $MERGETOOL -f restart_trc_0000.nc && nccopy -k 4 -d 9 restart_trc.nc $TMPDIR/$folder/restart_trc.nc && echo "   restart_trc.nc" || echo '   restart_trc.nc failed!!!' 
  for file in ORCA2_1*.nc
  do
     nccopy -k 4 -d 9 $file $TMPDIR/$folder/$file && rm $file && echo "   " $file || echo "   " $file 'failed!!!'
  done
  rsync -a ocean.output.bz2 $TMPDIR/$folder
  chmod g+rX $TMPDIR/$folder/*
fi
echo "Syncing to $ARCHIVEDIR"
ssh ceto7 "mkdir -p $ARCHIVEDIR/"
# sync with -a without preserving group and ownser:
ssh ceto7 "rsync -rlptD --exclude=*.tmp --exclude=*.nc.* $TMPDIR/* $ARCHIVEDIR/" && rm -rf $TMPDIR/$folder
echo "Done."

