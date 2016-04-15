#!/bin/bash 
I=$2
J=$3
IJ=$4
y=$5
m=$6
nit0=$7

RUNDIR=/work/n01/n01/slwa/NEMO//src/NEMO_V3.6_STABLE_top_bdy/NEMO/NEMOGCM/CONFIG/XIOS_AMM7_top/EXP01
if [ $m -lt 10 ] ; then
  month=0$m
else
  month=$m
fi

echo qsub -v I=$I,J=$J,IJ=$IJ,m=$m,y=$y,nit0=$nit0  -o $RUNDIR/NE-AMM7-$IJ-$y-$month -N NE$y$month $1
qsub -v I=$I,J=$J,IJ=$IJ,m=$m,y=$y,nit0=$nit0 -o $RUNDIR/NE-AMM7-$IJ-$y-$month -N NE$y$month  $1



#usage
#  ./run_nemo.sh monthlyrun_restart.pbs 12 16 192 1995 1 1
#  ./run_nemo.sh annualrun.pbs 12 16 192 1995 1 1
# run perpetuates
