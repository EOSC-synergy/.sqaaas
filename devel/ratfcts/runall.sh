#!/bin/bash


MPIRUN="mpirun $OPENQXD_MPIRUN_OPTS"
NAME=ratfcts


if [ $# -ne 1 ] ; then
  echo "Usage: `basename $0` NPROC"
  exit 85
fi

declare -i NPROC
NPROC=$1

if [ "$NPROC" -le 0 ] ; then
   echo "Usage: `basename $0` NPROC"
   echo "with NPROC>0 !!!"
   exit 86
fi



if [ -z "$OPENQXD_CHECK_LOGDIR" ] ; then
   echo "Please define the environment variable OPENQXD_CHECK_LOGDIR"
   exit -1
fi

if [ ! -d "$OPENQXD_CHECK_LOGDIR" ] ; then
   echo "The environment variable OPENQXD_CHECK_LOGDIR must contain the absolute path of a valid directory"
   exit -1
fi

if [[ "$OPENQXD_CHECK_LOGDIR" != /* ]]; then
   echo "The environment variable OPENQXD_CHECK_LOGDIR must contain the absolute path of a valid directory"
   exit -1
fi


LDIR=${OPENQXD_CHECK_LOGDIR}/${NAME}
if [ ! -d "$LDIR" ] ; then
   mkdir $LDIR
fi

cd `dirname $0`


for IF in 1 2 3 4 5 6 7 ; do
   cp check1.${IF}.in check1.in
   $MPIRUN -np ${NPROC} ./check1
   mv check1.log ${LDIR}/check1-if${IF}-np${NPROC}.log
done
cp check1.1.in check1.in
