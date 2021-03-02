#!/bin/bash


MPIRUN="mpirun $OPENQXD_MPIRUN_OPTS"
NAME=cstates


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



for IBC in 0 1 2 3 ; do
for ICS in 1 2 3 ; do

   ${MPIRUN} -np ${NPROC} ./check1 -bc ${IBC} -cs ${ICS}
   mv check1.log ${LDIR}/check1-bc${IBC}-cs${ICS}.log

   ${MPIRUN} -np ${NPROC} ./check2 -bc ${IBC} -cs ${ICS}
   mv check2.log ${LDIR}/check2-bc${IBC}-cs${ICS}.log

   ${MPIRUN} -np ${NPROC} ./check2b -bc ${IBC} -cs ${ICS}
   mv check2b.log ${LDIR}/check2b-bc${IBC}-cs${ICS}.log

   ${MPIRUN} -np ${NPROC} ./check3 -bc ${IBC} -cs ${ICS}
   mv check3.log ${LDIR}/check3-bc${IBC}-cs${ICS}.log

done
done


for IBC in 0 1 2 3 ; do
for ICS in 0 1 2 3 ; do

   ${MPIRUN} -np ${NPROC} ./check4 -bc ${IBC} -cs ${ICS}
   mv check4.log ${LDIR}/check4-bc${IBC}-cs${ICS}.log

   ${MPIRUN} -np ${NPROC} ./check5 -bc ${IBC} -cs ${ICS}
   mv check5.log ${LDIR}/check5-bc${IBC}-cs${ICS}.log

done
done



for IF in 1 2 3 4 5 6 7 8 9 10 11 12 ; do

   cp check6.${IF}.in check6.in
   ${MPIRUN} -np ${NPROC} ./check6 
   mv check6.log ${LDIR}/check6-${IF}.log

done
cp check6.12.in check6.in
