#!/bin/bash


NAME=update


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



if [ -z "$QCDQED_CHECK_LOGDIR" ] ; then
   echo "Please define the environment variable QCDQED_CHECK_LOGDIR"
   exit -1
fi

if [ ! -d "$QCDQED_CHECK_LOGDIR" ] ; then
   echo "The environment variable QCDQED_CHECK_LOGDIR must contain the absolute path of a valid directory"
   exit -1
fi

if [[ "$QCDQED_CHECK_LOGDIR" != /* ]]; then
   echo "The environment variable QCDQED_CHECK_LOGDIR must contain the absolute path of a valid directory"
   exit -1
fi


LDIR=${QCDQED_CHECK_LOGDIR}/${NAME}
if [ ! -d "$LDIR" ] ; then
   mkdir $LDIR
fi

cd `dirname $0`



for IF in 1 2 3 ; do
   cp check1.${IF}.in check1.in
   mpirun -np ${NPROC} ./check1
   mv check1.log ${LDIR}/check1-if${IF}.log
done
cp check1.1.in check1.in

for IF in `seq 1 9` ; do
   cp check2.${IF}.in check2.in
   mpirun -np ${NPROC} ./check2
   mv check2.log ${LDIR}/check2-if${IF}.log
done
cp check2.1.in check2.in

for IF in `seq 1 9` ; do
   cp check3.${IF}.in check3.in
   mpirun -np ${NPROC} ./check3
   mv check3.log ${LDIR}/check3-if${IF}.log
done
cp check3.1.in check3.in

mpirun -np ${NPROC} ./check4
mv check4.log ${LDIR}/check4.log

for BC in 0 1 2 3 ; do
   mpirun -np ${NPROC} ./check5 -bc ${BC}
   mv check5.log ${LDIR}/check5-bc${BC}.log
done

for BC in 0 1 2 3 ; do
   mpirun -np ${NPROC} ./check6 -bc ${BC}
   mv check6.log ${LDIR}/check6-bc${BC}.log
done
