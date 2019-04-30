#!/bin/bash


NAME=dirac

ALL="1 2 3 4 5 6 7 8 9"
WITHBC="1 2 3 4 5 6 7 8 9"
WITHCF="1 2 3 4 5 6 7 8 9"


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


for CN in $ALL ; do
   bc[$CN]=0
   cf[$CN]=0
done

for CN in $WITHBC ; do
   bc[$CN]=1
done

for CN in $WITHCF ; do
   cf[$CN]=1
done

for CN in $ALL ; do
   if [ ${bc[$CN]} -eq 1 ] && [ ${cf[$CN]} -eq 1 ] ; then
      for BC in 0 1 2 3 ; do
         for FL in 1 2 3 ; do
            mpirun -np ${NPROC} ./check${CN} -bc ${BC} -gg ${FL}
            mv check${CN}.log ${LDIR}/check${CN}-bc${BC}-gg${FL}-np${NPROC}.log
         done
      done
   elif [ ${bc[$CN]} -eq 1 ] && [ ${cf[$CN]} -eq 0 ] ; then
      for BC in 0 1 2 3 ; do
         mpirun -np ${NPROC} ./check${CN} -bc ${BC}
         mv check${CN}.log ${LDIR}/check${CN}-bc${BC}-np${NPROC}.log
      done
   elif [ ${bc[$CN]} -eq 0 ] && [ ${cf[$CN]} -eq 1 ] ; then
      for FL in 1 2 3 ; do
         mpirun -np ${NPROC} ./check${CN} -gg ${FL}
         mv check${CN}.log ${LDIR}/check${CN}-gg${FL}-np${NPROC}.log
      done
   elif [ ${bc[$CN]} -eq 0 ] && [ ${cf[$CN]} -eq 0 ] ; then
      mpirun -np ${NPROC} ./check${CN}
      mv check${CN}.log ${LDIR}/check${CN}-np${NPROC}.log
   fi
done
