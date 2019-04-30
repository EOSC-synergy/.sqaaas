#!/bin/bash


MPIRUN="mpirun $OPENQXD_MPIRUN_OPTS"
NAME=archive


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


for BC in 0 1 2 3 ; do
   for CS in 0 1 2 3 ; do
      for GG in 1 2 3 ; do
         $MPIRUN -np ${NPROC} ./check1 -bc ${BC} -gg ${GG} -cs ${CS}
         mv check1.log ${LDIR}/check1-bc${BC}-cs${CS}-gg${GG}-np${NPROC}.log
      done
   done
done

for BC in 0 1 2 3 ; do
   for CS in 0 1 2 3 ; do
      $MPIRUN -np ${NPROC} ./check2 -bc ${BC} -gg 1 -cs ${CS}
      mv check2.log ${LDIR}/check2-bc${BC}-cs${CS}-gg1-np${NPROC}.log
      $MPIRUN -np ${NPROC} ./check3 -bc ${BC} -gg 1 -cs ${CS}
      mv check3.log ${LDIR}/check3-bc${BC}-cs${CS}-gg1-gg1-np${NPROC}.log

      $MPIRUN -np ${NPROC} ./check2 -bc ${BC} -gg 1 -cs ${CS}
      mv check2.log ${LDIR}/check2-bc${BC}-cs${CS}-gg1-np${NPROC}.log
      $MPIRUN -np ${NPROC} ./check3 -bc ${BC} -gg 3 -cs ${CS}
      mv check3.log ${LDIR}/check3-bc${BC}-cs${CS}-gg1-gg3-np${NPROC}.log

      $MPIRUN -np ${NPROC} ./check2 -bc ${BC} -gg 2 -cs ${CS}
      mv check2.log ${LDIR}/check2-bc${BC}-cs${CS}-gg2-np${NPROC}.log
      $MPIRUN -np ${NPROC} ./check3 -bc ${BC} -gg 2 -cs ${CS}
      mv check3.log ${LDIR}/check3-bc${BC}-cs${CS}-gg2-gg2-np${NPROC}.log

      $MPIRUN -np ${NPROC} ./check2 -bc ${BC} -gg 2 -cs ${CS}
      mv check2.log ${LDIR}/check2-bc${BC}-cs${CS}-gg2-np${NPROC}.log
      $MPIRUN -np ${NPROC} ./check3 -bc ${BC} -gg 3 -cs ${CS}
      mv check3.log ${LDIR}/check3-bc${BC}-cs${CS}-gg2-gg3-np${NPROC}.log

      $MPIRUN -np ${NPROC} ./check2 -bc ${BC} -gg 3 -cs ${CS}
      mv check2.log ${LDIR}/check2-bc${BC}-cs${CS}-gg3-np${NPROC}.log
      $MPIRUN -np ${NPROC} ./check3 -bc ${BC} -gg 3 -cs ${CS}
      mv check3.log ${LDIR}/check3-bc${BC}-cs${CS}-gg3-gg3-np${NPROC}.log
   done
done

for CN in 4 5 6 ; do
   $MPIRUN -np ${NPROC} ./check${CN}
   mv check${CN}.log ${LDIR}/check${CN}-np${NPROC}.log
done
