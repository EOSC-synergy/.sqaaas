#!/bin/bash


MPIRUN="mpirun $OPENQXD_MPIRUN_OPTS"
NAME=dirac

ALL="1 2 3 4 5 6 7 8 9"
WITHBC="1 2 3 4 5 6 7 8 9"
WITHGG="1 2 3 4 5 6 7 8 9"
WITHCS="1 2 3 4 5 6 7 8 9"


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


for CN in $ALL ; do
   bc[$CN]=0
   gg[$CN]=0
   cs[$CN]=0
done

for CN in $WITHBC ; do
   bc[$CN]=1
done

for CN in $WITHGG ; do
   gg[$CN]=1
done

for CN in $WITHCS ; do
   cs[$CN]=1
done

for CN in $ALL ; do
   if [ ${bc[$CN]} -eq 1 ] && [ ${gg[$CN]} -eq 1 ] && [ ${cs[$CN]} -eq 1 ] ; then
      for BC in 0 1 2 3 ; do
         for GG in 1 2 3 ; do
            for CS in 0 1 2 3 ; do
               $MPIRUN -np ${NPROC} ./check${CN} -bc ${BC} -gg ${GG} -cs ${CS}
               mv check${CN}.log ${LDIR}/check${CN}-bc${BC}-gg${GG}-cs${CS}-np${NPROC}.log
            done
         done
      done
   elif [ ${bc[$CN]} -eq 1 ] && [ ${gg[$CN]} -eq 0 ] && [ ${cs[$CN]} -eq 1 ] ; then
      for BC in 0 1 2 3 ; do
         for CS in 0 1 2 3 ; do
            $MPIRUN -np ${NPROC} ./check${CN} -bc ${BC} -cs ${CS}
            mv check${CN}.log ${LDIR}/check${CN}-bc${BC}-cs${CS}-np${NPROC}.log
         done
      done
   elif [ ${bc[$CN]} -eq 0 ] && [ ${gg[$CN]} -eq 1 ] && [ ${cs[$CN]} -eq 1 ] ; then
      for GG in 1 2 3 ; do
         for CS in 0 1 2 3 ; do
            $MPIRUN -np ${NPROC} ./check${CN} -gg ${GG} -cs ${CS}
            mv check${CN}.log ${LDIR}/check${CN}-gg${GG}-cs${CS}-np${NPROC}.log
         done
      done
   elif [ ${bc[$CN]} -eq 0 ] && [ ${gg[$CN]} -eq 0 ] && [ ${cs[$CN]} -eq 1 ] ; then
      for CS in 0 1 2 3 ; do
         $MPIRUN -np ${NPROC} ./check${CN} -cs ${CS}
         mv check${CN}.log ${LDIR}/check${CN}-cs${CS}-np${NPROC}.log
      done
   elif [ ${bc[$CN]} -eq 1 ] && [ ${gg[$CN]} -eq 1 ] && [ ${cs[$CN]} -eq 0 ] ; then
      for BC in 0 1 2 3 ; do
         for GG in 1 2 3 ; do
            $MPIRUN -np ${NPROC} ./check${CN} -bc ${BC} -gg ${GG}
            mv check${CN}.log ${LDIR}/check${CN}-bc${BC}-gg${GG}-np${NPROC}.log
         done
      done
   elif [ ${bc[$CN]} -eq 1 ] && [ ${gg[$CN]} -eq 0 ] && [ ${cs[$CN]} -eq 0 ] ; then
      for BC in 0 1 2 3 ; do
         $MPIRUN -np ${NPROC} ./check${CN} -bc ${BC}
         mv check${CN}.log ${LDIR}/check${CN}-bc${BC}-np${NPROC}.log
      done
   elif [ ${bc[$CN]} -eq 0 ] && [ ${gg[$CN]} -eq 1 ] && [ ${cs[$CN]} -eq 0 ] ; then
      for GG in 1 2 3 ; do
         $MPIRUN -np ${NPROC} ./check${CN} -gg ${GG}
         mv check${CN}.log ${LDIR}/check${CN}-gg${GG}-np${NPROC}.log
      done
   elif [ ${bc[$CN]} -eq 0 ] && [ ${gg[$CN]} -eq 0 ] && [ ${cs[$CN]} -eq 0 ] ; then
      $MPIRUN -np ${NPROC} ./check${CN}
      mv check${CN}.log ${LDIR}/check${CN}-np${NPROC}.log
   fi
done
