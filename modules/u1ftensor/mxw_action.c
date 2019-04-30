
/*******************************************************************************
*
* File mxw_action.c
*
* Copyright (C) 2016, 2017 Nazario Tantalo
*
* Based on openQCD-1.6/modules/tcharge/ym_action.c
* Copyright (C) 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the Maxwell action using the symmetric U(1) field tensor.
*
* The externally accessible functions are
*
*   double mxw_action(void)
*     Returns the Maxwell action S (w/o prefactor 1/e^2) of the double-precision 
*     gauge field, using a symmetric expression for the gauge-field tensor.
*
*   double mxw_action_slices(double *asl)
*     Computes the sum asl[t] of the Maxwell action density (w/o prefactor 1/e^2) 
*     of the double-precision gauge field at time t=0,1,...,N0-1 
*     (where N0=NPROC0*L0). The program returns the total action.
*
* Notes:
*
* The Maxwell action density s(x) is defined by
*
*  s(x)=(1/4)*sum_{mu nu} [Ahat_{mu nu}(x)]^2
*
* where Ahat_{mu nu}(x) are the components of the symmetric U(1) field tensor 
* returned by the program u1ftensor() [u1ftensor.c]. At the boundaries of the 
* lattice (if any), the action density is set to zero. The total action S is 
* the sum of s(x) over all points x with time component in the range
*
*  0<x0<NPROC0*L0-1        (open bc),
*
*  0<x0<NPROC0*L0          (SF and open-SF bc),
*
*  0<=x0<NPROC0*L0         (periodic bc).
*
* The programs in this module perform global operations and must be called
* simultaneously on all MPI processes.
*
*******************************************************************************/

#define MXW_ACTION_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "linalg.h"
#include "u1ftensor.h"
#include "global.h"


#define N0 (NPROC0*L0)

static int isx[L0],init=0;
static double asl0[N0];
static double **ft;


static double density(int ix)
{
   double sm;

   sm=ft[0][ix]*ft[0][ix]+ft[1][ix]*ft[1][ix]+ft[2][ix]*ft[2][ix]+
      ft[3][ix]*ft[3][ix]+ft[4][ix]*ft[4][ix]+ft[5][ix]*ft[5][ix];

   return sm;
}


double mxw_action(void)
{
   int bc,ix,t,tmx;
   double S;

   if (init==0)
   {
      for (t=0;t<L0;t++)
         isx[t]=init_hsum(1);

      init=1;
   }

   ft=u1ftensor();
   bc=bc_type();
   if (bc==0)
      tmx=N0-1;
   else
      tmx=N0;
   reset_hsum(isx[0]);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (((t>0)&&(t<tmx))||(bc==3))
      {
         S=density(ix);
         add_to_hsum(isx[0],&S);
      }
   }

   if (NPROC>1)
      global_hsum(isx[0],&S);
   else
      local_hsum(isx[0],&S);

   return 0.5*S;
}


double mxw_action_slices(double *asl)
{
   int bc,ix,t,t0,tmx;
   double S;

   if (init==0)
   {
      for (t=0;t<L0;t++)
         isx[t]=init_hsum(1);

      init=1;
   }

   ft=u1ftensor();
   bc=bc_type();
   if (bc==0)
      tmx=N0-1;
   else
      tmx=N0;
   t0=cpr[0]*L0;

   for (t=0;t<L0;t++)
      reset_hsum(isx[t]);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (((t>0)&&(t<tmx))||(bc==3))
      {
         t-=t0;
         S=density(ix);
         add_to_hsum(isx[t],&S);
      }
   }

   for (t=0;t<N0;t++)
      asl0[t]=0.0;

   for (t=0;t<L0;t++)
   {
      local_hsum(isx[t],&S);
      asl0[t+t0]=0.5*S;
   }

   if (NPROC>1)
   {
      MPI_Reduce(asl0,asl,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(asl,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
         asl[t]=asl0[t];
   }

   S=0.0;

   for (t=0;t<N0;t++)
      S+=asl[t];

   return S;
}

#undef N0
