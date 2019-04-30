
/*******************************************************************************
*
* File u1fluxes.c
*
* Copyright (C) 2016, 2017 Nazario Tantalo
*
* Based on openQCD-1.6/modules/tcharge/tcharge.c
* Copyright (C) 2010-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the electromagnetic fluxes using the symmetric U(1) 
* field tensor.
*
* The externally accessible functions are
*
*   double u1fluxes(int munu)
*     Calculates the "field-theoretic" electromagnetic fluxes of the global
*     double-precision U(1) gauge field, using a symmetric expression for the
*     gauge-field tensor. The function returns Phi_{mu nu} where the values of
*     mu and nu corresponding to munu= 0,1,2,3,4,5 are (mu nu)= {0 1}, {0 2},
*     {0 3}, {2 3}, {3 1}, {1 2}.
*
*   double u1fluxes_slices(int munu,double *Phisl)
*     Calculates the "field-theoretic" electromagnetic fluxes Phi_{mu nu}(x0) of
*     the global double-precision U(1) gauge field, using a symmetric expression
*     for the gauge-field tensor. The variables Phisl[x0] for
*     x0=0,...,NPROC0*L0-1 store the fluxes Phi_{mu nu}(x0) where the values of
*     mu and nu corresponding to munu= 0,1,2,3,4,5 are (mu nu)= {0 1}, {0 2},
*     {0 3}, {2 3}, {3 1}, {1 2}. In the case of open b.c. Phisl[NPROC0*L0-1]=0.
*     The function returns Phi_{mu nu}.
*
* Notes:
*
* The U(1) elctromagnetic fluxes are defined as the sum over all points
* x=(x0,vecx) with time component in the range
*
*  0<x0<TMX,    TMX=NPROC0*L0-1       (open bc),
*
*  0<x0<TMX,    TMX=NPROC0*L0         (SF and open-SF bc),
*
*  0<=x0<TMX,   TMX=NPROC0*L0         (periodic bc).
*
* of the field tensor Ahat_{mu nu}(x), returned by  the program u1ftensor() 
* [u1ftensor.c], with appropriate normalization factors, 
*
*  Phi_{mu nu}(x0) = (2*Pi*Lp_{munu})^{-1} * sum_{vecx} Ahat_{mu nu}(x) ,
*
*  Phi_{mu nu} = sum_{x0} Phi_{mu nu}(x0) ,
*
* where
*
*  Lp_{01} = NPROC2*L2*NPROC3*L3
*  Lp_{02} = NPROC1*L1*NPROC3*L3
*  Lp_{03} = NPROC1*L1*NPROC2*L2
*
*  Lp_{23} = TMX*NPROC1*L1
*  Lp_{31} = TMX*NPROC2*L2
*  Lp_{23} = TMX*NPROC1*L1
*
* In case of C* boundary conditions the sum is performed only over the points
* in the fundamental domain, and NPROC1 is replaced by NPROC1/2 in the above
* formulae.
*
* The programs in this module perform global communications and must be
* called simultaneously on all processes.
*
*******************************************************************************/

#define U1FLUXES_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "tcharge.h"
#include "u1ftensor.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int isx[L0],init=0;
static double msl0[N0];
static double **ft;


double u1fluxes(int munu)
{
   int bc,tmx;
   int ix,t;
   double twopi,Phi,norm[6];

   check_global_int("u1fluxes",1,munu);   
   error((munu<0)||(munu>5),1,"u1fluxes [u1fluxes.c]","munu has to be in [0,5]");

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

   norm[0]=(double)(N2*N3);
   norm[1]=(double)(N1*N3);
   norm[2]=(double)(N1*N2);
   norm[3]=(double)(tmx*N1);
   norm[4]=(double)(tmx*N2);
   norm[5]=(double)(tmx*N3);
   
   if (bc_cstar()!=0)
   {
      norm[1]*=0.5;
      norm[2]*=0.5;
      norm[3]*=0.5;
   }

   twopi=8.0*atan(1.0);

   reset_hsum(isx[0]);

   if ((bc_cstar()==0)||(cpr[1]<NPROC1/2))
   {
      for (ix=0;ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (((t>0)&&(t<tmx))||(bc==3))
            add_to_hsum(isx[0],ft[munu]+ix);
      }
   }

   if (NPROC>1)
      global_hsum(isx[0],&Phi);
   else
      local_hsum(isx[0],&Phi);

   Phi/= twopi*norm[munu];

   return Phi;
}


double u1fluxes_slices(int munu,double *Phisl)
{
   int bc,tmx;
   int ix,t,t0;
   double twopi,Phi,norm[6];

   check_global_int("u1fluxes_slices",1,munu);   
   error((munu<0)||(munu>5),1,"u1fluxes_slices [u1fluxes.c]","munu has to be in [0,5]");

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

   norm[0]=(double)(N2*N3);
   norm[1]=(double)(N1*N3);
   norm[2]=(double)(N1*N2);
   norm[3]=(double)(tmx*N1);
   norm[4]=(double)(tmx*N2);
   norm[5]=(double)(tmx*N3);

   if (bc_cstar()!=0)
   {
      norm[1]*=0.5;
      norm[2]*=0.5;
      norm[3]*=0.5;
   }

   twopi=8.0*atan(1.0);

   for (t=0;t<L0;t++)
      reset_hsum(isx[t]);

   t0=cpr[0]*L0;

   if ((bc_cstar()==0)||(cpr[1]<NPROC1/2))
   {
      for (ix=0;ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (((t>0)&&(t<tmx))||(bc==3))
         {
            t-=t0;
            add_to_hsum(isx[t],ft[munu]+ix);
         }
      }
   }

   for (t=0;t<N0;t++)
      msl0[t]=0.0;

   for (t=0;t<L0;t++)
   {
      local_hsum(isx[t],&Phi);
      msl0[t+t0]=Phi/(twopi*norm[munu]);
   }

   if (NPROC>1)
   {
      MPI_Reduce(msl0,Phisl,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(Phisl,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
         Phisl[t]=msl0[t];
   }

   Phi=0.0;   
   for (t=0;t<N0;t++)
      Phi+=Phisl[t];

   return Phi;
}

#undef N0
#undef N1
#undef N2
#undef N3
