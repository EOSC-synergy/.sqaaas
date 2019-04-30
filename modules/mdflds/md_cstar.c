/*******************************************************************************
*
* File md_cstar.c
*
* Copyright (C) 2016,2017 Patrick Fritzsch, Agostino Patella
*
* Implementation of cstar boundary conditions for the algebra (momentum, force
* and field tensor) fields.
*
* The externally accessible functions are
*
*   void orbi_cpy_su3mom(void)
*     The SU(3) momentum field on each secondary process is set to be
*     equal to the cstar-transform of the field on the mirror process.
*
*   void orbi_cpy_u1mom(void)
*     The U(1) momentum field on each secondary process is set to be
*     equal to the cstar-transform of the field on the mirror process.
*
* Notes:
*
* The mirror process is defined as follows. If cpr[k] for k=0,1,2,3 are the
* Cartesian coordinates of the running process, the coordinates of the mirror
* process are mpr[k] defined as mpr[k]=cpr[k] for k=0,2,3 and
* mpr[1]=(cpr[1]+NPROC1/2)%NPROC1
*
* A process is said to be primary if cpr[1]<NPROC1/2, and secondary otherwise.
*
* The cstar transform of an algebra matrix is defined as its complex conjugate.
* The components of a su3 algebra matrix transform as
*     Ak = -Ak for k=1,2,4,6,8
*     Ak = Ak  for k=3,5,7
* The u1 algebra element transforms as
*     A = -A
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes. We assume that C* boundary
* conditions are present and the topology has been set up already. We do 
* not additionally check it again.
*
*******************************************************************************/

#define MD_CSTAR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "linalg.h"
#include "mdflds.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static void cstar_su3alg_dble(const int vol,su3_alg_dble* f3)
{
  su3_alg_dble *fm; 

  fm=f3+vol;
  for (;f3<fm;f3++)
  {    
     (*f3).c1 *= -1;
     (*f3).c2 *= -1;
     (*f3).c4 *= -1;
     (*f3).c6 *= -1;
     (*f3).c8 *= -1;
  }
}


static void cstar_double(int vol,double* f1)
{
  double *fm; 

  fm=f1+vol;
  for (;f1<fm;f1++)  (*f1) *= -1.0;
}


void orbi_cpy_su3mom(void)
{
   int mirror, tag;
   MPI_Status stat;
   su3_alg_dble *mom;
   
   if(bc_cstar()>0) {
      mom=mdflds()->su3mom;
      mirror=get_mirror_rank();
      tag=mpi_tag();
      if(cpr[1]<NPROC1/2)
      {
         MPI_Send(mom,8*4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(mom,8*4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD,&stat);
         cstar_su3alg_dble(4*VOLUME,mom);
      }
   }
}


void orbi_cpy_u1mom(void) {
   int mirror, tag;
   MPI_Status stat;
   double *mom;
   
   if(bc_cstar()>0)
   {
      mom=mdflds()->u1mom;
      mirror=get_mirror_rank();
      tag=mpi_tag();
      if(cpr[1]<NPROC1/2)
      {
         MPI_Send(mom,4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD);
      }
      else
      {
         MPI_Recv(mom,4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD,&stat);
         cstar_double(4*VOLUME,mom);
      }
   }
}
