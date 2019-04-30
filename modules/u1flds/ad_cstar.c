
/*******************************************************************************
*
* File ad_cstar.c
*
* Copyright (C) 2016,2017 Patrick Fritzsch, Agostino Patella
*
* Implementation of cstar boundary conditions for the U(1) gauge field.
*
* The externally accessible functions are
*
*   void orbi_cpy_ad(void)
*     The non-compact U(1) gauge field on each secondary process is set to be
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
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define AD_CSTAR_C

#include <string.h>
#include <stdio.h>
#include "global.h"
#include "mpi.h"
#include "u1flds.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static void cstar_double(int vol,double* ad)
{
  double *am; 

  am=ad+vol;
  for (;ad<am;ad++)  (*ad) *= -1.0;
}

void orbi_cpy_ad(void)
{
   int mirror, tag;
   MPI_Status stat;
   double *ad;
   
   if(bc_cstar()>0) {
      ad=adfld();
      mirror=get_mirror_rank();
      tag=mpi_tag();
      if(cpr[1]<NPROC1/2) {
         MPI_Send(ad,4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD);
      } else {
         MPI_Recv(ad,4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD,&stat);
         cstar_double(4*VOLUME,ad);
      }
      set_flags(UPDATED_AD);
   }
}
