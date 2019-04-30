
/*******************************************************************************
*
* File ud_cstar.c
*
* Copyright (C) 2016,2017 Patrick Fritzsch, Agostino Patella
*
* Implementation of cstar boundary conditions for the SU(3) gauge field.
*
* The externally accessible functions are
*
*   void orbi_cpy_ud(void)
*     The SU(3) gauge field on each secondary process is set to be
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
* The cstar-transform of a group matrix is defined as its complex conjugate.
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define UD_CSTAR_C

#include <string.h>
#include <stdio.h>
#include "global.h"
#include "mpi.h"
#include "su3.h"
#include "uflds.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static void cstar_su3_dble(int vol,su3_dble* u)
{
  su3_dble *um; 

  um=u+vol;
  for (;u<um;u++)
  {    
     (*u).c11.im *= -1.0;
     (*u).c12.im *= -1.0;
     (*u).c13.im *= -1.0;
     (*u).c21.im *= -1.0;
     (*u).c22.im *= -1.0;
     (*u).c23.im *= -1.0;
     (*u).c31.im *= -1.0;
     (*u).c32.im *= -1.0;
     (*u).c33.im *= -1.0;
  }
}


void orbi_cpy_ud(void)
{
   int mirror, tag;
   MPI_Status stat;
   su3_dble *ud;
   
   if(bc_cstar()>0) {
      ud=udfld();
      mirror=get_mirror_rank();
      tag=mpi_tag();
      if(cpr[1]<NPROC1/2) {
         MPI_Send(ud,18*4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD);
      } else {
         MPI_Recv(ud,18*4*VOLUME,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD,&stat);
         cstar_su3_dble(4*VOLUME,ud);
      }
      set_flags(UPDATED_UD);
   }
}
