
/*******************************************************************************
*
* File cstar.c
*
* Copyright (C) 2016,2017 Patrick Fritzsch, Agostino Patella
*
* Implementation of cstar boundary conditions for the gauge fields.
*
* The externally accessible functions are
*
*   int get_mirror_rank(void)
*     It returns the MPI rank of the mirror process. The mirror process is
*     defined as follows. If cpr[k] for k=0,1,2,3 are the Cartesian coordinates
*     of the running process, the coordinates of the mirror process are mpr[k]
*     defined as mpr[k]=cpr[k] for k=0,2,3 and mpr[1]=(cpr[1]+NPROC1/2)%NPROC1.
*
* Notes:
*
* A process is said to be primary if cpr[1]<NPROC1/2, and secondary otherwise.
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define CSTAR_C

#include <string.h>
#include "global.h"
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"


int get_mirror_rank(void)
{
   static int mirror=-1;
   int mpr[4];

   error(iup[0][0]==0,1,"get_mirror_rank [cstar.c]",
         "Geometry arrays must be set first");

   if(mirror==-1)
   {
      mirror=ipr_global(cpr);
      if(bc_cstar()!=0)
      {
         mpr[0]=cpr[0];
         mpr[1]=(cpr[1]+NPROC1/2)%NPROC1;
         mpr[2]=cpr[2];
         mpr[3]=cpr[3];
         mirror=ipr_global(mpr);
      }
   }
   
   return mirror;
}
