
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the inversion program inv_nabla_sq().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "u1flds.h"
#include "u1ftensor.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "cstates.h"
#include "global.h"
#include "gflds_utils.h"
#include "sflds_utils.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,ix,t;
   double phi[2],phi_prime[2];
   double d,d1[NPROC0*L0],d1sum[NPROC0*L0];
   double *phi1,*phi2,*phi3;
   double relres;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);
      printf("\n");
      printf("Check of the inversion program inv_nabla_sq()\n");
      printf("---------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check5.c]",
                     "Syntax: check5 [-cs <cstar>] [-bc <type>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check5.c]",
                     "Syntax: check5 [-cs <cstar>] [-bc <type>]");
   }

   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();

   phi1=amalloc(3*(VOLUME+BNDRY)*sizeof(*phi1),4);
   error(phi1==NULL,1,"main [check5.c]",
         "Unable to allocate service array");

   phi2=phi1+VOLUME+BNDRY;
   phi3=phi2+VOLUME+BNDRY;

   random_dvec(VOLUME,phi1);

   for (t=0;t<NPROC0*L0;++t)
      d1[t]=0.0;

   for(ix=0;ix<VOLUME;++ix)
   {
      t=global_time(ix);
      d1[t]+=phi1[ix];
   }

   for (t=0;t<NPROC0*L0;++t)
      d1[t]/=(double)(L1*L2*L3*NPROC1*NPROC2*NPROC3);

   if (NPROC>1)
   {
      MPI_Reduce(d1,d1sum,NPROC0*L0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(d1sum,NPROC0*L0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<NPROC0*L0;t++)
      d1sum[t]=d1[t];
   }

   for(ix=0;ix<VOLUME;++ix)
   {
      t=global_time(ix);
      phi1[ix]-=d1sum[t];
   }

   relres=inv_nabla_sq(phi2,phi1);
   nabla_sq_dvec(phi2,phi3);


   muladd_assign_dvec(VOLUME,-1.0,phi1,phi3);
   d=norm_square_dvec(VOLUME,1,phi3)/norm_square_dvec(VOLUME,1,phi1);

   if (my_rank==0)
   {   
      printf("|{1- nabla^2/nabla^2} phi| / |phi| = %.2e\n",sqrt(d));
      printf("(should smaller than %.2e)\n\n",relres);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
