
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency between U(1) and SU(3) plaquette sums.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)



int main(int argc,char *argv[])
{
   int my_rank,t,bc,cs,size;
   double su3phi[2],su3phi_prime[2];
   double u1phi,u1phi_prime;
   double qel;
   complex_dble *u1d,*u1db,*u1dm;
   su3_dble *ud;
   double u1plaq,u3plaq,u1sl[N0],u3sl[N0];
   double d,dmax,dmax_all;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);

      printf("\n");
      printf("Consistency between U(1) and SU(3) plaquette sums\n");
      printf("-------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-cs <cstar>]");
   }
   
   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   su3phi[0]=0.0;
   su3phi[1]=0.0;
   su3phi_prime[0]=0.0;
   su3phi_prime[1]=0.0;
   u1phi=0.573;
   u1phi_prime=-1.827;
   set_bc_parms(bc,cs,su3phi,su3phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   qel=1.111;
   set_u1lat_parms(0,4.348,1.0/qel,0.0,7.0,0.0,0.0,0);
   print_lat_parms();

   start_ranlux(0,12345);
   geometry();

   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
      size+=3;

   random_ad();
   u1db=u1dfld(EXT);
   u1dm=u1db+size;
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for (u1d=u1db;u1d<u1dm;u1d++)
   {
      (*ud).c11=(*u1d);
      (*ud).c22=(*u1d);
      (*ud).c33=(*u1d);
      ud++;
   }
   set_flags(UPDATED_UD);
   
   
   u1plaq=u1_plaq_sum_dble(0);
   u3plaq=plaq_sum_dble(0);
   dmax=fabs(u1plaq-u3plaq/3.);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between u1_plaq_sum_dble(0) and plaq_sum_dble(0), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   
   u1plaq=u1_plaq_sum_dble(1);
   u3plaq=plaq_sum_dble(1);
   dmax=fabs(u1plaq-u3plaq/3.);
   dmax_all=dmax;

   if (my_rank==0)
   {
      printf("Consistency between u1_plaq_sum_dble(1) and plaq_sum_dble(1), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   
   u1plaq=u1_plaq_wsum_dble(0);
   u3plaq=plaq_wsum_dble(0);
   dmax=fabs(u1plaq-u3plaq/3.);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between u1_plaq_wsum_dble(0) and plaq_wsum_dble(0), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   
   u1plaq=u1_plaq_wsum_dble(1);
   u3plaq=plaq_wsum_dble(1);
   dmax=fabs(u1plaq-u3plaq/3.);
   dmax_all=dmax;

   if (my_rank==0)
   {
      printf("Consistency between u1_plaq_wsum_dble(1) and plaq_wsum_dble(1), dev = %.2e\n",dmax_all);
      printf("\n");
   }


   u1plaq=u1_plaq_action_slices(u1sl);
   u3plaq=plaq_action_slices(u3sl);
   dmax=0.0;
   for (t=0;t<N0;t++)
   {
      d=fabs(2.*qel*qel*u1sl[t]-u3sl[t]/3.);
      if (d>dmax)
         dmax=d;
   }
   dmax_all=dmax;

   if (my_rank==0)
   {
      printf("Consistency between u1_plaq_action_slices and plaq_action_slices, dev = %.2e\n",dmax_all);
   }
   
   dmax=fabs(2.*qel*qel*u1plaq-u3plaq/3.);
   dmax_all=dmax;

   if (my_rank==0)
   {
      printf("Consistency between return values of u1_plaq_action_slices and plaq_action_slices, dev = %.2e\n",dmax_all);
      printf("\n");
   }


   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
