
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2016, 2017 Agostino Patella
* 
* Based on openQCD-1.6/devel/uflds/check1.c
* Copyright (C) 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of the comm buffers of the compact U(1) field.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "su3fcts.h"
#include "u1flds.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)



int main(int argc,char *argv[])
{
   int my_rank,bc,cs,ie;
   size_t size;
   double d,dmax,dmax_all;
   double su3phi[2],su3phi_prime[2];
   double u1phi,u1phi_prime;
   complex_dble *u1d,*u1db,*u1dm;
   su3_dble *ud;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);

      printf("\n");
      printf("Calculation of the comm buffers of the compact U(1) field\n");
      printf("---------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>] [-cs <cstar>]");
   }
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   su3phi[0]=0.573;
   su3phi[1]=-0.573;
   su3phi_prime[0]=-1.827;
   su3phi_prime[1]=1.827;
   u1phi=0.573;
   u1phi_prime=-1.827;
   set_bc_parms(bc,cs,su3phi,su3phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   start_ranlux(0,123456);
   geometry();
   
   random_ad();
   
   u1db=u1dfld(LOC);
   u1dm=u1db+4*VOLUME;
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for (u1d=u1db;u1d<u1dm;u1d++)
   {
      (*ud).c11=(*u1d);
      (*ud).c22.re=(*u1d).re;
      (*ud).c22.im=-(*u1d).im;
      if(((*u1d).re!=0)||((*u1d).im!=0.0))
         (*ud).c33.re=1.0;
      ud++;
   }
   set_flags(UPDATED_UD);
   copy_bnd_ud();


   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
      size+=3;

   dmax=0.0;
   u1db=u1dfld(EXT);
   u1dm=u1db+size;
   ud=udfld();
   for (u1d=u1db;u1d<u1dm;u1d++)
   {
      d=fabs((*ud).c11.re-(*u1d).re);
      if (d>dmax)
         dmax=d;
      d=fabs((*ud).c11.im-(*u1d).im);
      if (d>dmax)
         dmax=d;
      ud++;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between u1d field and its boundaries, dev = %.2e\n",dmax_all);
      printf("\n");
   }

   print_flags();

   ie=check_bc(1.0e-15);
   error_root(ie==0,1,"main [check2.c]","U(1) boundary conditions are not consistent with SU(3) boundary conditions");


   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();
   exit(0);
}
