
/*******************************************************************************
*
* File check7.c
*
* Copyright (C) 2016 Nazario Tantalo
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks between U(1) and SU(3) field tensors
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "u1ftensor.h"
#include "tcharge.h"
#include "global.h"
#include "gflds_utils.h"


#define N0 (NPROC0*L0)

static int isx[L0],init=0;
static double **u1ft;
static u3_alg_dble **ft;


void compare(double *diff)
{
   int bc,tmx;
   int ix,t,munu;
   double tdiff,iqel;
   u1lat_parms_t u1lp;
   
   u1lp=u1lat_parms();
   iqel=u1lp.invqel;

   if (init==0)
   {
      for (t=0;t<L0;t++)
	 isx[t]=init_hsum(1);
      
      init=1;
   }
   
   bc=bc_type();
   if (bc==0)
      tmx=N0-1;
   else
      tmx=N0;

   for (munu=0;munu<6;munu++)
   {
      reset_hsum(isx[0]);
      reset_hsum(isx[1]);

      for (ix=0;ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (((t>0)&&(t<tmx))||(bc==3))
         {
            tdiff=fabs(u1ft[munu][ix]/iqel-ft[munu][ix].c1);
            add_to_hsum(isx[0],&tdiff);

            tdiff=fabs(u1ft[munu][ix]/iqel);
            add_to_hsum(isx[1],&tdiff);
         }
      }

      if (NPROC>1)
      {
         global_hsum(isx[0],diff+munu);
         global_hsum(isx[1],&tdiff);
      }
      else
      {
         local_hsum(isx[0],diff+munu);
         local_hsum(isx[1],&tdiff);
      }

      diff[munu]/=tdiff;

   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc,munu;
   double phi[2],phi_prime[2];
   double diff[6],mxw,ym,iqel2;
   complex_dble *u1d,*u1db,*u1dm,cone;
   su3_dble *ud;
   FILE *flog=NULL;
   u1lat_parms_t u1lp;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check7.log","w",stdout);
      printf("\n");
      printf("Consistency checks between U(1) and SU(3) field tensors\n");
      printf("-------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check7.c]",
                    "Syntax: check7 [-bc <type>]");
   }

   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();

   set_u1lat_parms(0,1.0/137.0,5.0,0.0,7.0,0.0,0.0);
   print_lat_parms();

   u1lp=u1lat_parms();
   iqel2 =u1lp.invqel;
   iqel2*=u1lp.invqel;

   start_ranlux(0,12345);
   geometry();

   cone.re=1.0;
   cone.im=0.0;

   random_ad();
   u1db=u1dfld(LOC);
   u1dm=u1db+4*VOLUME;
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for (u1d=u1db;u1d<u1dm;u1d++)
   {
      (*ud).c11=(*u1d);
      (*ud).c22=cone;
      (*ud).c33=cone;
      ud++;
   }
   set_flags(UPDATED_UD);

   u1ft=u1ftensor();
   ft=ftensor();
   
   mxw=mxw_action();
   ym=ym_action();

   mxw/=iqel2;
   ym=3.0*ym/4.0;

   compare(diff);
   
   if (my_rank==0)
   {
      printf("\n");
      printf("dev_munu= {sum_x abs(u1ft[munu][x]-ft[munu][x].c1)} / {sum_x abs(u1ft[munu][x])}\n\n");
      printf("Relative deviation (munu,dev_munu)                  = ");
      for (munu=0;munu<6;munu++)
	 printf("(%d, %.1e) ",munu,diff[munu]);
      printf("\n");

      printf("Relative deviation between Maxwell and Y.M. actions = %.1e\n\n",fabs(1.0-mxw/ym));
      fclose(flog);
   }


   MPI_Finalize();
   exit(0);
}

#undef N0
