
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

static double **u1ft;
static u3_alg_dble **ft;


double compare(void)
{
   int ix,munu;
   double d,diff,iqel;
   u1lat_parms_t u1lp;
   
   u1lp=u1lat_parms();
   iqel=u1lp.invqel;

   diff=0.0;
   for (munu=0;munu<6;munu++)
   {
      for (ix=0;ix<VOLUME;ix++)
      {
         d=fabs(u1ft[munu][ix]/iqel-ft[munu][ix].c1);
         if(d>diff)
            diff=d;
      }
   }
   
   d=diff;
   MPI_Reduce(&d,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   return diff;
}


int main(int argc,char *argv[])
{
   int my_rank,bc;
   double su3phi[2],su3phi_prime[2];
   double u1phi,u1phi_prime;
   double diff,mxw,ym,iqel2;
   complex_dble *u1d,*u1db,*u1dm;
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

   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   su3phi[0]=0.573;
   su3phi[1]=-0.573;
   su3phi_prime[0]=-1.827;
   su3phi_prime[1]=1.827;
   u1phi=0.573;
   u1phi_prime=-1.827;
   set_bc_parms(bc,0,su3phi,su3phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   set_u1lat_parms(0,1.0/137.0,5.0,0.0,7.0,0.0,0.0,0);
   print_lat_parms();

   u1lp=u1lat_parms();
   iqel2 =u1lp.invqel;
   iqel2*=u1lp.invqel;

   start_ranlux(0,12345);
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

   u1ft=u1ftensor();
   ft=ftensor();
   
   mxw=mxw_action();
   ym=ym_action();

   mxw/=iqel2;
   ym=ym/4.0;

   diff=compare();
   
   if (my_rank==0)
   {
      printf("\n");
      printf("Maximum deviation between SU(3) and U(1) field tensors = %.1e\n\n",diff);
      printf("Relative deviation between Maxwell and Y.M. actions = %.1e\n\n",fabs(1.0-mxw/ym));
      fclose(flog);
   }


   MPI_Finalize();
   exit(0);
}

#undef N0
