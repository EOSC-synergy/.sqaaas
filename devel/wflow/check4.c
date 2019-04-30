
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency between U(1) and SU(3) Wilson flows.
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
#include "mdflds.h"
#include "wflow.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)



int main(int argc,char *argv[])
{
   int my_rank,bc,i,n;
   double alpha,beta,qel2,phi[2],phi_prime[2],eps;
   double *ad;
   su3_dble *ud;
   double d,dmax,dmax_all;
   FILE *flog=NULL;
   FILE *fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      fin=freopen("check4.in","r",stdin);

      printf("\n");
      printf("Consistency between compact U(1) and SU(3) gauge forces and actions\n");
      printf("-------------------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("n","%d\n",&n);
      read_line("eps","%lf",&eps);
      fclose(fin);

      printf("n = %d\n",n);
      printf("eps = %.3e\n\n",eps);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>]");
   }

   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();

   alpha=.013;
   qel2=2.45;
   beta=1./(16.*atan(1)*qel2*alpha);
   set_u1lat_parms(0,alpha,1.0/sqrt(qel2),0.0,0.482,0.87,0.57);
   set_su3lat_parms(beta,0.482,0.87,0.57);
   print_lat_parms();

   start_ranlux(0,12345);
   geometry();
   alloc_wf3d(1);
   alloc_wf1d(1);
   
   
   random_ad();
   ad=adfld();
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for(i=0;i<4*VOLUME;i++)
   {
      ud[i].c11.re=ud[i].c22.re=cos(ad[i]);
      ud[i].c12.re=sin(ad[i]);
      ud[i].c21.re=-ud[i].c12.re;
      ud[i].c33.re=1.0;
   }
   set_flags(UPDATED_UD);
   
   fwd_su3_euler(n,eps);
   fwd_u1_euler(n,eps*qel2);
   
   dmax=0.0;
   for(i=0;i<4*VOLUME;i++)
   {
      d=fabs(ud[i].c11.re-cos(ad[i]));
      if(d>dmax) dmax=d;
      d=fabs(ud[i].c12.re-sin(ad[i]));
      if(d>dmax) dmax=d;
   }
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between fwd_su3_euler() and fwd_u1_euler(), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   
   random_ad();
   ad=adfld();
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for(i=0;i<4*VOLUME;i++)
   {
      ud[i].c11.re=ud[i].c22.re=cos(ad[i]);
      ud[i].c12.re=sin(ad[i]);
      ud[i].c21.re=-ud[i].c12.re;
      ud[i].c33.re=1.0;
   }
   set_flags(UPDATED_UD);
   
   fwd_su3_rk2(n,eps);
   fwd_u1_rk2(n,eps*qel2);
   
   dmax=0.0;
   for(i=0;i<4*VOLUME;i++)
   {
      d=fabs(ud[i].c11.re-cos(ad[i]));
      if(d>dmax) dmax=d;
      d=fabs(ud[i].c12.re-sin(ad[i]));
      if(d>dmax) dmax=d;
   }
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between fwd_su3_rk2() and fwd_u1_rk2(), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   
   random_ad();
   ad=adfld();
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for(i=0;i<4*VOLUME;i++)
   {
      ud[i].c11.re=ud[i].c22.re=cos(ad[i]);
      ud[i].c12.re=sin(ad[i]);
      ud[i].c21.re=-ud[i].c12.re;
      ud[i].c33.re=1.0;
   }
   set_flags(UPDATED_UD);
   
   fwd_su3_rk3(n,eps);
   fwd_u1_rk3(n,eps*qel2);
   
   dmax=0.0;
   for(i=0;i<4*VOLUME;i++)
   {
      d=fabs(ud[i].c11.re-cos(ad[i]));
      if(d>dmax) dmax=d;
      d=fabs(ud[i].c12.re-sin(ad[i]));
      if(d>dmax) dmax=d;
   }
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between fwd_su3_rk3() and fwd_u1_rk3(), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
