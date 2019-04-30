
/*******************************************************************************
*
* File check14.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency between compact U(1) and SU(3) gauge forces and actions.
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
#include "forces.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)



int main(int argc,char *argv[])
{
   int my_rank,bc,sf,cs,i;
   double alpha,beta,qel2,phi[2],phi_prime[2];
   double *ad;
   su3_dble *ud;
   double u1act,u3act;
   double d,dmax,dmax_all;
   double *frc1;
   su3_alg_dble *frc3;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check14.log","w",stdout);

      printf("\n");
      printf("Consistency between compact U(1) and SU(3) gauge forces and actions\n");
      printf("-------------------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check14.c]",
                    "Syntax: check14 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");

      sf=find_opt(argc,argv,"-sf");

      if (sf!=0)
         error_root(sscanf(argv[sf+1],"%d",&sf)!=1,1,"main [check14.c]",
                    "Syntax: check14 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check14.c]",
                    "Syntax: check14 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");
   }
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&sf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,sf,cs,phi,phi_prime);
   print_bc_parms();

   alpha=.013;
   qel2=2.45;
   beta=1./(16.*atan(1)*qel2*alpha);
   set_u1lat_parms(0,alpha,1.0/sqrt(qel2),0.0,0.482,0.87,0.57);
   set_su3lat_parms(beta,0.482,0.87,0.57);
   print_lat_parms();

   start_ranlux(0,12345);
   geometry();
   
   /*
   cm3x3_zero(1,M);
   d=.21;
   M[0].c12.re=d;
   M[0].c21.re=-M[0].c12.re;
   cm3x3_unity(1,M+1);
   prod2su3alg(M,M+1,&X);
   if (my_rank==0)
   {
      printf("T3 = ( %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e )\n\n",
      X.c1/d,X.c2/d,X.c3/d,X.c4/d,X.c5/d,X.c6/d,X.c7/d,X.c8/d);
   }
   
   expXsu3(1.0,&X,M+1);
   if (my_rank==0)
   {
      printf("                 | ( %.2e  %.2e ) ( %.2e  %.2e ) ( %.2e  %.2e ) |\n",
             M[1].c11.re,M[1].c11.im,M[1].c12.re,M[1].c12.im,M[1].c13.re,M[1].c13.im);
      printf("exp ( d * T3 ) = | ( %.2e  %.2e ) ( %.2e  %.2e ) ( %.2e  %.2e ) |\n",
             M[1].c21.re,M[1].c21.im,M[1].c22.re,M[1].c22.im,M[1].c23.re,M[1].c23.im);
      printf("                 | ( %.2e  %.2e ) ( %.2e  %.2e ) ( %.2e  %.2e ) |\n\n",
             M[1].c31.re,M[1].c31.im,M[1].c32.re,M[1].c32.im,M[1].c33.re,M[1].c33.im);
      
      printf("cos(d) = %.2e\n",cos(d));
      printf("sin(d) = %.2e\n\n",sin(d));
   }
   */

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
   
   
   /****************************************************************************
   
   A = sum_{a=1}^8 Aa Ta
   (A,B) = -2 tr AB = 12 (A1 B1 + A2 B2) - 6 (A1 B2 + A2 B1)
                      + 4 (A3 B3 + A4 B4 + A5 B5 + A6 B6 + A7 B7 + A8 B8)
   
   d/dt f(e^{t*X} U) |_{t=0} = ( X , d f(U) )
   d f(U) = sum_{a=1}^8 Ta da f(U)
   
   in particular
   (T3,T3) = 4
   d/dt f(e^{t*T3} U) |_{t=0} = ( T3 , d f(U) ) = 4 d3 f(U)
   

        |  0 1 0 |
   T3 = | -1 0 0 |
        |  0 0 0 |
   
                   | cos A   sin A   0 |
   U = exp(A*T3) = | -sin A  cos A   0 |
                   | 0       0       1 |
   
   tr U = 1 + 2 cos A
   Re tr U = 1 + 2 Re exp(iA)
   
   S1(A) = 1 - Re exp(i Ap)
   S3(U) = 1 - 1/3 Re tr Up
   
   S3(exp(A*X3)) = 1 - 1/3 Re tr exp(A*X3)
                 = 1 - 1/3 ( 1 + 2 Re exp(i Ap) )
                 = 2/3 S1(A)
   
   d^3 S3(U) = d/dt S3(exp(t*T3) U) |_{t=0}
   4 d^3 S3(exp(A*T3)) = d/dt S3(exp(t*T3) exp(A*T3)) |_{t=0}
                       = d/dA [S3(exp(A*T3))]
                       = 2/3 d/dA S1(A)
   
   beta_QCD = 6/g2
   beta_QED = 1/(qel2 e2)
   
   we have chosen beta_QCD = beta_QED i.e. g2/6 = qel2 e2
   
   plaq_su3frc = g2 F3 = g2/6 F1 = qel2 e2 F1 = qel2 plaq_u1frc
   
   ****************************************************************************/
   
   u1act=action6(0);
   u3act=action0(0);
   dmax=fabs(2.*u1act-3.*u3act);
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between action0(0) and action6(0), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   
   force0(3.92);
   force6(3.92);
   frc1=mdflds()->u1frc;
   frc3=mdflds()->su3frc;
   dmax=0.0;
   for(i=0;i<4*VOLUME;i++)
   {
      d=fabs(frc1[i]-6.*frc3[i].c3);
      if(d>dmax) dmax=d;
   }
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between force0() and force6(), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   plaq_su3frc();
   plaq_u1frc();
   frc1=mdflds()->u1frc;
   frc3=mdflds()->su3frc;
   dmax=0.0;
   for(i=0;i<4*VOLUME;i++)
   {
      d=fabs(qel2*frc1[i]-frc3[i].c3);
      if(d>dmax) dmax=d;
   }
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between plaq_su3frc() and plaq_u1frc(), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   
   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
