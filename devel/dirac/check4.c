
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge covariance of Dw_dble().
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
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"
#include "gflds_utils.h"
#include "sflds_utils.h"



int main(int argc,char *argv[])
{
   int my_rank,i,bc,cf,q;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   double mu,d;
   spinor_dble **psd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      printf("\n");
      printf("Gauge covariance of Dw_dble() (random fields)\n");
      printf("---------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check4.c]",
                  "Syntax: check4 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check4.c]",
                  "Syntax: check4 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      }
      else
         q=-3;
   }

   MPI_Bcast(&cf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&q,1,MPI_INT,0,MPI_COMM_WORLD);
   set_flds_parms(cf,0);
   print_flds_parms();
   if(gauge()==1) q=0;

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   if (gauge()==2)
   {
      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
   }   
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();
   alloc_wsd(5);
   psd=reserve_wsd(5);

   mu=0.0376;
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   if ((gauge()&1)!=0) su3csw=0.95;
   if ((gauge()&2)!=0) u1csw=0.8;
   if (bc_type()!=3)
   {
      cF[0]=1.301;
      cF[1]=0.789;
   }
   if (bc_cstar()==0)
   {
      theta[0]=0.35;
      theta[1]=-1.25;
      theta[2]=0.78;
   }
   set_dirac_parms9(q,-0.0123,su3csw,u1csw,cF[0],cF[1],
                    theta[0],theta[1],theta[2]);
   print_dirac_parms();

   random_g();
   random_gflds();
   sw_term(NO_PTS);

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[4]);
   bnd_sd2zero(ALL_PTS,psd[4]);
   Dw_dble(mu,psd[0],psd[1]);
   mulr_spinor_add_dble(VOLUME,psd[4],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[4]);
   error(d!=0.0,1,"main [check4.c]","Dw_dble() changes the input field");

   Dw_dble(mu,psd[0],psd[4]);
   mulr_spinor_add_dble(VOLUME,psd[4],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[4]);
   error(d!=0.0,1,"main [check4.c]","Action of Dw_dble() depends "
         "on the boundary values of the input field");

   assign_sd2sd(VOLUME,psd[1],psd[4]);
   bnd_sd2zero(ALL_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME,psd[4],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[4]);
   error(d!=0.0,1,"main [check4.c]",
         "Dw_dble() does not vanish at global time 0 and NPROC0*L0-1 ");

   transform_sd(psd[0],psd[2]);
   transform_gflds();
   sw_term(NO_PTS);
   Dw_dble(mu,psd[2],psd[3]);
   transform_sd(psd[1],psd[2]);

   mulr_spinor_add_dble(VOLUME,psd[3],psd[2],-1.0);
   d=norm_square_dble(VOLUME,1,psd[3])/norm_square_dble(VOLUME,1,psd[0]);

   if (my_rank==0)
   {
      printf("Normalized difference = %.2e\n",sqrt(d));
      printf("(should be around 1*10^(-15) or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
