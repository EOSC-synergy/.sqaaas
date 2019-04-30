
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2005, 2008, 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Hermiticity of Dw_dble() and comparison with Dwee_dble(),..,Dwhat_dble().
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
#include "utils.h"
#include "u1flds.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"
#include "gflds_utils.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,i,cf,q;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   double mu,d;
   complex_dble z1,z2;
   spinor_dble **psd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);
      printf("\n");
      printf("Hermiticity of Dw_dble() and comparison with Dwee_dble(),..,"
             "Dwhat_dble()\n");
      printf("------------------------------------------------------------"
             "------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      printf("For this test to pass, the calculated differences\n");
      printf("should be at most 1*10^(-15) or so\n\n");

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check6.c]",
                  "Syntax: check6 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check6.c]",
                  "Syntax: check6 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
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

   random_gflds();
   sw_term(NO_PTS);

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   Dw_dble(mu,psd[0],psd[2]);
   mulg5_dble(VOLUME,psd[2]);
   Dw_dble(-mu,psd[1],psd[3]);
   mulg5_dble(VOLUME,psd[3]);

   z1=spinor_prod_dble(VOLUME,1,psd[0],psd[3]);
   z2=spinor_prod_dble(VOLUME,1,psd[2],psd[1]);

   d=sqrt((z1.re-z2.re)*(z1.re-z2.re)+
          (z1.im-z2.im)*(z1.im-z2.im));
   d/=sqrt((double)(12*NPROC)*(double)(VOLUME));

   if (my_rank==0)
      printf("Deviation from gamma5-Hermiticity             = %.1e\n",d);

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[1]);
   assign_sd2sd(VOLUME,psd[2],psd[3]);
   Dwee_dble(mu,psd[1],psd[2]);

   bnd_sd2zero(EVEN_PTS,psd[0]);
   mulr_spinor_add_dble(VOLUME,psd[1],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[1]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwee_dble() changes the input field in unexpected ways");

   mulr_spinor_add_dble(VOLUME/2,psd[2]+(VOLUME/2),psd[3]+(VOLUME/2),-1.0);
   assign_sd2sd(VOLUME/2,psd[2],psd[4]);
   bnd_sd2zero(EVEN_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME/2,psd[2],psd[4],-1.0);
   d=norm_square_dble(VOLUME,1,psd[2]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwee_dble() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[1]);
   assign_sd2sd(VOLUME,psd[2],psd[3]);
   Dwoo_dble(mu,psd[1],psd[2]);

   bnd_sd2zero(ODD_PTS,psd[0]);
   mulr_spinor_add_dble(VOLUME,psd[1],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[1]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwoo_dble() changes the input field in unexpected ways");

   mulr_spinor_add_dble(VOLUME/2,psd[2],psd[3],-1.0);
   assign_sd2sd(VOLUME/2,psd[2]+(VOLUME/2),psd[4]+(VOLUME/2));
   bnd_sd2zero(ODD_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME/2,psd[2]+(VOLUME/2),psd[4]+(VOLUME/2),-1.0);
   d=norm_square_dble(VOLUME,1,psd[2]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwoo_dble() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[1]);
   assign_sd2sd(VOLUME,psd[2],psd[3]);
   Dwoe_dble(psd[1],psd[2]);

   bnd_sd2zero(EVEN_PTS,psd[0]);
   mulr_spinor_add_dble(VOLUME,psd[1],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[1]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwoe_dble() changes the input field in unexpected ways");

   mulr_spinor_add_dble(VOLUME/2,psd[2],psd[3],-1.0);
   assign_sd2sd(VOLUME/2,psd[2]+(VOLUME/2),psd[4]+(VOLUME/2));
   bnd_sd2zero(ODD_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME/2,psd[2]+(VOLUME/2),psd[4]+(VOLUME/2),-1.0);
   d=norm_square_dble(VOLUME,1,psd[2]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwoe_dble() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[1]);
   assign_sd2sd(VOLUME,psd[2],psd[3]);
   Dweo_dble(psd[1],psd[2]);

   bnd_sd2zero(ODD_PTS,psd[0]);
   mulr_spinor_add_dble(VOLUME,psd[1],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[1]);

   error(d!=0.0,1,"main [check6.c]",
         "Dweo_dble() changes the input field in unexpected ways");

   mulr_spinor_add_dble(VOLUME/2,psd[2]+(VOLUME/2),psd[3]+(VOLUME/2),-1.0);
   assign_sd2sd(VOLUME/2,psd[2],psd[4]);
   bnd_sd2zero(EVEN_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME/2,psd[2],psd[4],-1.0);
   d=norm_square_dble(VOLUME,1,psd[2]);

   error(d!=0.0,1,"main [check6.c]",
         "Dweo_dble() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[1]);
   assign_sd2sd(VOLUME,psd[2],psd[3]);
   Dwhat_dble(mu,psd[1],psd[2]);

   bnd_sd2zero(EVEN_PTS,psd[0]);
   mulr_spinor_add_dble(VOLUME,psd[1],psd[0],-1.0);
   d=norm_square_dble(VOLUME,1,psd[1]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwhat_dble() changes the input field in unexpected ways");

   mulr_spinor_add_dble(VOLUME/2,psd[2]+(VOLUME/2),psd[3]+(VOLUME/2),-1.0);
   assign_sd2sd(VOLUME/2,psd[2],psd[4]);
   bnd_sd2zero(EVEN_PTS,psd[4]);
   mulr_spinor_add_dble(VOLUME/2,psd[2],psd[4],-1.0);
   d=norm_square_dble(VOLUME,1,psd[2]);

   error(d!=0.0,1,"main [check6.c]",
         "Dwhat_dble() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[2]);
   Dw_dble(mu,psd[0],psd[1]);
   Dwee_dble(mu,psd[2],psd[3]);
   set_sd2zero(VOLUME/2,psd[0]);
   mulr_spinor_add_dble(VOLUME/2,psd[0],psd[3],-1.0);
   Dweo_dble(psd[2],psd[0]);
   set_sd2zero(VOLUME/2,psd[3]);
   mulr_spinor_add_dble(VOLUME/2,psd[3],psd[0],-1.0);

   Dwoo_dble(mu,psd[2],psd[3]);
   Dwoe_dble(psd[2],psd[4]);
   mulr_spinor_add_dble(VOLUME/2,psd[3]+(VOLUME/2),psd[4]+(VOLUME/2),1.0);

   mulr_spinor_add_dble(VOLUME,psd[3],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[3])/norm_square_dble(VOLUME,1,psd[1]);
   d=sqrt(d);

   if (my_rank==0)
      printf("Deviation of Dw_dble() from Dwee_dble(),..    = %.1e\n",d);

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(NSPIN,psd[0],psd[1]);
   Dwhat_dble(mu,psd[0],psd[2]);

   Dwoe_dble(psd[1],psd[1]);
   Dwee_dble(mu,psd[1],psd[1]);
   Dwoo_dble(0.0,psd[1],psd[1]);
   Dweo_dble(psd[1],psd[1]);

   mulr_spinor_add_dble(VOLUME/2,psd[1],psd[2],-1.0);
   d=norm_square_dble(VOLUME/2,1,psd[1])/norm_square_dble(VOLUME/2,1,psd[2]);
   d=sqrt(d);

   if (my_rank==0)
      printf("Deviation of Dwhat_dble() from Dwee_dble(),.. = %.1e\n",d);

   for (i=0;i<4;i++)
      random_sd(NSPIN,psd[i],1.0);

   assign_sd2sd(VOLUME,psd[0],psd[2]);

   set_tm_parms(1);
   Dw_dble(mu,psd[0],psd[1]);
   set_tm_parms(0);

   Dwee_dble(mu,psd[2],psd[3]);
   mulr_spinor_add_dble(VOLUME/2,psd[1],psd[3],-1.0);
   Dweo_dble(psd[2],psd[1]);
   Dwoe_dble(psd[2],psd[3]);
   mulr_spinor_add_dble(VOLUME/2,psd[1]+(VOLUME/2),psd[3]+(VOLUME/2),-1.0);
   Dwoo_dble(0.0,psd[2],psd[3]);
   mulr_spinor_add_dble(VOLUME/2,psd[1]+(VOLUME/2),psd[3]+(VOLUME/2),-1.0);
   d=norm_square_dble(VOLUME,1,psd[1])/norm_square_dble(VOLUME,1,psd[2]);
   d=sqrt(d);

   if (my_rank==0)
   {
      printf("Check of Dw_dble()|eoflg=1                    = %.1e\n\n",d);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
