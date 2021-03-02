
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2008, 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Hermiticity of Dw() and comparison with Dwee(),..,Dwhat().
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
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"
#include "gflds_utils.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,i,cf,q;
   float mu,d;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   complex z1,z2;
   spinor **ps;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      printf("\n");
      printf("Hermiticity of Dw() and comparison with Dwee(),..,Dwhat()\n");
      printf("---------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      printf("For this test to pass, the calculated differences\n");
      printf("should be at most 1*10^(-5) or so\n\n");

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>] [-q <echarge>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check3.c]",
                  "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check3.c]",
                  "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>] [-q <echarge>]");
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
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();
   alloc_ws(5);
   ps=reserve_ws(5);

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
   assign_swd2sw();

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   Dw(mu,ps[0],ps[2]);
   mulg5(VOLUME,ps[2]);
   Dw(-mu,ps[1],ps[3]);
   mulg5(VOLUME,ps[3]);

   z1=spinor_prod(VOLUME,1,ps[0],ps[3]);
   z2=spinor_prod(VOLUME,1,ps[2],ps[1]);

   d=(float)(sqrt((double)((z1.re-z2.re)*(z1.re-z2.re)+
                           (z1.im-z2.im)*(z1.im-z2.im))));
   d/=(float)(sqrt((double)(12*NPROC)*(double)(VOLUME)));

   if (my_rank==0)
      printf("Deviation from gamma5-Hermiticity    = %.1e\n",d);

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_s2s(VOLUME,ps[2],ps[3]);
   Dwee(mu,ps[1],ps[2]);

   bnd_s2zero(EVEN_PTS,ps[0]);
   mulr_spinor_add(VOLUME,ps[1],ps[0],-1.0f);
   d=norm_square(VOLUME,1,ps[1]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwee() changes the input field in unexpected ways");

   mulr_spinor_add(VOLUME/2,ps[2]+(VOLUME/2),ps[3]+(VOLUME/2),-1.0f);
   assign_s2s(VOLUME/2,ps[2],ps[4]);
   bnd_s2zero(EVEN_PTS,ps[4]);
   mulr_spinor_add(VOLUME/2,ps[2],ps[4],-1.0f);
   d=norm_square(VOLUME,1,ps[2]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwee() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_s2s(VOLUME,ps[2],ps[3]);
   Dwoo(mu,ps[1],ps[2]);

   bnd_s2zero(ODD_PTS,ps[0]);
   mulr_spinor_add(VOLUME,ps[1],ps[0],-1.0f);
   d=norm_square(VOLUME,1,ps[1]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwoo() changes the input field in unexpected ways");

   mulr_spinor_add(VOLUME/2,ps[2],ps[3],-1.0f);
   assign_s2s(VOLUME/2,ps[2]+(VOLUME/2),ps[4]+(VOLUME/2));
   bnd_s2zero(ODD_PTS,ps[4]);
   mulr_spinor_add(VOLUME/2,ps[2]+(VOLUME/2),ps[4]+(VOLUME/2),-1.0f);
   d=norm_square(VOLUME,1,ps[2]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwoo() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_s2s(VOLUME,ps[2],ps[3]);
   Dwoe(ps[1],ps[2]);

   bnd_s2zero(EVEN_PTS,ps[0]);
   mulr_spinor_add(VOLUME,ps[1],ps[0],-1.0f);
   d=norm_square(VOLUME,1,ps[1]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwoe() changes the input field in unexpected ways");

   mulr_spinor_add(VOLUME/2,ps[2],ps[3],-1.0f);
   assign_s2s(VOLUME/2,ps[2]+(VOLUME/2),ps[4]+(VOLUME/2));
   bnd_s2zero(ODD_PTS,ps[4]);
   mulr_spinor_add(VOLUME/2,ps[2]+(VOLUME/2),ps[4]+(VOLUME/2),-1.0f);
   d=norm_square(VOLUME,1,ps[2]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwoe() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_s2s(VOLUME,ps[2],ps[3]);
   Dweo(ps[1],ps[2]);

   bnd_s2zero(ODD_PTS,ps[0]);
   mulr_spinor_add(VOLUME,ps[1],ps[0],-1.0f);
   d=norm_square(VOLUME,1,ps[1]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dweo() changes the input field in unexpected ways");

   mulr_spinor_add(VOLUME/2,ps[2]+(VOLUME/2),ps[3]+(VOLUME/2),-1.0f);
   assign_s2s(VOLUME/2,ps[2],ps[4]);
   bnd_s2zero(EVEN_PTS,ps[4]);
   mulr_spinor_add(VOLUME/2,ps[2],ps[4],-1.0f);
   d=norm_square(VOLUME,1,ps[2]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dweo() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[1]);
   assign_s2s(VOLUME,ps[2],ps[3]);
   Dwhat(mu,ps[1],ps[2]);

   bnd_s2zero(EVEN_PTS,ps[0]);
   mulr_spinor_add(VOLUME,ps[1],ps[0],-1.0f);
   d=norm_square(VOLUME,1,ps[1]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwhat() changes the input field in unexpected ways");

   mulr_spinor_add(VOLUME/2,ps[2]+(VOLUME/2),ps[3]+(VOLUME/2),-1.0f);
   assign_s2s(VOLUME/2,ps[2],ps[4]);
   bnd_s2zero(EVEN_PTS,ps[4]);
   mulr_spinor_add(VOLUME/2,ps[2],ps[4],-1.0f);
   d=norm_square(VOLUME,1,ps[2]);

   error(d!=0.0f,1,"main [check3.c]",
         "Dwhat() changes the output field where it should not");

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[2]);
   Dw(mu,ps[0],ps[1]);
   Dwee(mu,ps[2],ps[3]);
   set_s2zero(VOLUME/2,ps[0]);
   mulr_spinor_add(VOLUME/2,ps[0],ps[3],-1.0f);
   Dweo(ps[2],ps[0]);
   set_s2zero(VOLUME/2,ps[3]);
   mulr_spinor_add(VOLUME/2,ps[3],ps[0],-1.0f);

   Dwoo(mu,ps[2],ps[3]);
   Dwoe(ps[2],ps[4]);
   mulr_spinor_add(VOLUME/2,ps[3]+(VOLUME/2),ps[4]+(VOLUME/2),1.0f);
   mulr_spinor_add(VOLUME,ps[3],ps[1],-1.0f);
   d=norm_square(VOLUME,1,ps[3])/norm_square(VOLUME,1,ps[1]);
   d=(float)(sqrt((double)(d)));

   if (my_rank==0)
      printf("Deviation of Dw() from Dwee(),..     = %.1e\n",d);

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(NSPIN,ps[0],ps[1]);
   Dwhat(mu,ps[0],ps[2]);

   Dwoe(ps[1],ps[1]);
   Dwee(mu,ps[1],ps[1]);
   Dwoo(0.0,ps[1],ps[1]);
   Dweo(ps[1],ps[1]);

   mulr_spinor_add(VOLUME/2,ps[1],ps[2],-1.0f);
   d=norm_square(VOLUME/2,1,ps[1])/norm_square(VOLUME/2,1,ps[2]);
   d=(float)(sqrt((double)(d)));

   if (my_rank==0)
      printf("Deviation of Dwhat() from Dwee(),..  = %.1e\n",d);

   for (i=0;i<4;i++)
      random_s(NSPIN,ps[i],1.0f);

   assign_s2s(VOLUME,ps[0],ps[2]);

   set_tm_parms(1);
   Dw(mu,ps[0],ps[1]);
   set_tm_parms(0);

   Dwee(mu,ps[2],ps[3]);
   mulr_spinor_add(VOLUME/2,ps[1],ps[3],-1.0f);
   Dweo(ps[2],ps[1]);
   Dwoe(ps[2],ps[3]);
   mulr_spinor_add(VOLUME/2,ps[1]+(VOLUME/2),ps[3]+(VOLUME/2),-1.0f);
   Dwoo(0.0f,ps[2],ps[3]);
   mulr_spinor_add(VOLUME/2,ps[1]+(VOLUME/2),ps[3]+(VOLUME/2),-1.0f);
   d=norm_square(VOLUME,1,ps[1])/norm_square(VOLUME,1,ps[2]);
   d=(float)(sqrt((double)(d)));

   if (my_rank==0)
   {
      printf("Check of Dw()|eoflg=1                = %.1e\n\n",d);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
