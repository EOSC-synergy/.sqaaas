
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2007, 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Direct check of Aw_dble() and Aw().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "sflds.h"
#include "vflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "dfl.h"
#include "little.h"
#include "global.h"
#include "gflds_utils.h"


static void random_basis(int Ns)
{
   int i;
   spinor **ws;

   ws=reserve_ws(Ns);

   for (i=0;i<Ns;i++)
   {
      random_s(VOLUME,ws[i],1.0f);
      bnd_s2zero(ALL_PTS,ws[i]);
   }

   dfl_subspace(ws);
   release_ws();
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cf,q;
   int bs[4],Ns,nb,nv;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2];
   double mu,dev;
   complex **wv,z;
   complex_dble **wvd,zd;
   spinor **ws;
   spinor_dble **wsd;
   FILE *fin=NULL,*flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Direct check of Aw_dble() and Aw()\n");
      printf("----------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      read_line("Ns","%d",&Ns);
      fclose(fin);

      printf("bs = %d %d %d %d\n",bs[0],bs[1],bs[2],bs[3]);
      printf("Ns = %d\n\n",Ns);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-gg <gauge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check3.c]",
                  "Syntax: check3 [-bc <type>] [-gg <gauge>]");
      else
         cf=1;
   }

   MPI_Bcast(&cf,1,MPI_INT,0,MPI_COMM_WORLD);
   set_flds_parms(cf,0);
   print_flds_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
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

   q=3;
   if(gauge()==1) q=0;
   mu=0.0376;
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   if ((gauge()&1)!=0) su3csw=0.95;
   if ((gauge()&2)!=0) u1csw=0.8;
   if (bc_type()!=3)
   {
      cF[0]=1.301;
      cF[1]=0.789;
   }
   set_dirac_parms9(q,-0.0123,su3csw,u1csw,cF[0],cF[1],0.0,0.0,0.0);
   print_dirac_parms();

   set_dfl_parms(bs,Ns);

   start_ranlux(0,123456);
   geometry();

   alloc_ws(Ns+2);
   alloc_wsd(2);
   alloc_wv(3);
   alloc_wvd(3);

   ws=reserve_ws(2);
   wsd=reserve_wsd(2);
   wv=reserve_wv(3);
   wvd=reserve_wvd(3);
   nb=VOLUME/(bs[0]*bs[1]*bs[2]*bs[3]);
   nv=Ns*nb;

   random_gflds();
   random_basis(Ns);
   set_Aw(mu);
   sw_term(NO_PTS);
   assign_swd2sw();

   random_vd(nv,wvd[0],1.0);
   Aw_dble(wvd[0],wvd[1]);
   dfl_vd2sd(wvd[0],wsd[0]);
   Dw_dble(mu,wsd[0],wsd[1]);
   dfl_sd2vd(wsd[1],wvd[2]);

   zd.re=-1.0;
   zd.im=0.0;
   mulc_vadd_dble(nv,wvd[2],wvd[1],zd);
   dev=vnorm_square_dble(nv,1,wvd[2])/vnorm_square_dble(nv,1,wvd[1]);

   if (my_rank==0)
      printf("Relative deviation (Aw_dble) = %.1e\n",sqrt(dev));

   random_v(nv,wv[0],1.0f);
   Aw(wv[0],wv[1]);
   dfl_v2s(wv[0],ws[0]);
   Dw((float)(mu),ws[0],ws[1]);
   dfl_s2v(ws[1],wv[2]);

   z.re=-1.0f;
   z.im=0.0f;
   mulc_vadd(nv,wv[2],wv[1],z);
   dev=(double)(vnorm_square(nv,1,wv[2])/vnorm_square(nv,1,wv[1]));

   if (my_rank==0)
   {
      printf("Relative deviation (Aw)      = %.1e\n\n",sqrt(dev));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
