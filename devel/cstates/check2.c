
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge covariance of mul_cfactor().
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
#include "cstates.h"
#include "global.h"
#include "gflds_utils.h"
#include "sflds_utils.h"

static int lsize[4]={L0,L1,L2,L3};
static int faces[4]={FACE0,FACE1,FACE2,FACE3};

static int q,init=0,*ix0b;
static double *gmirror;


static void set_gmirror(int mu)
{
   int i,j,ix,iy,iz;
   int x[4],xm[4];
   int tag,mirror_rank;
   double *sbuf,*rbuf;
   double* g;
   MPI_Status stat;

   if (init==0)
   {
      gmirror=amalloc(BNDRY*sizeof(*gmirror),ALIGN);
      error(gmirror==NULL,1,"set_lines [check2.c]",
            "Unable to allocate service array");

      ix0b=amalloc(VOLUME*sizeof(*ix0b),ALIGN);
      error(ix0b==NULL,1,"set_lines [check2.c]",
            "Unable to allocate service array");
   }

   g=g1tr();

   for(i=0;i<4;++i)
      xm[i]=lsize[i];
   xm[mu]=1;

   i=0;
   for (x[0]=0;x[0]<xm[0];++x[0])
   for (x[1]=0;x[1]<xm[1];++x[1])
   for (x[2]=0;x[2]<xm[2];++x[2])
   for (x[3]=0;x[3]<xm[3];++x[3])
   {
      ix=ipt[x[3]+L3*x[2]+L2*L3*x[1]+L1*L2*L3*x[0]];

      ix0b[ix]=i;
      gmirror[i]=g[ix];

      iz=ix;
      for(j=1;j<lsize[mu];++j)
      {
         iy=iup[iz][mu];
         ix0b[iy]=i;
         iz=iy;
      }
      ++i;
   }

   mirror_rank=get_mirror_rank();

   sbuf=gmirror;
   rbuf=gmirror+faces[mu];

   tag=mpi_tag();
   if (cpr[1]<NPROC1/2)
   {
      MPI_Send(sbuf,faces[mu],MPI_DOUBLE,mirror_rank,tag,MPI_COMM_WORLD);
      MPI_Recv(rbuf,faces[mu],MPI_DOUBLE,mirror_rank,tag,MPI_COMM_WORLD,&stat);
   }
   else
   {
      MPI_Recv(rbuf,faces[mu],MPI_DOUBLE,mirror_rank,tag,MPI_COMM_WORLD,&stat);
      MPI_Send(sbuf,faces[mu],MPI_DOUBLE,mirror_rank,tag,MPI_COMM_WORLD);
   }
}


void transform_gfactor(int mu,spinor_dble *pk,spinor_dble *pl)
{
   int ix,i;
   double arg,*g;
   complex_dble zx,sc,*prc,*psc;

   g=g1tr();

   for (ix=0;ix<VOLUME;ix++)
   {
      arg=2.0*g[ix]-gmirror[ix0b[ix]]-gmirror[faces[mu]+ix0b[ix]];

      arg*=0.5*q;

      zx.re=cos(arg);
      zx.im=sin(arg);

      #define cmpl_mul(a,b,c) \
      (a).re=(b).re*(c).re-(b).im*(c).im; \
      (a).im=(b).re*(c).im+(b).im*(c).re

      psc=(complex_dble*)(pk+ix);
      prc=(complex_dble*)(pl+ix);
      for(i=0;i<12;i++)
      {
         sc=psc[i];
         cmpl_mul(prc[i],zx,sc);
      }

      #undef cmpl_mul
   }
}


int main(int argc,char *argv[])
{
   int nu,my_rank,bc,cs;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   double d;
   spinor_dble **psd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      printf("\n");
      printf("Gauge covariance of mul_cfactor()\n");
      printf("---------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
      error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                  "Syntax: check2 -cs <cstar> [-bc <type>] [-q <echarge>]");

      cs=find_opt(argc,argv,"-cs");
      error_root(cs==0 || sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check2.c]",
                  "Syntax: check2 -cs <cstar> [-bc <type>] [-q <echarge>]");

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check2.c]",
         "Syntax: check2 -cs <cstar> [-bc <type>] [-q <echarge>]");
      }
      else
         q=-8;
   }

   set_flds_parms(2,0);
   print_flds_parms();

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
   alloc_wsd(4);
   psd=reserve_wsd(4);

   MPI_Bcast(&q,1,MPI_INT,0,MPI_COMM_WORLD);
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   set_dirac_parms9(q,-0.0123,su3csw,u1csw,cF[0],cF[1],
   theta[0],theta[1],theta[2]);
   print_dirac_parms();


   for (nu=1;nu<=bc_cstar();++nu)
   {
      random_g();
      random_gflds();

      random_sd(NSPIN,psd[0],1.0);

      set_gmirror(nu);
      mul_cfactor(0,1,nu,psd[0],psd[1]);
      transform_gflds();
      mul_cfactor(0,1,nu,psd[0],psd[2]);
      transform_gfactor(nu,psd[2],psd[3]);

      mulr_spinor_add_dble(VOLUME,psd[3],psd[1],-1.0);
      d=norm_square_dble(VOLUME,1,psd[3])/norm_square_dble(VOLUME,1,psd[0]);

      if (my_rank==0)
      {   
         printf("mu = %2d\n",nu);
         printf("Normalized difference = %.2e\n",sqrt(d));
         printf("(should be around 1*10^(-15) or so)\n\n");
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
