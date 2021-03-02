
/*******************************************************************************
*
* File check2b.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge covariance with orbifold constraint of mul_cfactor().
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


static void cstar_double(int vol,double* ad)
{
   double *am;

   am=ad+vol;
   for (;ad<am;ad++)  (*ad) *= -1.0;
}


static void orbi_cpy_gtr(void)
{
   int mirror, tag;
   MPI_Status stat;
   double *g;

   if(bc_cstar()>0) {
      g=g1tr();
      mirror=get_mirror_rank();
      tag=mpi_tag();
      if(cpr[1]<NPROC1/2) {
         MPI_Send(g,NSPIN,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD);
      } else {
         MPI_Recv(g,NSPIN,MPI_DOUBLE,mirror,tag,MPI_COMM_WORLD,&stat);
         cstar_double(NSPIN,g);
      }
   }
}


static void pack_and_send_g1buf(void)
{
   static double *gbuf=NULL;
   static int nfc[8],ofs[8],np[8];
   int ifc,ib,ix;
   int saddr,raddr;
   int nbf,tag;
   double *sbuf,*rbuf,*g;
   MPI_Status stat;

   g=g1tr();

   if (BNDRY!=0 && gbuf==NULL)
   {
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);
      error(gbuf==NULL,1,"pack_and_send_gbuf [check2b]",
            "Unable to allocate memory space for gbuf");
   }

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (ifc=0;ifc<8;ifc++)
   {
      for (ib=0;ib<nfc[ifc];ib++)
      {
         ix=map[ofs[ifc]+ib];
         gbuf[ofs[ifc]+ib]=g[ix];
      }
   }

   np[0]=np[2]=np[4]=np[6]=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   if ((bc_cstar()>=2)&&(cpr[1]>=NPROC1/2))
      np[4]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+NPROC2+cpr[3])&0x1;
   if ((bc_cstar()>=3)&&(cpr[1]>=NPROC1/2))
      np[6]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+cpr[3]+NPROC3)&0x1;
   for(ifc=0;ifc<4;++ifc)
      np[2*ifc+1]=np[2*ifc];

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=nfc[ifc];

      if (nbf>0)
      {
      tag=mpi_tag();
         saddr=npr[ifc^0x1];
         raddr=npr[ifc];
         sbuf=gbuf+ofs[ifc];
         rbuf=g+VOLUME+ofs[ifc];

         if (np[ifc]==0)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


static void cstar_random_g1(void)
{
   int ix,t,bc;
   double *g1x;
   double twopi;

   bc=bc_type();

   g1x=g1tr();
   twopi=4*atan(1.);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t>0)||(bc!=1))
      {
         ranlxd(g1x,1);
         (*g1x)=((*g1x)-.5)*twopi;
      }
      else
         (*g1x)=0.;

      g1x+=1;
   }

   orbi_cpy_gtr();

   if (BNDRY>0)
      pack_and_send_g1buf();
}


int main(int argc,char *argv[])
{
   int nu,my_rank,bc,cs,q;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   double d;
   spinor_dble **psd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2b.log","w",stdout);
      printf("\n");
      printf("Gauge covariance with orbifold constraint of mul_cfactor()\n");
      printf("----------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2b.c]",
                     "Syntax: check2b -cs <cstar> [-bc <type>] [-q <echarge>]");

      cs=find_opt(argc,argv,"-cs");
      error_root(cs==0 || sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check2b.c]",
                  "Syntax: check2b -cs <cstar> [-bc <type>] [-q <echarge>]");

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check2b.c]",
                     "Syntax: check2b -cs <cstar> [-bc <type>] [-q <echarge>]");
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
      random_gflds();
      orbi_cpy_ad();

      cstar_random_g1();

      random_sd(NSPIN,psd[0],1.0);
      mul_cfactor(0,1,nu,psd[0],psd[1]);

      transform_gflds();
      mul_cfactor(0,1,nu,psd[0],psd[2]);
      transform_sd(psd[2],psd[3]);

      mulr_spinor_add_dble(VOLUME,psd[3],psd[1],-1.0);
      d=norm_square_dble(VOLUME,1,psd[3])/norm_square_dble(VOLUME,1,psd[0]);

      if (my_rank==0)
      {   
         printf("mu = %2d\n",nu);
         printf("Normalized difference = %.2e\n",sqrt(d));
         printf("(should be around 1*10^(-15) or so)\n\n");
      }
   }

   random_gflds();
   orbi_cpy_ad();

   cstar_random_g1();

   random_sd(NSPIN,psd[0],1.0);
   mul_cfactor_muaverage(0,1,psd[0],psd[1]);

   transform_gflds();
   mul_cfactor_muaverage(0,1,psd[0],psd[2]);
   transform_sd(psd[2],psd[3]);

   mulr_spinor_add_dble(VOLUME,psd[3],psd[1],-1.0);
   d=norm_square_dble(VOLUME,1,psd[3])/norm_square_dble(VOLUME,1,psd[0]);

   if (my_rank==0)
   {   
      printf("mu_average\n");
      printf("Normalized difference = %.2e\n",sqrt(d));
      printf("(should be around 1*10^(-15) or so)\n\n");
   }


   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
