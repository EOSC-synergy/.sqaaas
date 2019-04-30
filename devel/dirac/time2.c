
/*******************************************************************************
*
* File time2.c
*
* Copyright (C) 2005, 2008, 2011-2013, 2016 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of Dw_dble() and Dwhat_dble().
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
#include "hflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,gg,count,nt;
   int i,nflds,qhat;
   double phi[2],phi_prime[2],theta[3],cF[2],su3csw,u1csw;
   double mu,wt1,wt2,wdt;
   spinor_dble **psd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time2.log","w",stdout);

      printf("\n");
      printf("Timing of Dw_dble() and Dwhat_dble()\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      if (NPROC>1)
         printf("There are %d MPI processes\n",NPROC);
      else
         printf("There is 1 MPI process\n");

      if ((VOLUME*sizeof(double))<(64*1024))
      {
         printf("The local size of the gauge field is %d KB\n",
                (int)((72*VOLUME*sizeof(double))/(1024)));
         printf("The local size of a quark field is %d KB\n",
                (int)((24*VOLUME*sizeof(double))/(1024)));
      }
      else
      {
         printf("The local size of the gauge field is %d MB\n",
                (int)((72*VOLUME*sizeof(double))/(1024*1024)));
         printf("The local size of a quark field is %d MB\n",
                (int)((24*VOLUME*sizeof(double))/(1024*1024)));
      }

#if (defined x64)
#if (defined AVX)
      printf("Using AVX instructions\n");
#else
      printf("Using SSE3 instructions and 16 xmm registers\n");
#endif
#if (defined P3)
      printf("Assuming SSE prefetch instructions fetch 32 bytes\n");
#elif (defined PM)
      printf("Assuming SSE prefetch instructions fetch 64 bytes\n");
#elif (defined P4)
      printf("Assuming SSE prefetch instructions fetch 128 bytes\n");
#else
      printf("SSE prefetch instructions are not used\n");
#endif
#endif
      printf("\n");

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time2.c]",
                    "Syntax: time2 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [time2.c]",
                    "Syntax: time2 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      gg=find_opt(argc,argv,"-gg");

      if (gg!=0)
         error_root(sscanf(argv[gg+1],"%d",&gg)!=1,1,"main [time2.c]",
                    "Syntax: time2 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");
      else
         gg=1;
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);
   
   set_flds_parms(gg,0);
   
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,0,cs,phi,phi_prime);

   qhat=0;
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   if ((gauge()&1)!=0) su3csw=0.95;
   if ((gauge()&2)!=0)
   {
      qhat=-3;
      u1csw=0.8;
   }
   if (bc_type()!=3)
   {
      cF[0]=0.9012;
      cF[1]=1.2034;
   }
   if (bc_cstar()==0)
   {
      theta[0]=0.35;
      theta[1]=-1.25;
      theta[2]=0.78;
   }
   set_dirac_parms9(qhat,-0.0123,su3csw,u1csw,cF[0],cF[1],
                    theta[0],theta[1],theta[2]);

   print_flds_parms();
   print_bc_parms();
   print_dirac_parms();

   start_ranlux(0,12345);
   geometry();

   mu=0.0785;

   if ((gauge()&1)!=0) random_ud();
   if ((gauge()&2)!=0) random_ad();
   sw_term(NO_PTS);
   hdfld();

   nflds=(int)((4*1024*1024)/(VOLUME*sizeof(double)))+1;
   if ((nflds%2)==1)
      nflds+=1;
   alloc_wsd(nflds);
   psd=reserve_wsd(nflds);

   for (i=0;i<nflds;i++)
      random_sd(VOLUME,psd[i],1.0);

   nt=(int)(1.0e6/(double)(nflds*VOLUME));
   if (nt<2)
      nt=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         for (i=0;i<nflds;i+=2)
            Dw_dble(mu,psd[i],psd[i+1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=4.0e6*wdt/((double)(nt)*(double)(nflds*VOLUME));

   if (my_rank==0)
   {
      printf("Time per lattice point for Dw_dble():\n");
      printf("%4.3f micro sec (%d Mflops)\n\n",wdt,(int)(1920.0/wdt));
   }

   nt=(int)(1.0e6/(double)(nflds*VOLUME));
   if (nt<2)
      nt=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<nt;count++)
      {
         for (i=0;i<nflds;i+=2)
            Dwhat_dble(mu,psd[i],psd[i+1]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      nt*=2;
   }

   wdt=4.0e6*wdt/((double)(nt)*(double)(nflds*VOLUME));

   if (my_rank==0)
   {
      printf("Time per lattice point for Dwhat_dble():\n");
      printf("%4.3f micro sec (%d Mflops)\n\n",wdt,(int)(1908.0/wdt));
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
