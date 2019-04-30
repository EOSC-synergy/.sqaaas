
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2005, 2008-2013, 2016 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of plaq_*frc(), sw_*frc(), hop_*frc().
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
#include "mdflds.h"
#include "forces.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,n,count;
   double phi[2],phi_prime[2],theta[3],cF[2];
   double wt1,wt2,wdt;
   FILE *flog=NULL;
   spinor_dble **wsd;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);

      printf("\n");
      printf("Timing of plaq_*frc(), sw_*frc(), hop_*frc()\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   
   set_flds_parms(3,0);
   
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);

   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
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
   set_dirac_parms9(-3,-0.0123,0.95,0.8,cF[0],cF[1],theta[0],theta[1],theta[2]);

   print_flds_parms();
   print_bc_parms();
   print_dirac_parms();

   start_ranlux(0,12345);
   geometry();

   alloc_wsd(2);
   wsd=reserve_wsd(2);

   random_ud();
   random_ad();
   random_sd(VOLUME,wsd[0],1.0);
   random_sd(VOLUME,wsd[1],1.0);
   bnd_sd2zero(ALL_PTS,wsd[0]);
   bnd_sd2zero(ALL_PTS,wsd[1]);

   plaq_u1frc();
   plaq_su3frc();
   set_su3frc2zero();
   set_u1frc2zero();
   set_xt2zero();
   add_prod2xt(-0.5,wsd[0],wsd[1]);
   add_prod2xv(-0.5,wsd[0],wsd[1]);
   sw_su3frc(1.0);
   hop_su3frc(1.0);
   sw_u1frc(1.0);
   hop_u1frc(1.0);

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         plaq_su3frc();
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
   {
      printf("Time per link:\n");
      printf("plaq_su3frc():      %4.3f usec\n",wdt);
   }

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         plaq_u1frc();
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
   {
      printf("plaq_u1frc():       %4.3f usec\n",wdt);
   }

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;
   set_xt2zero();

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         add_prod2xt(0.0,wsd[0],wsd[1]);
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
      printf("add_prod2xt():      %4.3f usec\n",wdt);

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;
   set_xv2zero();

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         add_prod2xv(0.0,wsd[0],wsd[1]);
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
      printf("add_prod2xv():      %4.3f usec\n",wdt);

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;
   set_su3frc2zero();

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         sw_su3frc(0.0);
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
      printf("sw_su3frc():        %4.3f usec\n",wdt);

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;
   set_u1frc2zero();

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         sw_u1frc(0.0);
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
      printf("sw_u1frc():         %4.3f usec\n",wdt);

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;
   set_su3frc2zero();

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         hop_su3frc(0.0);
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
      printf("hop_su3frc():       %4.3f usec\n",wdt);

   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;
   set_u1frc2zero();

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
         hop_u1frc(0.0);
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
   {
      printf("hop_u1frc():        %4.3f usec\n\n",wdt);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
