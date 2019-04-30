
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/devel/forces/time1.c
* Copyright (C) 2005, 2008-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of hdfld() and ud1fld(LOC).
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
#include "hflds.h"
#include "global.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,gg,n,count,qhat,qhatmax;
   double phi[2],phi_prime[2],theta[3];
   double wt1,wt2,wdt;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);

      printf("\n");
      printf("Timing of hdfld() and ud1fld(LOC)\n");
      printf("---------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      gg=find_opt(argc,argv,"-gg");

      if (gg!=0)
         error_root(sscanf(argv[gg+1],"%d",&gg)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");
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

   print_flds_parms();
   print_bc_parms();

   theta[0]=theta[1]=theta[2]=0.0;
   if (bc_cstar()==0)
   {
      theta[0]=0.35;
      theta[1]=-1.25;
      theta[2]=0.78;
   }

   start_ranlux(0,12345);
   geometry();

   if ((gauge()&1)!=0) random_ud();
   if ((gauge()&2)!=0) random_ad();

   hdfld();
   
   if (gauge()==1)
   {
      set_dirac_parms9(0,-0.0123,0.0,0.0,0.0,0.0,theta[0],theta[1],theta[2]);

      n=(int)(3.0e6/(double)(4*VOLUME));
      if (n<2)
         n=2;
      wdt=0.0;

      while (wdt<5.0)
      {
         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();
         for (count=0;count<n;count++)
         {
            set_flags(UPDATED_AD);
            hdfld();
         }
         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();

         wdt=wt2-wt1;
         n*=2;
      }

      wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

      if (my_rank==0)
      {
         printf("Time per link\n");
         printf("hdfld():      %4.6f usec\n",wdt);
      }
   }
   else
   {
      n=(int)(3.0e6/(double)(4*VOLUME));
      if (n<2)
         n=2;
      wdt=0.0;

      while (wdt<5.0)
      {
         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();
         for (count=0;count<n;count++)
         {
            set_flags(UPDATED_AD);
            u1dfld(LOC);
         }
         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();

         wdt=wt2-wt1;
         n*=2;
      }

      wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

      if (my_rank==0)
      {
         printf("\nu1dfld(LOC):      %4.6f usec\n\n\n",wdt);
         fflush(flog);
      }

      if (my_rank==0)
      {
         printf("The three quoted values in each of the following lines correspond\n");
         printf("to the time needed for:\n");
         printf("1) calculating hdfld from new fundamental fields;\n");
         printf("2) calculating hdfld assuming that u1dfld is up-to-date;\n");
         printf("3) retrieving the hdfld pointer assuming that hdfld is up-to-date.\n");
         printf("Each line corresponds to a different value of the electric charge.\n\n");
         printf("Time per link, hdfld():\n");
      }

      qhatmax=0;
      if((gauge()&2)!=0) qhatmax=8;
      for(qhat=0;qhat<=qhatmax;qhat++)
      {
         set_dirac_parms9(qhat,-0.0123,0.0,0.0,0.0,0.0,theta[0],theta[1],theta[2]);

         n=(int)(3.0e6/(double)(4*VOLUME));
         if (n<2)
            n=2;
         wdt=0.0;

         while (wdt<5.0)
         {
            MPI_Barrier(MPI_COMM_WORLD);
            wt1=MPI_Wtime();
            for (count=0;count<n;count++)
            {
               set_flags(UPDATED_AD);
               hdfld();
            }
            MPI_Barrier(MPI_COMM_WORLD);
            wt2=MPI_Wtime();

            wdt=wt2-wt1;
            n*=2;
         }

         wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

         if (my_rank==0)
         {
            printf("q = %+.2d:      %4.6f",qhat,wdt);
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
            {
               set_flags(ERASED_HD);
               hdfld();
            }
            MPI_Barrier(MPI_COMM_WORLD);
            wt2=MPI_Wtime();

            wdt=wt2-wt1;
            n*=2;
         }

         wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

         if (my_rank==0)
         {
            printf("      %4.6f",wdt);
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
            {
               hdfld();
            }
            MPI_Barrier(MPI_COMM_WORLD);
            wt2=MPI_Wtime();

            wdt=wt2-wt1;
            n*=2;
         }

         wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

         if (my_rank==0)
         {
            printf("      %4.6f usec\n",wdt);
            fflush(flog);
         }
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
