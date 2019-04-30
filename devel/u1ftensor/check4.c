
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2016 Nazario Tantalo
*               2017 Agostino Patella
*
* Based on openQCD-1.6/devel/tcharge/check4.c
* Copyright (C) 2010, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the gauge and translation invariance of the U(1) fluxes.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "u1ftensor.h"
#include "global.h"
#include "gflds_utils.h"


static void random_vec(int *svec)
{
   int mu,bs[4];
   double r[4];

   bs[0]=(NPROC0*L0);
   bs[1]=(NPROC1*L1);
   bs[2]=(NPROC2*L2);
   bs[3]=(NPROC3*L3);

   ranlxd(r,4);

   for (mu=0;mu<4;mu++)
   {
      svec[mu]=(int)((double)(bs[mu])*r[mu]);
      if (svec[mu]>(bs[mu]/2))
         svec[mu]-=bs[mu];
   }

   MPI_Bcast(svec,4,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int my_rank,i,n,bc,s[4];
   double phi[2],phi_prime[2];
   double d,dmax1[6],dmax2[6];
   double Phi1[6],Phi2[6],phi1,phi2;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      printf("\n");
      printf("Gauge and translation invariance of the U(1) fluxes\n");
      printf("---------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>]");
   }

   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   set_u1lat_parms(0,1.0/137.0,5.0,0.0,7.0,0.0,0.0,0);
   print_lat_parms();

   start_ranlux(0,12345);
   geometry();

   for(n=0;n<6;++n)
   {
      dmax1[n]=0.0;
      dmax2[n]=0.0;
   }

   for (i=0;i<8;i++)
   {
      random_gflds();
         
      for(n=0;n<6;++n)
	 Phi1[n]=u1fluxes(n);
      random_vec(s);
      if (bc!=3)      
         s[0]=0;
      shift_gflds(s);
      for(n=0;n<6;++n)
	 Phi2[n]=u1fluxes(n);

      for(n=0;n<6;++n)
      {
	 d=fabs(Phi1[n]-Phi2[n])/fabs(Phi1[n]);
	 if ((d>dmax1[n]) && (fabs(Phi1[n])> 1.0e-12))
	    dmax1[n]=d;
      }

      random_g();
      transform_gflds();
      for(n=0;n<6;++n)
	 Phi2[n]=u1fluxes(n);
      
      for(n=0;n<6;++n)
      {
	 d=fabs(Phi1[n]-Phi2[n])/fabs(Phi1[n]);
	 if ((d>dmax2[n]) && (fabs(Phi1[n])> 1.0e-12))
	    dmax2[n]=d;
      
	 phi1=Phi1[n];
	 phi2=Phi2[n];

	 MPI_Bcast(&phi1,1,MPI_INT,0,MPI_COMM_WORLD);
	 MPI_Bcast(&phi2,1,MPI_INT,0,MPI_COMM_WORLD);

	 error((phi1!=Phi1[n])||(phi2!=Phi2[n]),1,"main [check4.c]",
	       "Fluxes are not globally the same");
      }
   }
   
   print_flags();
   
   if (my_rank==0)
   {
      printf("Translation invariance = ");
      for(n=0;n<6;++n)
	 printf("%.2e  ",dmax1[n]);
      printf("\n");
      printf("Gauge invariance       = ");
      for(n=0;n<6;++n)
	 printf("%.2e  ",dmax2[n]);
      printf("\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
