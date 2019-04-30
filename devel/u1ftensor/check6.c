
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2016 Nazario Tantalo
*               2017 Agostino Patella
*
* Based on openQCD-1.6/devel/tcharge/check6.c
* Copyright (C) 2010, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program u1fluxes_slices().
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
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "u1ftensor.h"
#include "global.h"
#include "gflds_utils.h"

#define N0 (NPROC0*L0)

static int bc;
static double Phi1[6],Phi2[6],Phi3[6],Phi[6*N0],Phi0[6*N0];


int main(int argc,char *argv[])
{
   int my_rank,i,imax,t,n,cs;
   double phi[2],phi_prime[2];   
   double dev,dev2;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);
      
      printf("\n");
      printf("Check of the program u1fluxes_slices()\n");
      printf("--------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);      

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-cs <cstar>]");
   }
   
   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,cs,phi,phi_prime);
   print_bc_parms();

   set_u1lat_parms(0,1.0/137.0,5.0,0.0,7.0,0.0,0.0);
   print_lat_parms();

   start_ranlux(0,123456);   
   geometry();

   imax=10;
   
   for (i=0;i<imax;i++)
   {
      random_gflds();

      for (n=0;n<6;n++)
      {      
         Phi1[n]=u1fluxes(n);
         Phi2[n]=u1fluxes_slices(n,Phi+n*N0);
      }

      for (n=0;n<6;n++)
      {
	 dev=Phi1[n];
	 dev2=Phi1[n];
	 for (t=0;t<N0;t++)
	 {
	    dev-=Phi[t+n*N0];
	    Phi0[t+n*N0]=Phi[t+n*N0];
	 }
	 Phi3[n]=Phi2[n];
	 dev2-=Phi2[n];

	 dev=fabs(dev)/fabs(Phi1[n]);
	 dev2=fabs(dev2)/fabs(Phi1[n]);
      
	 if (my_rank==0)
	 {
	    printf("i=%3d, Phi[%d]=%12.4e, dev=%7.1e, dev2=%7.1e, Phi[%d][0...%d]=%10.2e",
		   i+1,n,Phi1[n],dev,dev2,n,N0-1,Phi0[0+n*N0]);

	    for (t=1;t<N0;t++)
	       printf(", %10.2e",Phi0[t+n*N0]);

	    printf("\n");
	 }

	 MPI_Bcast(Phi0+n*N0,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);

	 for (t=0;t<N0;t++)
	 {      
	    if ((Phi[t+n*N0]-Phi0[t+n*N0])!=0.0)
	       break;
	 }

	 error(t!=N0,1,"main [check6.c]",
	       "Fluxes slices are not globally the same");

	 MPI_Bcast(Phi3+n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	 error((Phi3[n]-Phi2[n])!=0.0,1,"main [check6.c]",
	       "Fluxes are not globally the same");

      }

      if (my_rank==0)
	 printf("\n\n");
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();    
   exit(0);
}

#undef N0
