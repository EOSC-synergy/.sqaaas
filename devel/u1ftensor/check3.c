
/*******************************************************************************
*
* File check3.c
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
* Check of the program mxw_action_slices().
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
static double A1,A2,A[N0],A0[N0];


int main(int argc,char *argv[])
{
   int my_rank,i,imax,t;
   double phi[2],phi_prime[2];   
   double dev;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      
      printf("\n");
      printf("Check of the program mxw_action_slices()\n");
      printf("----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);      

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }
   
   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();

   set_u1lat_parms(0,1.0/137.0,5.0,0.0,7.0,0.0,0.0);
   print_lat_parms();

   start_ranlux(0,123456);   
   geometry();

   imax=10;
   
   for (i=0;i<imax;i++)
   {
      random_gflds();

      A1=mxw_action();
      A2=mxw_action_slices(A);
      dev=fabs(A1-A2)/A1;
      
      for (t=0;t<N0;t++)
      {
         A2-=A[t];
         A0[t]=A[t];
      }

      dev+=fabs(A2)/A1;
      
      if (my_rank==0)
      {
         printf("i=%3d, A=%12.6e, dev=%6.1e, A[0...%d]=%8.2e",
                i+1,A1,dev,N0-1,A[0]);

         for (t=1;t<N0;t++)
            printf(", %8.2e",A[t]);

         printf("\n");
      }

      MPI_Bcast(A0,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);

      for (t=0;t<N0;t++)
      {      
         if ((A[t]-A0[t])!=0.0)
            break;
      }

      error(t!=N0,1,"main [check3.c]",
            "Action slices are not globally the same");
   }

   if (my_rank==0)
   {    
      printf("\n");
      fclose(flog);
   }   

   MPI_Finalize();    
   exit(0);
}

#undef N0
