
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2017 Agostino Patella
* 
* Based on openQCD-1.6/devel/uflds/check6.c
* Copyright (C) 2005, 2009, 2010, 2011 Martin Luescher, Filippo Palombi
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks of the programs in the module dalg
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "linalg.h"
#include "global.h"

#define NMOM 3200033



int main(int argc,char *argv[])
{
   int my_rank,n;
   double dev,dmax,dmax_all;
   double nsq1,nsq2,sprod1,sprod2;
   double sm;
   double rn,cij,eij;
   double var,var_all;
   double *X,*Y,*x;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);

      printf("\n");
      printf("Checks of the programs in the module dalg\n");
      printf("-----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      printf("Number of momenta: %d\n\n",NMOM);
   }

   start_ranlux(0,123456);
   geometry();

   X=amalloc(2*NMOM*sizeof(*X),4);
   error((X==NULL),1,
         "main [check6.c]","Unable to allocate field arrays");
   Y=X+NMOM;

   set_dvec2zero(NMOM,X);
   dmax=0.0;

   for (n=0;n<NMOM;n++)
   {
      dev=fabs(X[n]);
      if (dev>dmax)
         dmax=dev;

      X[n]=1.0;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Check of set_dvec2zero():\n\n");
      printf("max|X| = %.1e (should be 0.0)\n\n",dmax_all);
   }

   dmax=fabs(norm_square_dvec(NMOM,1,X)-(double)(NMOM*NPROC));
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Check of norm_square_dvec():\n\n");
      printf("Element count = %.1e (should be 0.0)\n\n",dmax_all);
   }

   sm=0.0;
   dmax=0.0;

   for (n=0;n<NMOM;n++)
   {
      x=X+n;

      random_dvec(1,x);

      nsq1=norm_square_dvec(1,0,x);

      nsq2=(*x)*(*x);

      dev=fabs(1.0-nsq2/nsq1);
      if (dev>dmax)
         dmax=dev;

      sm+=nsq2;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("|1.0-(X,X)/||X||^2| = %.1e (single elements)\n",dmax_all);
      printf("(should be less than %.1e or so)\n\n",DBL_EPSILON*sqrt(8.0));
   }

   dmax=fabs(1.0-sm/norm_square_dvec(NMOM,0,X));
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("|1.0+2*(X,X)/||X||^2| = %.1e (whole vector)\n",dmax_all);
      printf("(should be less than %.1e or so)\n\n",
             DBL_EPSILON*sqrt((double)(NMOM)));
   }

   random_dvec(NMOM,X);
   random_dvec(NMOM,Y);

   nsq1=norm_square_dvec(NMOM,1,X);
   nsq2=norm_square_dvec(NMOM,1,Y);
   sprod1=scalar_prod_dvec(NMOM,1,X,Y);

   for (n=0;n<NMOM;n++)
   {
      X[n]+=Y[n];
   }

   sprod2=0.5*(norm_square_dvec(NMOM,1,X)-nsq1-nsq2);

   if (my_rank==0)
   {
      printf("Check of scalar_prod_alg(): %.1e\n",
             fabs(sprod1-sprod2)/sqrt(nsq1*nsq2));
      printf("(should be less than %.1e or so)\n\n",
             DBL_EPSILON*sqrt((double)(NMOM)*(double)(NPROC)));
   }

   random_dvec(NMOM,X);

   var=0.0;

   for (n=0;n<NMOM;n++)
   {
      var+=X[n]*X[n];
   }

   MPI_Reduce(&var,&var_all,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Check of random_dvec():\n\n");
      dmax=0.0;
      rn=1.0/((double)(NMOM)*(double)(NPROC));

      cij=1.0;
      eij=sqrt(2.0*rn);

      var_all*=rn;
      printf("<b*b> = % .4e, deviation = %.1e+-%.1e\n\n",
             var_all,fabs(var_all-cij),eij);
   }

   rn=-1.2345;
   random_dvec(NMOM,X);
   random_dvec(NMOM,Y);

   nsq1=norm_square_dvec(NMOM,1,X);
   nsq2=norm_square_dvec(NMOM,1,Y);
   sprod1=scalar_prod_dvec(NMOM,1,X,Y);

   muladd_assign_dvec(NMOM,rn,X,Y);
   sm=norm_square_dvec(NMOM,1,Y)-nsq2-rn*rn*nsq1-2.0*rn*sprod1;
   sm=fabs(sm)/nsq1;

   if (my_rank==0)
   {
      printf("Check of muladd_assign_dvec(): %.1e\n",sm);
      printf("(should be less than 1.0e-15 or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
