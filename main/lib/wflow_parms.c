
/*******************************************************************************
*
* File wflow_parms.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/main/qcd1.c
*      and openQCD-1.6/main/ym1.c
* Copyright (C) 2011-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "global.h"
#include "mpi.h"
#include "main_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static int init=0;
static wflow_parms_t wfp={0,0,0,0,0.0};


wflow_parms_t wflow_parms(void)
{
   error_root(!init,1,"wflow_parms [wflow_parms.c]",
              "Wilson flow parameters are not initialized");
   return wfp;
}


void read_wflow_parms(void)
{
   int my_rank,nstep,dnms,flint;
   double eps;
   char line[NAME_SIZE];

   error_root(init,1,"read_wflow_parms [wflow_parms.c]",
              "Wilson flow parameters can be initialized only once");
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      find_section("Wilson flow");
      read_line("integrator","%s",line);
      read_line("eps","%lf",&eps);
      read_line("nstep","%d",&nstep);
      read_line("dnms","%d",&dnms);

      if (strcmp(line,"EULER")==0)
         flint=0;
      else if (strcmp(line,"RK2")==0)
         flint=1;
      else if (strcmp(line,"RK3")==0)
         flint=2;
      else
         error_root(1,1,"read_wflow_parms [wflow_parms.c]","Unknown integrator");

      error_root((dnms<1)||(nstep<dnms)||((nstep%dnms)!=0),1,
                 "read_wflow_parms [wflow_parms.c]",
                 "nstep must be a multiple of dnms");
   }

   MPI_Bcast(&flint,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dnms,1,MPI_INT,0,MPI_COMM_WORLD);

   wfp.dn=dnms;
   wfp.nn=nstep/dnms;
   wfp.tmax=N0;
   wfp.eps=eps;
   wfp.flint=flint;
   
   init=1;
}


void check_wflow_parms(FILE *fpar)
{
   error_root(!init,1,"check_wflow_parms [wflow_parms.c]",
              "Wilson flow parameters are not initialized");
   check_little_int("read_wflow_parms [wflow_parms.c]",fpar,3,
                     wfp.flint,wfp.nn*wfp.dn,wfp.dn);
   check_little_dble("read_wflow_parms [wflow_parms.c]",fpar,1,wfp.eps);
}


void write_wflow_parms(FILE *fpar)
{
   error_root(!init,1,"write_wflow_parms [wflow_parms.c]",
              "Wilson flow parameters are not initialized");
   write_little_int(1,fpar,3,wfp.flint,wfp.nn*wfp.dn,wfp.dn);
   write_little_dble(1,fpar,1,wfp.eps);
}


void print_wflow_parms(void)
{
   int n;

   error_root(!init,1,"print_wflow_parms [wflow_parms.c]",
              "Wilson flow parameters are not initialized");
   
   printf("Wilson flow:\n");
   if (wfp.flint==0)
      printf("Euler integrator\n");
   else if (wfp.flint==1)
      printf("2nd order RK integrator\n");
   else
      printf("3rd order RK integrator\n");
   n=fdigits(wfp.eps);
   printf("eps = %.*f\n",IMAX(n,1),wfp.eps);
   printf("nstep = %d\n",wfp.dn*wfp.nn);
   printf("dnms = %d\n\n",wfp.dn);
}
