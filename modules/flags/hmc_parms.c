
/*******************************************************************************
*
* File hmc_parms.c
*
* Copyright (C) 2009, 2010, 2011, 2013 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic HMC parameters
*
* The externally accessible functions are
*
*   hmc_parms_t set_hmc_parms(int nact,int *iact,int npf,int nmu,
*                             double *mu,int nlv,double tau)
*     Sets some basic parameters of the HMC algorithm. The parameters are
*
*       nact        Number of terms in the total action
*                     
*       iact        Indices iact[i] of the action terms (i=0,..,nact-1) 
*
*       npf         Number of pseudo-fermion fields on which the action
*                   depends
*
*       nmu         Number of twisted mass parameters on which the
*                   pseudo-fermion actions and forces depend
*
*       mu          Twisted masses mu[i] (i=0,..,nmu-1)
*
*       nlv         Number of levels of the molecular-dynamics integrator
*
*       tau         Molecular-dynamics trajectory length
*
*     The total action must include the gauge action, but pseudo-fermion
*     actions are optional and the momentum action is treated separately.
*     The program returns a structure that contains the parameters listed
*     above.
*
*   hmc_parms_t hmc_parms(void)
*     Returns a structure containing the current values of the parameters
*     listed above.
*
*   void print_hmc_parms(void)
*     Prints the HMC parameters to stdout on MPI process 0.
*
*   void write_hmc_parms(FILE *fdat)
*     Writes the HMC parameters to the file fdat on MPI process 0.
*
*   void check_hmc_parms(FILE *fdat)
*     Compares the HMC parameters with the values stored on the file fdat
*     on MPI process 0, assuming the latter were written to the file by
*     the program write_hmc_parms().
*
* Notes:
*
* To ensure the consistency of the data base, the parameters must be set
* simultaneously on all processes. The type hmc_parms_t is defined in the
* in the file flags.h.
*
*******************************************************************************/

#define HMC_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

static hmc_parms_t hmc={0,0,0,0,NULL,0.0,NULL};


hmc_parms_t set_hmc_parms(int nact,int *iact,int npf,int nmu,
                          double *mu,int nlv,double tau)
{
   int i;

   check_global_int("set_hmc_parms",4,nact,npf,nmu,nlv);
   check_global_dble("set_hmc_parms",1,tau);
   check_global_intarray("set_hmc_parms",nact,iact);
   check_global_dblearray("set_hmc_parms",nmu,mu);
   
   error_root((npf<0)||(nlv<1),1,"set_hmc_parms [hmc_parms.c]",
              "Improper number of pseudo-fermion fields or integrator levels");

   error_root((hmc.nlv>0)&&(npf!=hmc.npf),1,"set_hmc_parms [hmc_parms.c]",
              "Number of pseudo-fermion fields may be set only once");
   
   if (nact!=hmc.nact)
   {
      if (hmc.iact!=NULL)
      {
         free(hmc.iact);
         hmc.iact=NULL;
      }

      if (nact>0)
      {
         hmc.iact=malloc(nact*sizeof(int));
         error(hmc.iact==NULL,1,"set_hmc_parms [hmc_parms.c]",
               "Unable to allocate parameter array");
      }
   }

   if (nmu!=hmc.nmu)
   {
      if (hmc.mu!=NULL)
      {
         free(hmc.mu);
         hmc.mu=NULL;
      }

      if (nmu>0)
      {
         hmc.mu=malloc(nmu*sizeof(double));
         error(hmc.mu==NULL,2,"set_hmc_parms [hmc_parms.c]",
               "Unable to allocate parameter array");
      }
   }

   hmc.nact=nact;
   hmc.npf=npf;
   hmc.nmu=nmu;
   hmc.nlv=nlv;
   hmc.tau=tau;

   for (i=0;i<nact;i++)
      hmc.iact[i]=iact[i];

   for (i=0;i<nmu;i++)
      hmc.mu[i]=mu[i];
   
   return hmc;
}


hmc_parms_t hmc_parms(void)
{
   return hmc;
}


void print_hmc_parms(void)
{
   int my_rank,n,i;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      printf("HMC parameters:\n");
      printf("actions =");
      for (i=0;i<hmc.nact;i++)
         printf(" %d",hmc.iact[i]);
      printf("\n");
      printf("npf = %d\n",hmc.npf);      
      if (hmc.nmu>0)
      {
         printf("mu =");
         for (i=0;i<hmc.nmu;i++)
         {
            n=fdigits(hmc.mu[i]);
            printf(" %.*f",IMAX(n,1),hmc.mu[i]);
         }
         printf("\n");
      }
      printf("nlv = %d\n",hmc.nlv);
      n=fdigits(hmc.tau);
      printf("tau = %.*f\n\n",IMAX(n,1),hmc.tau);
   }
}


void write_hmc_parms(FILE *fdat)
{
   write_little_int(fdat,4,hmc.nact,hmc.npf,hmc.nmu,hmc.nlv);
   write_little_intarray(fdat,hmc.nact,hmc.iact);
   write_little_dble(fdat,1,hmc.tau);
   write_little_dblearray(fdat,hmc.nmu,hmc.mu);
}


void check_hmc_parms(FILE *fdat)
{
   check_fpar_int("check_hmc_parms",fdat,4,hmc.nact,hmc.npf,hmc.nmu,hmc.nlv);
   check_fpar_intarray("check_hmc_parms",fdat,hmc.nact,hmc.iact);
   check_fpar_dble("check_hmc_parms",fdat,1,hmc.tau);
   check_fpar_dblearray("check_hmc_parms",fdat,hmc.nmu,hmc.mu);
}
