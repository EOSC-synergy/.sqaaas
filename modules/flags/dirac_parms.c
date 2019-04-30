
/*******************************************************************************
*
* File dirac_parms.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/modules/flags/lat_parms.c
* Copyright (C) 2009-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Dirac-operator parameters.
*
* The externally accessible functions are
*
*   dirac_parms_t set_dirac_parms9(int qhat,double m0,
*                                  double su3csw,double u1csw,
*                                  double cF,double cF_prime,
*                                  double th1,double th2,double th3)
*     Sets the parameters of the Dirac operator. The adjustable parameters are
*
*       m0             Base Wilson mass.
*
*       qhat           Integer electric charge (in units of the elementary one).
*
*       su3csw         Coefficient of the SU(3) Sheikholeslami-Wohlert term.
*
*       u1csw          Coefficient of the U(1) Sheikholeslami-Wohlert term.
*
*       cF,cF_prime    Fermion action improvement coefficients at time 0 and T,
*                      respectively.
*
*       th1,th2,th3    Angles specifying the phase-periodic boundary
*                      conditions for the quark fields in direction 1,2,3.       
*
*     The return value is a structure that contains all above parameters.
*
*   dirac_parms_t dirac_parms_t set_dirac_parms1(dirac_parms_t *par)
*     Sets the parameters of the Dirac operator, copying them from the structure
*     (*par).
*
*   dirac_parms_t dirac_parms(void)
*     Returns the parameters currently set for the Dirac operator.
*
*   void print_dirac_parms(void);
*     Prints the parameters of the Dirac operator to stdout on MPI process 0.
*
*   tm_parms_t set_tm_parms(int eoflg)
*     Sets the twisted-mass flag. The parameter is
*
*       eoflg          If the flag is set (eoflg!=0), the twisted-mass term
*                      in the Dirac operator, the SAP preconditioner and the
*                      little Dirac operator is turned off on the odd lattice
*                      sites.
*
*     The return value is a structure that contains the twisted-mass flag.
*
*   tm_parms_t tm_parms(void)
*     Returns a structure containing the twisted-mass flag.
*
* Notes:
*
* To ensure the consistency of the data base, the parameters must be set
* simultaneously on all processes. The data type dirac_parms_t is defined in the
* file flags.h.
*
*******************************************************************************/

#define DIRAC_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static dirac_parms_t dp={0,DBL_MAX,0.0,0.0,{0.0,0.0},{0.0,0.0,0.0}};
static tm_parms_t tm={0};


dirac_parms_t set_dirac_parms9(int qhat,double m0,
                               double su3csw,double u1csw,
                               double cF,double cF_prime,
                               double th1,double th2,double th3)
{
   check_global_int("set_dirac_parms",1,qhat);
   check_global_dble("set_dirac_parms",8,m0,su3csw,u1csw,cF,cF_prime,
                                         th1,th2,th3);

   error_root(((gauge()&1)==0)&&(su3csw!=0.0),1,
              "set_dirac_parms9 [dirac_parms.c]",
              "su3csw must be 0 if SU(3) is not dynamical");

   error_root(((gauge()&2)==0)&&(u1csw!=0.0),1,
              "set_dirac_parms9 [dirac_parms.c]",
              "u1csw must be 0 if U(1) is not dynamical");

   error_root(((gauge()&2)==0)&&(qhat!=0),1,
              "set_dirac_parms9 [dirac_parms.c]",
              "qhat must be 0 if U(1) is not dynamical");

   error_root((bc_type()==3)&&((cF!=0.0)||(cF_prime!=0.0)),1,
              "set_dirac_parms9 [dirac_parms.c]",
              "cF and cF' must be 0 with periodic boundary conditions");

   error_root((bc_cstar()!=0)&&((th1!=0.0)||(th2!=0.0)||(th3!=0.0)),1,
              "set_dirac_parms9 [dirac_parms.c]",
              "theta[k] must be 0 with C-periodic boundary conditions");

   if ((su3csw!=dp.su3csw)||(u1csw!=dp.u1csw)||
       (cF!=dp.cF[0])||(cF_prime!=dp.cF[1]))
   {
      set_flags(ERASED_SW);
      set_flags(ERASED_SWD);
      set_grid_flags(SAP_BLOCKS,ERASED_SW);
      set_flags(ERASED_AW);
      set_flags(ERASED_AWHAT);
   }

   if (qhat!=dp.qhat)
   {
      set_flags(ERASED_H);
      set_flags(ERASED_HD);
      set_grid_flags(SAP_BLOCKS,ERASED_HBGR);
      set_flags(ERASED_SW);
      set_flags(ERASED_SWD);
      set_grid_flags(SAP_BLOCKS,ERASED_SW);
      set_flags(ERASED_AW);
      set_flags(ERASED_AWHAT);
   }
   
   if (m0!=dp.m0)
   {
      set_flags(ERASED_SW);
      set_flags(ERASED_SWD);
      set_grid_flags(SAP_BLOCKS,ERASED_SW);
      set_flags(ERASED_AWHAT);
   }
   
   if ((th1!=dp.theta[0])||(th2!=dp.theta[1])||(th3!=dp.theta[2]))
   {
      set_flags(ERASED_H);
      set_flags(ERASED_HD);
      set_grid_flags(SAP_BLOCKS,ERASED_HBGR);
      set_flags(ERASED_AW);
      set_flags(ERASED_AWHAT);
   }

   dp.qhat=qhat;
   dp.m0=m0;
   dp.su3csw=su3csw;
   dp.u1csw=u1csw;
   dp.cF[0]=cF;
   dp.cF[1]=cF_prime;
   dp.theta[0]=th1;
   dp.theta[1]=th2;
   dp.theta[2]=th3;

   return dp;
}


dirac_parms_t set_dirac_parms1(dirac_parms_t *par)
{
   return set_dirac_parms9((*par).qhat,(*par).m0,(*par).su3csw,
                           (*par).u1csw,(*par).cF[0],(*par).cF[1],
                           (*par).theta[0],(*par).theta[1],(*par).theta[2]);
}


dirac_parms_t dirac_parms(void)
{
   return dp;
}


void print_dirac_parms(void)
{
   int my_rank,n;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      printf("Dirac-operator parameters:\n");

      printf("qhat = %d\n",dp.qhat);
      if (dp.m0==DBL_MAX)
      {
         printf("m0 = infty\n");
      }
      else
      {
         n=fdigits(dp.m0);
         printf("m0 = %.*f\n",IMAX(n,1),dp.m0);
      }
      n=fdigits(dp.su3csw);
      printf("su3csw = %.*f\n",IMAX(n,1),dp.su3csw);
      n=fdigits(dp.u1csw);
      printf("u1csw = %.*f\n",IMAX(n,1),dp.u1csw);
      n=fdigits(dp.cF[0]);
      printf("cF = %.*f\n",IMAX(n,1),dp.cF[0]);
      n=fdigits(dp.cF[1]);
      printf("cF' = %.*f\n",IMAX(n,1),dp.cF[1]);
      n=fdigits(dp.theta[0]);
      printf("theta = %.*f,",IMAX(n,1),dp.theta[0]);
      n=fdigits(dp.theta[1]);
      printf("%.*f,",IMAX(n,1),dp.theta[1]);
      n=fdigits(dp.theta[2]);
      printf("%.*f\n",IMAX(n,1),dp.theta[2]);
      printf("\n");
      fflush(stdout);
   }
}


tm_parms_t set_tm_parms(int eoflg)
{
   check_global_int("set_tm_parms",1,eoflg);

   if (eoflg!=tm.eoflg)
      set_flags(ERASED_AWHAT);

   tm.eoflg=eoflg;

   return tm;
}


tm_parms_t tm_parms(void)
{
   return tm;
}
