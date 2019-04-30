
/*******************************************************************************
*
* File ratfcts.c
*
* Copyright (C) 2012 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Rational function coefficients data base
*
* The externally accessible functions are
*
*   ratfct_t ratfct(int *irat)
*     Returns a structure containing the coefficients of the rational
*     function specified by the integers irat[3] (see the notes).
*
* Notes:
*
* Let R(x) be the [n,n]-rational approximation of x^w in the range
* ra^2<=x<=rb^2 which minimizes the relative error
*
*   delta = max_x | 1 - x^(-w) R(x) |
*
* The rational function is represented as
*
*   R(x) = A prod_{j=0}^{n-1} (x+Nu[j]^2)/(x+Mu[j]^2)
*
* The rational function is retrieved from the paramter database. In particular
* ratfct(irat) uses the rational approximation returned by rat_parms(irat[0]).
* The dictionary between the symbols in the previous formulae and the member
* variables of rat_parms(irat[0]) follows.
*
*   rp=rat_parms(irat[0]);
*   
*   w     :=: (double)(rp.power[0])/(double)(rp.power[1])
*   n     :=: rp.degree
*   A     :=: rp.A
*   delta :=: rp.delta
*   Nu[j] :=: rp.nu[j]   for j=0,...,n-1
*   Mu[j] :=: rp.mu[j]   for j=0,...,n-1
*
* Setting
*
*   k := irat[1]
*   l := irat[2]
*   
* the function ratfct(irat) returns a structure that identifies the rational
* function constructed by selecting a subset of terms from R(x), i.e.
*
*   r(x) = prod_{j=k}^{l} (x+Nu[j]^2)/(x+Mu[j]^2)
*        = prod_{j=0}^{l-k} (x+nu[j]^2)/(x+mu[j]^2)
*
* its partial fraction decomposition
*
*   r(y) = 1 + sum_{j=0}^{l-k} rmu[j]/(x^2+mu[j]^2)
*
* and the following auxiliary partial fraction decomposition 
*
*   prod_{j=0}^{l-k} (x+i*mu[j])/(x+i*nu[j])
*
*       = 1 + i* sum_{j=0}^{l-k} rnu[j]/(x+i*nu[j])
*
* The dictionary between the symbols in the previous formulae and the member
* variables of ratfct(irat) follows.
*
*   rf=ratfct(irat)
*
*   A             :=: rf.A
*   delta         :=: rf.delta
*   l-k+1         :=: rf.np
*   Nu[j+k]=nu[j] :=: rf.nu[j]    for j=0,...,np-1
*   Mu[j+k]=mu[j] :=: rf.mu[j]    for j=0,...,np-1
*   rnu[j]        :=: rf.rnu[j]   for j=0,...,np-1
*   rmu[j]        :=: rf.rmu[j]   for j=0,...,np-1
*
* Summarizing, the parameters of the ratfct function are
*
*   irat[0]       Index of the rational function in the
*                 parameter data base.
*
*   irat[1]       Lower end k of the selected coefficient range.       
*
*   irat[2]       Upper end l of the selected coefficient range.
*
* The members of the ratfct_t structure are
*
*   np        Number of poles of r(x) (that is, np=l-k+1).
*
*   power     Power to be approximated.
*
*   A         Normalization factor.
*
*   delta     Relative approximation error.
*
*   mu        Square-roots of the (-1)*zeroes of the denominator of r(x)
*
*   rmu       Array of the associated residues.
*
*   nu        Square-roots of the (-1)*zeroes of the numerator of r(x)
*
*   rnu       Array of the associated residues.
*
* The program in this module may perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define RATFCTS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "ratfcts.h"

#define IRMAX 32

static int init=0,irs,irats[IRMAX][3];
static ratfct_t rats[IRMAX]={{0,0.0,0.0,1.0,NULL,NULL,NULL,NULL}};


static void init_rat(void)
{
   int ir;
   
   for (ir=0;ir<IRMAX;ir++)
   {
      if (ir>0)
         rats[ir]=rats[0];

      irats[ir][0]=0;
      irats[ir][1]=0;
      irats[ir][2]=0;
   }

   irs=0;
   init=1;
}


static int fnd_rat(int *irat)
{
   int ir;

   for (ir=0;ir<irs;ir++)
   {
      if ((irat[0]==irats[ir][0])&&(irat[1]==irats[ir][1])&&
          (irat[2]==irats[ir][2]))
         return ir;
   }

   return irs;
}


static void alloc_rat(int *irat)
{
   int n,np,k,l;
   double *mu;
   rat_parms_t rp;
   
   error(irs==IRMAX,1,"alloc_rat [ratfcts.c]",
         "Attempt to define more than %d rational functions",IRMAX);

   rp=rat_parms(irat[0]);
   n=rp.degree;
   k=irat[1];
   l=irat[2];
   np=l-k+1;
   
   error((k<0)||(l<k)||(l>=n),1,"alloc_rat [ratfcts.c]",
         "Improper coefficient range or undefined rational function");

   mu=malloc(4*np*sizeof(*mu));   

   error(mu==NULL,1,"alloc_rat [ratfcts.c]",
         "Unable to allocate coefficient arrays");
   
   rats[irs].np=np;
   rats[irs].mu=mu;
   rats[irs].rmu=mu+np;
   rats[irs].nu=mu+2*np;
   rats[irs].rnu=mu+3*np;

   irats[irs][0]=irat[0];
   irats[irs][1]=irat[1];
   irats[irs][2]=irat[2];               
}


static void set_rat(int *irat)
{
   int np,k,l,i,j;
   double pmu,pnu;
   double *mu,*nu,*rmu,*rnu;
   rat_parms_t rp;

   rp=rat_parms(irat[0]);
   k=irat[1];
   l=irat[2];
   np=l-k+1;

   rats[irs].A=rp.A;
   rats[irs].delta=rp.delta;
   rats[irs].power=(1.0*rp.power[0])/(rp.power[1]);

   mu=rats[irs].mu;
   nu=rats[irs].nu;
   rmu=rats[irs].rmu;
   rnu=rats[irs].rnu;

   for (i=0;i<np;i++)
   {
      mu[i]=rp.mu[k+i];
      nu[i]=rp.nu[k+i];
   }

   for (i=0;i<np;i++)
   {  
      pmu=1.0;
      pnu=1.0;

      for (j=0;j<np;j++)
      {
         if (j!=i)
         {
            pmu*=((nu[j]*nu[j]-mu[i]*mu[i])/(mu[j]*mu[j]-mu[i]*mu[i]));
            pnu*=((mu[j]-nu[i])/(nu[j]-nu[i]));
         }
      }

      rmu[i]=(nu[i]*nu[i]-mu[i]*mu[i])*pmu;
      rnu[i]=(mu[i]-nu[i])*pnu;
   }

   irs+=1;
}


ratfct_t ratfct(int *irat)
{
   int ir;
   
   if (init==0)
      init_rat();

   ir=fnd_rat(irat);

   if (ir==irs)
   {
      alloc_rat(irat);
      set_rat(irat);
   }

   return rats[ir];
}
