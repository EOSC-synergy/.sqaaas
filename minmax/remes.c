
/*******************************************************************************
*
* File remes.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the rational approximation to f(x)= ( x + mu^2 )^{p/q}
*
* The externally accessible functions are
*
*   void alloc_remes(remes_t *rem)
*     the function must be called to allocate the basics arrays of the
*     remes_t structure before starting any calculation. The variable
*     rem->precision si used as input and must be set before calling the function.
*
*   void remes(remes_t *rem)
*     the function runs Remes's second algorithm until a MinMax approximation
*     of the function f(x)= ( x + mu^2 ^{p,q} with an error delta<= goal has been 
*     found. The rem variable passed to the function must be allocated by calling 
*     the function alloc_remes(). The following members of the *rem structure are
*     used as inputes and must be set:
*
*     rem->p,          the function to be approximated is f(x)= ( x + mu^2 )^{p/q}
*     rem->q,
*     rem->mu,
*     rem->ra,         the function is approximated in the range [ra^2,rb^2]
*     rem->rb,
*     rem->goal,       the final approximation is such that delta<= goal 
*     rem->relflag,    relflag=1, delta= relative error, 
*                      relflag=0, delta= absolute error
*     rem->verbose,    verbose=1, the openQcdQed output is written,
*                      verbose=0, the openQcdQed output is not written
*     rem->path        the results are written in the directory path/
*     rem->nstart      the starting point of the iteration
*                      if nstart>1, then also the variables
*                      rem->alpha and rem->beta are used   
*     
* Notes: 
*
* The routines:
*
*     ludcmp
*     lubksb
*     vandermonde
*     balanc
*     hqr
*     zrhqr
*
* are my multiple-precision implementation of the corresponding routines 
* given in
*
*    Press, W.H. et al. Numerical Recipes in C, 
*                       Cambridge University Press, Second Edition, 2002
*
* See the README file for further explanations.
* 
*******************************************************************************/

#define REMES_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpf2mpfr.h"
#include "utils.h"
#include "remes.h"

#define E2S0FACT     (1.15)
#define E1S0FACT     (0.85)
#define E2S1FACT     (1.1)
#define EE1FACT      (0.5)
#define EE2FACT      (0.6)
#define SOLVEIMAX    (200)
#define BRACKETIMAX  (100)
#define BRACKETBND   (0.55)
#define SECONDIMAX   (10000)
#define HQRIMAX      (200)
#define REMESSDX     (5.0e-2)
#define REMESEPS     (1.0e-50)
#define REMESREPS    (1.0e-18)


static void ludcmp(mpfr_t *a,int n,int *indx)
{
   int i,imax,j,k;
   mpfr_t big,dum,sum;
   mpfr_t temp,d;
   mpfr_t *vv;   

   imax=-1;

   vv=malloc(n*sizeof(mpfr_t));
   error(vv==NULL,1,"ludcmp [remes.c]","Unable to allocate auxiliary array");

   mpfr_inits(big,dum,sum,temp,d,NULL);
   for (i=0;i<n;i++) 
      mpfr_init(vv[i]);

   mpfr_set_d(d,1.0,ROUNDING);
   for (i=0;i<n;i++) 
   {
      mpfr_set_d(big,0.0,ROUNDING);
      for(j=0;j<n;j++)
      {
	 mpfr_abs(temp,a[i+n*j],ROUNDING);
	 if (mpfr_greater_p(temp,big)) 
	    mpfr_set(big,temp,ROUNDING);
      }

      error(mpfr_cmp_d(big,0.0)<=0,1,"ludcmp [remes.c]","Singular matrix in routine ludcmp");

      mpfr_d_div(vv[i],1.0,big,ROUNDING);
   }

   for (j=0;j<n;j++) 
   {
      for (i=0;i<j;i++) 
      {
	 mpfr_set(sum,a[i+n*j],ROUNDING);
	 for (k=0;k<i;k++) 
	 {
	    mpfr_mul(dum,a[i+n*k],a[k+n*j],ROUNDING);
	    mpfr_sub(sum,sum,dum,ROUNDING);
	 }
	 mpfr_set(a[i+n*j],sum,ROUNDING);
      }

      mpfr_set_d(big,0.0,ROUNDING);
      for (i=j;i<n;i++) 
      {
	 mpfr_set(sum,a[i+n*j],ROUNDING);
	 for (k=0;k<j;k++)
	 {
	    mpfr_mul(dum,a[i+n*k],a[k+n*j],ROUNDING);
	    mpfr_sub(sum,sum,dum,ROUNDING);
	 }
	 mpfr_set(a[i+n*j],sum,ROUNDING);

	 mpfr_abs(dum,sum,ROUNDING);
	 mpfr_mul(dum,dum,vv[i],ROUNDING);

	 if ( mpfr_greaterequal_p(dum,big) ) 
	 {
	    mpfr_set(big,dum,ROUNDING);
	    imax=i;
	 }
      }

      if (j != imax) 
      {
	 for (k=0;k<n;k++) 
	 {
	    mpfr_set(dum,a[imax+n*k],ROUNDING);
	    mpfr_set(a[imax+n*k],a[j+n*k],ROUNDING);
	    mpfr_set(a[j+n*k],dum,ROUNDING);
	 }
	 mpfr_neg(d,d,ROUNDING);
	 mpfr_set(vv[imax],vv[j],ROUNDING);
      }

      indx[j]=imax;

      if (j!=n-1) 
      {
	 mpfr_d_div(dum,1.0,a[j+n*j],ROUNDING);
	 for(i=j+1;i<n;i++) 
	    mpfr_mul(a[i+n*j],a[i+n*j],dum,ROUNDING);
      }
   }

   mpfr_clears(big,dum,sum,temp,d,NULL);
   for (i=0;i<n;i++) 
      mpfr_clear(vv[i]);

   free(vv);
}


static void lubksb(mpfr_t *a,int n,int *indx,mpfr_t *b)
{
   int i,ii,ip,j;
   mpfr_t sum,temp;

   mpfr_inits(sum,temp,NULL);

   ii=-1;

   for (i=0;i<n;i++) 
   {
      ip=indx[i];
      mpfr_set(sum,b[ip],ROUNDING);
      mpfr_set(b[ip],b[i],ROUNDING);
      if (ii>=0)
	 for (j=ii;j<i;j++)
	 {
	    mpfr_mul(temp,a[i+n*j],b[j],ROUNDING);
	    mpfr_sub(sum,sum,temp,ROUNDING);
	 } 
      else if (mpfr_cmp_d(sum,0.0)) 
	 ii=i;
      mpfr_set(b[i],sum,ROUNDING);
   }

   for (i=n-1;i>=0;i--) 
   {
      mpfr_set(sum,b[i],ROUNDING);
      for (j=i+1;j<n;j++)
      {
	 mpfr_mul(temp,a[i+n*j],b[j],ROUNDING);
	 mpfr_sub(sum,sum,temp,ROUNDING); 
      }
      mpfr_div(b[i],sum,a[i+n*i],ROUNDING);
   }

   mpfr_clears(sum,temp,NULL);
}


void vandermonde(mpfr_t x[],mpfr_t y[],int n,mpfr_t cof[])
{
   int k,j,i;
   mpfr_t t1,t2,phi,ff,b,*s;

   s=malloc(n*sizeof(mpfr_t));
   error(s==NULL,1,"vandermonde [remes.c]","Unable to allocate auxiliary array");

   mpfr_inits(t1,t2,phi,ff,b,NULL);

   for (i=0;i<n;i++) 
   {
      mpfr_init(s[i]);
      mpfr_set_d(s[i],0.0,ROUNDING);
      mpfr_set_d(cof[i],0.0,ROUNDING);
   }
   mpfr_neg(s[n-1],x[0],ROUNDING);

   for (i=1;i<n;i++) 
   {
      for (j=n-i-1;j<n-1;j++)
      {
	 mpfr_mul(t1,x[i],s[j+1],ROUNDING);
	 mpfr_sub(s[j],s[j],t1,ROUNDING);
      }
      mpfr_sub(s[n-1],s[n-1],x[i],ROUNDING);
   }

   for (j=0;j<n;j++) 
   {
      mpfr_set_d(phi,(double)(n),ROUNDING);
      for (k=n-1;k>=1;k--)
      {
	 mpfr_set_d(t1,(double)(k),ROUNDING);
	 mpfr_mul(t1,t1,s[k],ROUNDING);
	 mpfr_mul(t2,x[j],phi,ROUNDING);

	 mpfr_add(phi,t1,t2,ROUNDING);
      }

      mpfr_div(ff,y[j],phi,ROUNDING);

      mpfr_set_d(b,1.0,ROUNDING);
      for (k=n-1;k>=0;k--) 
      {
	 mpfr_mul(t1,b,ff,ROUNDING);
	 mpfr_add(cof[k],cof[k],t1,ROUNDING);

	 mpfr_mul(t2,x[j],b,ROUNDING);
	 mpfr_add(b,s[k],t2,ROUNDING);
      }
   }

   mpfr_clears(t1,t2,phi,ff,b,NULL);
   for (i=0;i<n;i++) 
      mpfr_clear(s[i]);

   free(s);
}


static void function(remes_t *rem,mpfr_t f,mpfr_t x)
{
   mpfr_t y;

   mpfr_init(y);

   mpfr_mul(y,rem->bma,x,ROUNDING);
   mpfr_add(y,y,rem->bpa,ROUNDING);

   mpfr_pow(f,y,rem->power,ROUNDING);

   mpfr_clear(y);
}


static void weight(remes_t *rem,mpfr_t w,mpfr_t x)
{
   if(rem->relflag==0)
      mpfr_set_d(w,1.0,ROUNDING);
   else
   {
      function(rem,w,x);
      mpfr_d_div(w,1.0,w,ROUNDING);
   }
}


static  void numerator(remes_t *rem,mpfr_t N,mpfr_t x)
{
   int i;
   mpfr_t Nx;

   mpfr_init(Nx);

   mpfr_set(N,rem->cn[rem->n],ROUNDING);
   for(i=rem->n-1;i>=0;--i)
   {
      mpfr_mul(Nx,N,x,ROUNDING);
      mpfr_add(N,Nx,rem->cn[i],ROUNDING);
   }

   mpfr_clear(Nx);
}


static void denominator(remes_t *rem,mpfr_t D,mpfr_t x)
{
   int i;
   mpfr_t Dx;

   mpfr_init(Dx);

   mpfr_set_d(D,0.0,ROUNDING);
   for(i=rem->d;i>0;--i)
   {
      mpfr_add(Dx,D,rem->cd[i],ROUNDING);
      mpfr_mul(D,Dx,x,ROUNDING);
   }
   mpfr_add_d(D,D,1.0,ROUNDING);

   mpfr_clear(Dx);
}


static void approximation(remes_t *rem,mpfr_t a,mpfr_t x)
{
   mpfr_t N,D;

   mpfr_inits(N,D,NULL);

   numerator(rem,N,x);
   denominator(rem,D,x);

   mpfr_div(a,N,D,ROUNDING);

   mpfr_clears(N,D,NULL);
}


static void deviation(remes_t *rem,mpfr_t e,mpfr_t x)
{
   mpfr_t w,f,a;

   mpfr_inits(w,f,a,NULL);

   function(rem,f,x);
   if(rem->relflag==0)
      mpfr_set_d(w,1.0,ROUNDING);
   else
      mpfr_d_div(w,1.0,f,ROUNDING);

   approximation(rem,a,x);

   mpfr_sub(e,f,a,ROUNDING);
   mpfr_mul(e,e,w,ROUNDING);
   
   mpfr_clears(w,f,a,NULL);
}


static void rfunction(remes_t *rem,mpfr_t f,mpfr_t x)
{
   mpfr_t y;

   mpfr_inits(y,NULL);

   mpfr_add(y,x,rem->mu2,ROUNDING);
   mpfr_pow(f,y,rem->power,ROUNDING);

   mpfr_clears(y,NULL);
}



static void Cset_n(mpfr_t *c,mpfr_t re)
{
   mpfr_set(c[0],re,ROUNDING);      
   mpfr_set_d(c[1],0.0,ROUNDING);      
}

static void Cadd(mpfr_t *c,mpfr_t *a,mpfr_t *b)
{
   mpfr_add(c[0],a[0],b[0],ROUNDING);      
   mpfr_add(c[1],a[1],b[1],ROUNDING);      
}

static void Cmul(mpfr_t *c,mpfr_t *a,mpfr_t *b)
{
   mpfr_t t1,t2,r1,r2;

   mpfr_inits(t1,t2,r1,r2,NULL);

   mpfr_mul(t1,a[0],b[0],ROUNDING);      
   mpfr_mul(t2,a[1],b[1],ROUNDING);      

   mpfr_sub(r1,t1,t2,ROUNDING);      

   mpfr_mul(t1,a[0],b[1],ROUNDING);      
   mpfr_mul(t2,a[1],b[0],ROUNDING);      

   mpfr_add(r2,t1,t2,ROUNDING);      

   mpfr_set(c[0],r1,ROUNDING);      
   mpfr_set(c[1],r2,ROUNDING);      

   mpfr_clears(t1,t2,r1,r2,NULL);
}

static  void rnumerator(remes_t *rem,mpfr_t N,mpfr_t x)
{
   int i;
   mpfr_t R[2],Nx[2],cx[2];

   mpfr_inits(R[0],R[1],
	      Nx[0],Nx[1],
	      cx[0],cx[1],
	      NULL
      );

   Cset_n(cx,x);   
   Cset_n(R,rem->cn[rem->n]);

   for(i=rem->n-1;i>=0;--i)
   {
      Cadd(Nx,cx,rem->rn+2*i);
      Cmul(R,R,Nx);
   }

   mpfr_set(N,R[0],ROUNDING);      

   mpfr_clears(R[0],R[1],
	      Nx[0],Nx[1],
	      cx[0],cx[1],
	      NULL
      );
}


static  void rdenominator(remes_t *rem,mpfr_t D,mpfr_t x)
{
   int i;
   mpfr_t R[2],Dx[2],cx[2];

   mpfr_inits(R[0],R[1],
	      Dx[0],Dx[1],
	      cx[0],cx[1],
	      NULL
      );

   Cset_n(cx,x);
   Cset_n(R,rem->cd[rem->d]);

   for(i=rem->d-1;i>=0;--i)
   {
      Cadd(Dx,cx,rem->rd+2*i);
      Cmul(R,R,Dx);
   }

   mpfr_set(D,R[0],ROUNDING);      

   mpfr_clears(R[0],R[1],
	      Dx[0],Dx[1],
	      cx[0],cx[1],
	      NULL
      );
}


static void rapproximation(remes_t *rem,mpfr_t a,mpfr_t x)
{
   mpfr_t N,D;

   mpfr_inits(N,D,NULL);

   rnumerator(rem,N,x);
   rdenominator(rem,D,x);

   mpfr_div(a,N,D,ROUNDING);

   mpfr_clears(N,D,NULL);
}


static void rdeviation(remes_t *rem,mpfr_t e,mpfr_t x)
{
   mpfr_t w,f,a;

   mpfr_inits(w,f,a,NULL);

   rfunction(rem,f,x);
   if(rem->relflag==0)
      mpfr_set_d(w,1.0,ROUNDING);
   else
      mpfr_d_div(w,1.0,f,ROUNDING);

   rapproximation(rem,a,x);

   mpfr_sub(e,f,a,ROUNDING);
   mpfr_mul(e,e,w,ROUNDING);
   
   mpfr_clears(w,f,a,NULL);
}


static void abs_deviation(remes_t *rem,mpfr_t e,mpfr_t x)
{
   deviation(rem,e,x);
   mpfr_abs(e,e,ROUNDING);      
}

static void search_emax(remes_t *rem,int nx,mpfr_t emax)
{
   int i;
   mpfr_t ei;

   mpfr_init(ei);

   mpfr_abs(emax,rem->e[0],ROUNDING);
   for(i=1;i<nx;++i)
   {
      mpfr_abs(ei,rem->e[i],ROUNDING);
      if (mpfr_greater_p(ei,emax))
	 mpfr_set(emax,ei,ROUNDING);
   }

   mpfr_clear(ei);
}


static void search_emin(remes_t *rem,int nx,mpfr_t emin)
{
   int i;
   mpfr_t ei;

   mpfr_init(ei);

   mpfr_abs(emin,rem->e[0],ROUNDING);
   for(i=1;i<nx;++i)
   {
      mpfr_abs(ei,rem->e[i],ROUNDING);
      if (mpfr_less_p(emin,ei))
	 mpfr_set(emin,ei,ROUNDING);
   }

   mpfr_clear(ei);
}


static void alloc_remes_nd(remes_t *rem,int n,int d)
{
   int i,nd2;

   nd2=n+d+2;

   rem->cn=malloc(11*nd2*sizeof(mpfr_t));
   error(rem->cn==NULL,1,"alloc_remes_nd [remes.c]","Unable to allocate remes structure");
 
   rem->cd =rem->cn+1*nd2;
   rem->x  =rem->cn+2*nd2;
   rem->y  =rem->cn+3*nd2;
   rem->dx =rem->cn+4*nd2;
   rem->e  =rem->cn+5*nd2;
   rem->se =rem->cn+6*nd2;
   rem->rn =rem->cn+7*nd2;
   rem->rd =rem->cn+9*nd2;

   for(i=0;i<11*nd2;++i)
      mpfr_init(rem->cn[i]);

   mpfr_init(rem->sdx);
   mpfr_init(rem->alpha);
   mpfr_init(rem->beta);
   mpfr_init(rem->eeps);
   mpfr_init(rem->eeps2);
   mpfr_init(rem->reeps);
   mpfr_init(rem->reeps2);
   mpfr_init(rem->ee);
   mpfr_init(rem->ra);
   mpfr_init(rem->rb);
   mpfr_init(rem->mu);
   mpfr_init(rem->mu2);
   mpfr_init(rem->a);
   mpfr_init(rem->b);
   mpfr_init(rem->bpa);
   mpfr_init(rem->bma);
   mpfr_init(rem->emax);
   mpfr_init(rem->emin);
   mpfr_init(rem->eprev);
   mpfr_init(rem->power);
   mpfr_init(rem->goal);
}


void alloc_remes(remes_t *rem)
{
   double precision;

   error(rem->precision<64,1,"alloc_remes [remes.c]",
	 "Parameter precision has not been set properly, precision>=64");

   precision=3.321928094*((double)(rem->precision));
   mpfr_set_default_prec(precision);

   alloc_remes_nd(rem,0,0);
}


static void free_remes(remes_t *rem)
{
   int i;

   for(i=0;i<11*rem->nd2;++i)
      mpfr_clear(rem->cn[i]);

   free(rem->cn);

   mpfr_clear(rem->sdx);
   mpfr_clear(rem->alpha);
   mpfr_clear(rem->beta);
   mpfr_clear(rem->eeps);
   mpfr_clear(rem->eeps2);
   mpfr_clear(rem->reeps);
   mpfr_clear(rem->reeps2);
   mpfr_clear(rem->ee);
   mpfr_clear(rem->ra);
   mpfr_clear(rem->rb);
   mpfr_clear(rem->mu);
   mpfr_clear(rem->mu2);
   mpfr_clear(rem->a);
   mpfr_clear(rem->b);
   mpfr_clear(rem->bpa);
   mpfr_clear(rem->bma);
   mpfr_clear(rem->emax);
   mpfr_clear(rem->emin);
   mpfr_clear(rem->eprev);
   mpfr_clear(rem->power);
   mpfr_clear(rem->goal);

}


static void copy_parameters(remes_t *next,remes_t *prev)
{
   int i;

   sprintf(next->path,"%s",prev->path);

   next->verbose=prev->verbose;

   next->relflag=prev->relflag;
   next->precision=prev->precision;
   next->p=prev->p;
   next->q=prev->q;

   next->nfit=prev->nfit;
   next->nstart=prev->nstart;

   mpfr_set(next->power,prev->power,ROUNDING);
   mpfr_set(next->goal,prev->goal,ROUNDING);

   mpfr_set(next->ra,prev->ra,ROUNDING);
   mpfr_set(next->rb,prev->rb,ROUNDING);

   mpfr_set(next->mu,prev->mu,ROUNDING);
   mpfr_set(next->mu2,prev->mu2,ROUNDING);

   mpfr_set(next->a,prev->a,ROUNDING);
   mpfr_set(next->b,prev->b,ROUNDING);

   mpfr_set(next->bpa,prev->bpa,ROUNDING);
   mpfr_set(next->bma,prev->bma,ROUNDING);

   mpfr_set(next->sdx,prev->sdx,ROUNDING);

   mpfr_set(next->alpha,prev->alpha,ROUNDING);
   mpfr_set(next->beta,prev->beta,ROUNDING);

   mpfr_set(next->eeps,prev->eeps,ROUNDING);
   mpfr_set(next->eeps2,prev->eeps2,ROUNDING);
   mpfr_set(next->reeps,prev->reeps,ROUNDING);
   mpfr_set(next->reeps2,prev->reeps2,ROUNDING);

   mpfr_set(next->emax,prev->emax,ROUNDING);   
   mpfr_set(next->emin,prev->emin,ROUNDING);   
   mpfr_set(next->eprev,prev->eprev,ROUNDING);   

   next->nx=prev->nx;

   for(i=0;i<next->nx;++i)
   {
      mpfr_set(next->x[i],prev->x[i],ROUNDING);
      mpfr_set(next->dx[i],prev->dx[i],ROUNDING);
      mpfr_set(next->e[i],prev->e[i],ROUNDING);
      mpfr_set(next->se[i],prev->se[i],ROUNDING);
   }
}


static void assign_nextden(remes_t *next,remes_t *prev)
{
   int i,n;
   
   copy_parameters(next,prev);

   next->n=prev->n;
   next->d=prev->d+1;

   next->nd=next->n+next->d;
   next->nd2=next->nd+2;

   n=next->d;
   for(i=0;i<n;++i)
      mpfr_set(next->cd[i],prev->cd[i],ROUNDING);
   mpfr_set_d(next->cd[n],0.0,ROUNDING);

   n=next->n+1;
   for(i=0;i<n;++i)
      mpfr_set(next->cn[i],prev->cn[i],ROUNDING);
}


static void assign_nextnum(remes_t *next,remes_t *prev)
{
   int i,n;

   copy_parameters(next,prev);

   next->n=prev->n+1;
   next->d=prev->d;

   next->nd=next->n+next->d;
   next->nd2=next->nd+2;

   n=next->d+1;
   for(i=0;i<n;++i)
      mpfr_set(next->cd[i],prev->cd[i],ROUNDING);

   n=next->n;
   for(i=0;i<n;++i)
      mpfr_set(next->cn[i],prev->cn[i],ROUNDING);
   mpfr_set_d(next->cn[n],0.0,ROUNDING);
}


static void assign_nextrem(remes_t *next,remes_t *prev)
{
   int i,n;

   if(next->cn!=NULL)
      free_remes(next);
   alloc_remes_nd(next,prev->n,prev->d);

   copy_parameters(next,prev);

   next->n=prev->n;
   next->d=prev->d;
   next->nd=next->n+next->d;
   next->nd2=next->nd+2;

   n=next->d+1;
   for(i=0;i<n;++i)
      mpfr_set(next->cd[i],prev->cd[i],ROUNDING);

   n=next->n+1;
   for(i=0;i<n;++i)
      mpfr_set(next->cn[i],prev->cn[i],ROUNDING);
}


static void set_nfit(remes_t *rem)
{
   mpfr_t a,len,lenm1;

   mpfr_inits(a,len,lenm1,NULL);

   mpfr_log(len,rem->emax,ROUNDING);
   mpfr_log(lenm1,rem->eprev,ROUNDING);

   mpfr_sub(rem->alpha,lenm1,len,ROUNDING);

   mpfr_mul_d(rem->beta,rem->alpha,(double)(rem->n),ROUNDING);
   mpfr_add(rem->beta,rem->beta,len,ROUNDING);
   
   mpfr_log(a,rem->goal,ROUNDING);
   mpfr_sub(a,rem->beta,a,ROUNDING);
   mpfr_div(a,a,rem->alpha,ROUNDING);

   rem->nfit=(int)(mpfr_get_d(a,ROUNDING))+1;

   printf("The error scales as exp(%.5f-n*%.5f) and the iteration will stop at\n\n\t n ~ %d.\n\n",
	  mpfr_get_d(rem->beta,ROUNDING),mpfr_get_d(rem->alpha,ROUNDING),rem->nfit);
   printf("To save time you may want to restart the program with the option:\n\n\t -nstart %d %.5f %.5f\n",
	  rem->nfit-1,mpfr_get_d(rem->alpha,ROUNDING),mpfr_get_d(rem->beta,ROUNDING)
      );
   printf("\n\n");

   mpfr_clears(a,len,lenm1,NULL);
}


static void starting_algorithm_cheby(remes_t *rem)
{
   int i;
   mpfr_t pio2n,x,n,k;
   remes_t next;   

   mpfr_set(rem->eprev,rem->emax,ROUNDING);

   alloc_remes_nd(&next,rem->n+1,rem->d+1);
   assign_nextden(&next,rem);
   assign_nextnum(&next,&next);
   assign_nextrem(rem,&next);

   mpfr_inits(pio2n,x,n,k,NULL);

   mpfr_set_d(n,(double)(rem->nd),ROUNDING);
   mpfr_set_d(pio2n,1.0,ROUNDING);
   mpfr_asin(pio2n,pio2n,ROUNDING);
   mpfr_div(pio2n,pio2n,n,ROUNDING);

   mpfr_set_d(rem->x[0],-1.0,ROUNDING);
   mpfr_set_d(rem->x[rem->nd2-1],1.0,ROUNDING);
   for(i=1;i<=rem->nd;++i)
   {
      mpfr_set_d(k,(double)(2*i-1),ROUNDING);
      mpfr_mul(x,k,pio2n,ROUNDING);
      mpfr_cos(x,x,ROUNDING);

      mpfr_set(rem->x[rem->nd-i+1],x,ROUNDING);
   }
   rem->nx=rem->nd2;

   search_emin(rem,rem->nx,rem->emin);
   search_emax(rem,rem->nx,rem->emax);

   if(rem->n>2 || rem->nstart>0)
   {
      mpfr_set_d(rem->emax,(double)(-rem->n),ROUNDING);
      mpfr_mul(rem->emax,rem->emax,rem->alpha,ROUNDING);
      mpfr_add(rem->emax,rem->emax,rem->beta,ROUNDING);
      mpfr_exp(rem->emax,rem->emax,ROUNDING);

   }

   mpfr_clears(pio2n,x,n,k,NULL);
}


static void second_algorithm_solve_fixed_E(remes_t *next,mpfr_t E)
{
   int i,j,n,*indx;
   mpfr_t sign,x,xx,f,w;
   mpfr_t *A,*b;

   n=next->nx-1;
   error(n!=next->n+next->d+1,1,"second_algorithm_solve [remes.c]","The iteration is not in proper stage");

   A=malloc((n*n+n)*sizeof(mpfr_t));
   error(A==NULL,1,"second_algorithm_solve [remes.c]","Unable to allocate auxiliary array");
   b=A+n*n;

   mpfr_inits(sign,x,xx,f,w,NULL);
   for (i=0;i<n*n+n;++i)
      mpfr_init(A[i]);

   indx=malloc(n*sizeof(int));
   error(indx==NULL,1,"second_algorithm_solve [remes.c]","Unable to allocate auxiliary array");

   mpfr_set_d(sign,-1.0,ROUNDING);
   for (i=0;i<n;++i)
   {
      mpfr_set(x,next->x[i],ROUNDING);

      function(next,f,x);
      weight(next,w,x);

      mpfr_mul(b[i],sign,E,ROUNDING);
      mpfr_div(b[i],b[i],w,ROUNDING);
      mpfr_sub(b[i],f,b[i],ROUNDING);

      mpfr_set_d(xx,1.0,ROUNDING);
      for (j=0;j<next->n+1;++j)
      {
	 mpfr_set(A[i+n*j],xx,ROUNDING);
	 mpfr_mul(xx,xx,x,ROUNDING);
      }
      
      mpfr_mul(xx,b[i],x,ROUNDING);
      mpfr_neg(xx,xx,ROUNDING);
      for (j=next->n+1;j<n;++j)
      {
	 mpfr_set(A[i+n*j],xx,ROUNDING);
	 mpfr_mul(xx,xx,x,ROUNDING);
      }
      
      mpfr_neg(sign,sign,ROUNDING);
   }      
   ludcmp(A,n,indx);
   lubksb(A,n,indx,b);

   for (j=0;j<next->n+1;++j)
      mpfr_set(next->cn[j],b[j],ROUNDING);

   for (j=next->n+1;j<n;++j)
      mpfr_set(next->cd[j-next->n],b[j],ROUNDING);

   for (i=0;i<next->nx;++i)
   {
      deviation(next,next->e[i],next->x[i]);
      mpfr_set_d(next->se[i],(double)mpfr_sgn(next->e[i]),ROUNDING);
   }

   mpfr_clears(sign,x,xx,f,w,NULL);
   for (i=0;i<n*n+n;++i)
      mpfr_clear(A[i]);

   free(A);
   free(indx);
}


static void second_algorithm_solve(remes_t *next)
{
   int i,l,imax,lmax,done,cmax,cmin;
   mpfr_t x,EE1,EE2,E1,E2,E3;
   mpfr_t sign,f1,f2,check,checks;
   mpfr_t max,min,re,rf,r;

   mpfr_inits(x,EE1,EE2,E1,E2,E3,
	      sign,f1,f2,check,checks,
	      max,min,re,rf,r,NULL);

   mpfr_set(x,next->x[next->nx-1],ROUNDING);

   mpfr_set_d(sign,1.0,ROUNDING);
   if(next->nx%2==1)
      mpfr_set_d(sign,-1.0,ROUNDING);

   if(next->sflag==0)
   {
      mpfr_mul(E2,next->se[0],next->emax,ROUNDING);
      mpfr_mul_d(E1,E2,E1S0FACT,ROUNDING);
      mpfr_mul_d(E2,E2,E2S0FACT,ROUNDING);

      mpfr_abs(max,next->emax,ROUNDING);
      mpfr_abs(min,next->emin,ROUNDING);      
   }
   else
   {
      mpfr_set(E1,next->ee,ROUNDING);
      mpfr_mul_d(E2,next->ee,E2S1FACT,ROUNDING);

      search_emax(next,next->nx,max);
      mpfr_set(min,next->ee,ROUNDING);
      mpfr_abs(min,min,ROUNDING);
   }
   mpfr_set(EE1,E1,ROUNDING);
   mpfr_set(EE2,E2,ROUNDING);

   done=0;
   imax=SOLVEIMAX;
   lmax=6;

   for(l=0;l<lmax;++l)
   {
      for(i=0;i<imax;++i)
      {
	 if(i==0)
	 {
	    second_algorithm_solve_fixed_E(next,E1);
	    deviation(next,f1,x);
	    mpfr_mul(r,sign,E1,ROUNDING);
	    mpfr_sub(f1,f1,r,ROUNDING);
	 }
	 else
	    mpfr_set(f1,f2,ROUNDING);

	 second_algorithm_solve_fixed_E(next,E2);
	 deviation(next,f2,x);
	 mpfr_mul(r,sign,E2,ROUNDING);
	 mpfr_sub(f2,f2,r,ROUNDING);
	 
	 mpfr_div(re,E1,E2,ROUNDING);
	 mpfr_neg(re,re,ROUNDING);
	 mpfr_add_d(re,re,1.0,ROUNDING);

	 mpfr_div(rf,f1,f2,ROUNDING);
	 mpfr_neg(rf,rf,ROUNDING);
	 mpfr_add_d(rf,rf,1.0,ROUNDING);

	 mpfr_div(r,re,rf,ROUNDING);
	 mpfr_neg(r,r,ROUNDING);
	 mpfr_add_d(r,r,1.0,ROUNDING);

	 mpfr_mul(E3,E2,r,ROUNDING);

	 mpfr_sub(check,E3,E2,ROUNDING);
	 mpfr_abs(check,check,ROUNDING);

	 mpfr_add(checks,E3,E2,ROUNDING);
	 mpfr_mul_d(checks,checks,0.5,ROUNDING);
	 mpfr_mul(checks,checks,next->reeps,ROUNDING);
	 mpfr_abs(checks,checks,ROUNDING);

	 mpfr_abs(r,E3,ROUNDING);
	 
	 cmax=1;
	 if (mpfr_greater_p(r,max) && next->sflag!=0)
	    cmax=0;
	 cmin=1;

	 if( (mpfr_less_p(check,next->eeps) || mpfr_less_p(check,checks)) && cmin && cmax )
	 {
	    done=1;
	    l=lmax;
	    i=imax;
	 }
	 else if( mpfr_less_p(check,next->eeps) && ( cmin==0 || cmax==0 ))
	    i=imax;
	 else
	 {
	    mpfr_set(E1,E2,ROUNDING);
	    mpfr_set(E2,E3,ROUNDING);
	 }
      }
      mpfr_mul_d(EE1,EE1,EE1FACT,ROUNDING);
      mpfr_mul_d(EE2,EE2,EE2FACT,ROUNDING);
      mpfr_set(E1,EE1,ROUNDING);
      mpfr_set(E2,EE2,ROUNDING);
   }
   error(done==0,1,"second_algorithm [remes.c]",
	 "Secant iteration did not converge in max num of iterations: try to rise the precision or p--> -p if the -abs option has not been set");

   mpfr_set(next->ee,E3,ROUNDING);
   second_algorithm_solve_fixed_E(next,E3);

   mpfr_clears(x,EE1,EE2,E1,E2,E3,
	       sign,f1,f2,check,checks,
	       max,min,re,rf,r,NULL);
}


static void bracket(remes_t *next,int i)
{
   int it,sec;
   mpfr_t dx;
   mpfr_t x0,x1,x2,b;
   mpfr_t x,xp,xm,xn;
   mpfr_t cx1,cx2,cx3;
   mpfr_t fx,fp,fm,fn;

   mpfr_inits(dx,
	      x0,x1,x2,b,
	      x,xp,xm,xn,
	      cx1,cx2,cx3,
	      fx,fp,fm,fn,
	      NULL);

   sec=1;

   mpfr_set(x,next->x[i],ROUNDING);
   mpfr_set(dx,next->dx[i],ROUNDING);

   mpfr_sub(xm,x,dx,ROUNDING);
   mpfr_add(xp,x,dx,ROUNDING);

   if ( mpfr_cmp_d(xm,-1.0)<=0 )
      mpfr_set_d(xm,-1.0,ROUNDING);

   if ( mpfr_cmp_d(xp,1.0)>=0 )
      mpfr_set_d(xp,1.0,ROUNDING);

   abs_deviation(next,fm,xm);
   abs_deviation(next,fx,x);
   abs_deviation(next,fp,xp);

   if ( mpfr_greater_p(fm,fp) )
   {
      mpfr_neg(dx,dx,ROUNDING);

      mpfr_set(xp,xm,ROUNDING);
      mpfr_set(fp,fm,ROUNDING);
            
      if(i==0)
	 mpfr_set_d(b,-1.0,ROUNDING);
      else
      {
	 mpfr_sub(cx1,next->x[i],next->x[i-1],ROUNDING);
	 mpfr_mul_d(cx1,cx1,BRACKETBND,ROUNDING);
	 mpfr_add(b,next->x[i-1],cx1,ROUNDING);
      }
   }
   else
   {
      if(i==next->nx-1)
	 mpfr_set_d(b,1.0,ROUNDING);
      else
      {
	 mpfr_sub(cx1,next->x[i+1],next->x[i],ROUNDING);
	 mpfr_mul_d(cx1,cx1,BRACKETBND,ROUNDING);
	 mpfr_sub(b,next->x[i+1],cx1,ROUNDING);
      }
   }

   mpfr_add(xn,xp,dx,ROUNDING);

   if ( mpfr_cmp_d(dx,0.0)>0 && mpfr_greaterequal_p(xn,b) )
      mpfr_set(xn,b,ROUNDING);

   if ( mpfr_cmp_d(dx,0.0)<0 && mpfr_lessequal_p(xn,b) )
      mpfr_set(xn,b,ROUNDING);

   abs_deviation(next,fn,xn);

   it=1;
   while ( mpfr_less_p(fp,fn) )
   {
      mpfr_set(xp,xn,ROUNDING);
      mpfr_set(fp,fn,ROUNDING);

      mpfr_add(xn,xp,dx,ROUNDING);

      if ( mpfr_cmp_d(dx,0.0)>0 && mpfr_greaterequal_p(xn,b) )
      {
	 sec=0;
	 mpfr_set(next->y[i],xp,ROUNDING);
	 break;
      }
      if ( mpfr_cmp_d(dx,0.0)<0 && mpfr_lessequal_p(xn,b) )
      {
	 sec=0;
	 mpfr_set(next->y[i],xp,ROUNDING);
	 break;
      }
      abs_deviation(next,fn,xn);

      ++it;
      if (it>BRACKETIMAX) 
      {
	 sec=0;
	 mpfr_set(next->y[i],xp,ROUNDING);
	 break;
      }
   }
   if(it==1)
   {
      sec=0;
      mpfr_set(next->y[i],xp,ROUNDING);
   }


   mpfr_set(x1,xp,ROUNDING);
   mpfr_set(x2,xn,ROUNDING);   
   mpfr_sub(x0,x1,dx,ROUNDING);

   if (sec)
   {
      abs_deviation(next,fm,x0);
      abs_deviation(next,fx,x1);
      abs_deviation(next,fp,x2);	 

      mpfr_mul(xm,x0,x0,ROUNDING);
      mpfr_mul(x,x1,x1,ROUNDING);
      mpfr_mul(xp,x2,x2,ROUNDING);

      mpfr_sub(cx3,xm,x,ROUNDING);
      mpfr_mul(cx1,fp,cx3,ROUNDING);

      mpfr_sub(cx3,x,xp,ROUNDING);
      mpfr_mul(cx3,cx3,fm,ROUNDING);
      mpfr_add(cx1,cx1,cx3,ROUNDING);

      mpfr_sub(cx3,xp,xm,ROUNDING);
      mpfr_mul(cx3,cx3,fx,ROUNDING);
      mpfr_add(cx1,cx1,cx3,ROUNDING);

      mpfr_set(xm,x0,ROUNDING);
      mpfr_set(x,x1,ROUNDING);
      mpfr_set(xp,x2,ROUNDING);

      mpfr_sub(cx3,xm,x,ROUNDING);
      mpfr_mul(cx2,fp,cx3,ROUNDING);

      mpfr_sub(cx3,x,xp,ROUNDING);
      mpfr_mul(cx3,cx3,fm,ROUNDING);
      mpfr_add(cx2,cx2,cx3,ROUNDING);

      mpfr_sub(cx3,xp,xm,ROUNDING);
      mpfr_mul(cx3,cx3,fx,ROUNDING);
      mpfr_add(cx2,cx2,cx3,ROUNDING);

      mpfr_div(cx3,cx1,cx2,ROUNDING);
      mpfr_mul_d(cx3,cx3,0.5,ROUNDING);

      mpfr_set(next->y[i],cx3,ROUNDING);	    
   }

   mpfr_clears(dx,
	       x0,x1,x2,b,
	       x,xp,xm,xn,
	       cx1,cx2,cx3,
	       fx,fp,fm,fn,
	       NULL);
}


static void second_algorithm_search(remes_t *next)
{
   int i,n;
   mpfr_t M,cx1,cx2;
   mpfr_t check,checks;

   mpfr_inits(M,cx1,cx2,check,checks,NULL);

   n=next->nx;

   for(i=0;i<n;++i)
      bracket(next,i);      

   for (i=0;i<n;++i)
   {
      mpfr_set(cx1,next->x[i],ROUNDING);	 
      mpfr_set(next->x[i],next->y[i],ROUNDING);	 
      mpfr_set(next->y[i],cx1,ROUNDING);	 

      deviation(next,next->e[i],next->x[i]);
      mpfr_set_d(next->se[i],(double)mpfr_sgn(next->e[i]),ROUNDING);      

      if( mpfr_cmp(next->x[i],next->y[i])!=0 )
      {
	 mpfr_sub(next->dx[i],next->x[i],next->y[i],ROUNDING);	 
	 mpfr_abs(next->dx[i],next->dx[i],ROUNDING);	 
	 mpfr_mul(next->dx[i],next->dx[i],next->sdx,ROUNDING);	 
      }
      else
	 mpfr_mul(next->dx[i],next->dx[i],next->sdx,ROUNDING);	 
   }

   search_emax(next,next->nx,M);   

   mpfr_abs(cx1,next->ee,ROUNDING);
   mpfr_sub(check,M,cx1,ROUNDING);
   mpfr_add(checks,M,cx1,ROUNDING);
   mpfr_mul(checks,checks,next->reeps2,ROUNDING);
   mpfr_mul_d(checks,checks,0.5,ROUNDING);

   if( mpfr_less_p(check,next->eeps2) || mpfr_less_p(check,checks) )
   {
      for (i=0;i<n;++i)
      {
	 mpfr_set(next->x[i],next->y[i],ROUNDING);	 
	 deviation(next,next->e[i],next->x[i]);
	 mpfr_set_d(next->se[i],(double)mpfr_sgn(next->e[i]),ROUNDING);      
      }
      for (i=1;i<n;++i)
      {
	 mpfr_mul(cx1,next->se[i],next->se[i-1],ROUNDING);      
	 error( (mpfr_cmp_d(cx1,0.0)>=0) || (mpfr_lessequal_p(next->x[i],next->x[i-1])) ,1,
		"second_algorithm_search [remes.c]","error extrema do not satisfy Chebyshev condition");
      }

      search_emin(next,next->nx,next->emin);
      search_emax(next,next->nx,next->emax);

      next->sflag=-1;
   }
   else
      next->sflag=1;

   mpfr_clears(M,cx1,cx2,check,checks,NULL);
}


static void second_algorithm(remes_t *rem)
{
   int i,n,it,itmax;

   itmax=SECONDIMAX;
   rem->sflag=0;
   n=rem->nx;
   for(i=0;i<n-1;++i)
   {
      mpfr_sub(rem->dx[i],rem->x[i],rem->x[i+1],ROUNDING);
      mpfr_mul(rem->dx[i],rem->dx[i],rem->sdx,ROUNDING);
      if(rem->n<3)
	 for(it=0;it<3;++it)
	    mpfr_mul(rem->dx[i],rem->dx[i],rem->sdx,ROUNDING);

      mpfr_abs(rem->dx[i],rem->dx[i],ROUNDING);
   }
   mpfr_set(rem->dx[n-1],rem->dx[0],ROUNDING);

   it=0;
   while(rem->sflag!=-1)
   {
      second_algorithm_solve(rem);
      second_algorithm_search(rem);
      ++it;
      error(it>itmax,1,"second_algorithm [remes.c]",
	    "remes algorithm did not converge in maximum number of iterations: try to rise the precision");
   }
   search_emax(rem,rem->nx,rem->emax);
}


static void balanc(mpfr_t *a,int n)
{
   int last,j,i;
   mpfr_t t1,t2,s,r,g,f,c,sqrdx,radix;

   mpfr_inits(t1,t2,s,r,g,f,c,sqrdx,radix,NULL);

   mpfr_set_d(radix,2.0,ROUNDING);
   mpfr_mul(sqrdx,radix,radix,ROUNDING);

   last=0;
   while (last == 0) 
   {
      last=1;
      for (i=0;i<n;i++) 
      {
	 mpfr_set_d(r,0.0,ROUNDING);
	 mpfr_set_d(c,0.0,ROUNDING);

	 for (j=0;j<n;j++)
	    if (j != i) 
	    {
	       mpfr_abs(t1,a[j+n*i],ROUNDING);
	       mpfr_abs(t2,a[i+n*j],ROUNDING);

	       mpfr_add(c,c,t1,ROUNDING);
	       mpfr_add(r,r,t2,ROUNDING);
	    }
	 if ( mpfr_get_d(c,ROUNDING) && mpfr_get_d(r,ROUNDING)) 
	 {
	    mpfr_div(g,r,radix,ROUNDING);
	    mpfr_set_d(f,1.0,ROUNDING);

	    mpfr_add(s,c,r,ROUNDING);

	    while (mpfr_less_p(c,g)) 
	    {
	       mpfr_mul(f,f,radix,ROUNDING);
	       mpfr_mul(c,c,sqrdx,ROUNDING);
	    }
	    mpfr_mul(g,r,radix,ROUNDING);

	    while (mpfr_greater_p(c,g)) 
	    {
	       mpfr_div(f,f,radix,ROUNDING);
	       mpfr_div(c,c,sqrdx,ROUNDING);
	    }

	    mpfr_add(t1,c,r,ROUNDING);
	    mpfr_div(t1,t1,f,ROUNDING);

	    mpfr_mul_d(t2,s,0.95,ROUNDING);

	    if (mpfr_less_p(t1,t2)) 
	    {
	       last=0;
	       mpfr_set_d(g,1.0,ROUNDING);
	       mpfr_div(g,g,f,ROUNDING);

	       for (j=0;j<n;j++) 
		  mpfr_mul(a[i+n*j],a[i+n*j],g,ROUNDING);

	       for (j=0;j<n;j++) 
		  mpfr_mul(a[j+n*i],a[j+n*i],f,ROUNDING);
	    }
	 }
      }
   }

   mpfr_clears(t1,t2,s,r,g,f,c,sqrdx,radix,NULL);
}


static void hqr(mpfr_t *a,int n,mpfr_t *roots)
{
   int nn,m,l,k,j,its,i,mmin;
   mpfr_t t1,t2,z,y,x,w,v,u,t,s,r,q,p,anorm;
   
   mpfr_inits(t1,t2,z,y,x,w,v,u,t,s,r,q,p,anorm,NULL);

   mpfr_abs(anorm,a[0],ROUNDING);
   for (i=2;i<=n;i++)
      for (j=(i-1);j<=n;j++)
      {
	 mpfr_abs(t1,a[(i-1)+n*(j-1)],ROUNDING);
	 mpfr_add(anorm,anorm,t1,ROUNDING);
      }

   nn=n;
   mpfr_set_d(t,0.0,ROUNDING);
   while (nn >= 1) 
   {
      its=0;
      do 
      {
	 for (l=nn;l>=2;l--) 
	 {
	    mpfr_abs(t1,a[(l-2)+n*(l-2)],ROUNDING);
	    mpfr_abs(t2,a[(l-1)+n*(l-1)],ROUNDING);
	    mpfr_add(s,t1,t2,ROUNDING);

	    if ( mpfr_cmp_d(s,0.0)==0 ) 
	       mpfr_set(s,anorm,ROUNDING);
	    
	    mpfr_abs(t1,a[(l-1)+n*(l-2)],ROUNDING);
	    mpfr_add(t1,t1,s,ROUNDING);
	    
	    if (mpfr_cmp(t1,s)==0) 
	       break;
	 }
	 mpfr_set(x,a[(nn-1)+n*(nn-1)],ROUNDING);
	 if (l == nn) 
	 {
	    mpfr_add(roots[2*(nn-1)],x,t,ROUNDING);
	    mpfr_set_d(roots[2*(nn-1)+1],0.0,ROUNDING);

	    nn--;
	 } 
	 else 
	 {
	    mpfr_set(y,a[(nn-2)+n*(nn-2)],ROUNDING);
	    mpfr_mul(w,a[(nn-1)+n*(nn-2)],a[(nn-2)+n*(nn-1)],ROUNDING);
	    if (l == (nn-1)) 
	    {
	       mpfr_sub(p,y,x,ROUNDING);
	       mpfr_mul_d(p,p,0.5,ROUNDING);

	       mpfr_mul(q,p,p,ROUNDING);
	       mpfr_add(q,q,w,ROUNDING);

	       mpfr_abs(z,q,ROUNDING);
	       mpfr_sqrt(z,z,ROUNDING);

	       mpfr_add(x,x,t,ROUNDING);

	       if (mpfr_cmp_d(q,0.0)>=0) 
	       {
		  
		  mpfr_abs(t1,z,ROUNDING);
		  if (mpfr_cmp_d(p,0.0)<0) 
		     mpfr_neg(t1,t1,ROUNDING);
		  
		  mpfr_add(z,p,t1,ROUNDING);

		  mpfr_add(roots[2*(nn-1)],x,z,ROUNDING);
		  mpfr_add(roots[2*(nn-2)],x,z,ROUNDING);

		  if (mpfr_get_d(z,ROUNDING))
		  {
		     mpfr_div(t1,w,z,ROUNDING);
		     mpfr_sub(roots[2*(nn-1)],x,t1,ROUNDING);
		  }
		  mpfr_set_d(roots[2*(nn-1)+1],0.0,ROUNDING);
		  mpfr_set_d(roots[2*(nn-2)+1],0.0,ROUNDING);
	       } 
	       else 
	       {
		  mpfr_add(roots[2*(nn-1)],x,p,ROUNDING);
		  mpfr_add(roots[2*(nn-2)],x,p,ROUNDING);

		  mpfr_set(roots[2*(nn-1)+1],z,ROUNDING);
		  mpfr_neg(roots[2*(nn-2)+1],z,ROUNDING);
	       }
	       nn-= 2;
	    } 
	    else 
	    {
	       error(its==HQRIMAX,1,"hqr [remes.c]","hqr did not converge in max num of iterations");
	       if ( (its+10)%10 == 0 ) 
	       {
		  mpfr_add(t,t,x,ROUNDING);
		  for (i=1;i<=nn;i++) 
		     mpfr_sub(a[(i-1)+n*(i-1)],a[(i-1)+n*(i-1)],x,ROUNDING);

		  mpfr_abs(t1,a[(nn-1)+n*(nn-2)],ROUNDING);
		  mpfr_abs(t2,a[(nn-2)+n*(nn-3)],ROUNDING);

		  mpfr_add(s,t1,t2,ROUNDING);

		  mpfr_mul_d(x,s,0.75,ROUNDING);
		  mpfr_mul_d(y,s,0.75,ROUNDING);

		  mpfr_mul(w,s,s,ROUNDING);
		  mpfr_mul_d(w,w,-0.4375,ROUNDING);
	       }
	       ++its;
	       for (m=(nn-2);m>=l;m--) 
	       {
		  mpfr_set(z,a[(m-1)+n*(m-1)],ROUNDING);

		  mpfr_sub(r,x,z,ROUNDING);
		  mpfr_sub(s,y,z,ROUNDING);

		  mpfr_mul(p,r,s,ROUNDING);
		  mpfr_sub(p,p,w,ROUNDING);
		  mpfr_div(p,p,a[m+n*(m-1)],ROUNDING);
		  mpfr_add(p,p,a[(m-1)+n*m],ROUNDING);

		  mpfr_sub(q,a[m+n*m],z,ROUNDING);
		  mpfr_sub(q,q,r,ROUNDING);
		  mpfr_sub(q,q,s,ROUNDING);

		  mpfr_set(r,a[(m+1)+n*m],ROUNDING);

		  mpfr_abs(t1,p,ROUNDING);
		  mpfr_abs(t2,q,ROUNDING);
		  mpfr_add(s,t1,t2,ROUNDING);
		  mpfr_abs(t2,r,ROUNDING);
		  mpfr_add(s,s,t2,ROUNDING);

		  mpfr_div(p,p,s,ROUNDING);
		  mpfr_div(q,q,s,ROUNDING);
		  mpfr_div(r,r,s,ROUNDING);

		  if (m == l) 
		     break;

		  mpfr_abs(t1,q,ROUNDING);
		  mpfr_abs(t2,r,ROUNDING);
		  mpfr_add(u,t1,t2,ROUNDING);
		  mpfr_abs(t2,a[(m-1)+n*(m-2)],ROUNDING);
		  mpfr_mul(u,u,t2,ROUNDING);

		  mpfr_abs(t1,a[m+n*m],ROUNDING);
		  mpfr_abs(t2,a[(m-2)+n*(m-2)],ROUNDING);
		  mpfr_add(v,t1,t2,ROUNDING);
		  mpfr_abs(t2,z,ROUNDING);
		  mpfr_add(v,v,t2,ROUNDING);
		  mpfr_abs(t2,p,ROUNDING);
		  mpfr_mul(v,v,t2,ROUNDING);

		  mpfr_add(t2,u,v,ROUNDING);

		  if (mpfr_cmp(t2,v)==0) 
		     break;
	       }
	       for (i=m+2;i<=nn;i++) 
	       {
		  mpfr_set_d(a[(i-1)+n*(i-3)],0.0,ROUNDING);
		  if (i != (m+2)) 
		     mpfr_set_d(a[(i-1)+n*(i-4)],0.0,ROUNDING);
	       }
	       for (k=m;k<=nn-1;k++) 
	       {
		  if (k != m) 
		  {
		     mpfr_set(p,a[(k-1)+n*(k-2)],ROUNDING);
		     mpfr_set(q,a[k+n*(k-2)],ROUNDING);
		     mpfr_set_d(r,0.0,ROUNDING);

		     if (k != (nn-1)) 
			mpfr_set(r,a[(k+1)+n*(k-2)],ROUNDING);

		     mpfr_abs(t1,r,ROUNDING);
		     mpfr_abs(t2,q,ROUNDING);
		     mpfr_add(x,t1,t2,ROUNDING);
		     mpfr_abs(t2,p,ROUNDING);
		     mpfr_add(x,x,t2,ROUNDING);

		     if (mpfr_cmp_d(x,0.0)!=0) 
		     {
			mpfr_div(p,p,x,ROUNDING);
			mpfr_div(q,q,x,ROUNDING);
			mpfr_div(r,r,x,ROUNDING);
		     }
		  }

		  mpfr_mul(t1,p,p,ROUNDING);
		  mpfr_mul(t2,q,q,ROUNDING);
		  mpfr_add(s,t1,t2,ROUNDING);
		  mpfr_mul(t2,r,r,ROUNDING);
		  mpfr_add(s,s,t2,ROUNDING);
		  mpfr_sqrt(s,s,ROUNDING);

		  if (mpfr_cmp_d(p,0.0)<0) 
		     mpfr_neg(s,s,ROUNDING);

		  if (mpfr_cmp_d(s,0.0)!=0) 
		  {
		     if (k == m) 
		     {
			if (l != m)
			   mpfr_neg(a[(k-1)+n*(k-2)],a[(k-1)+n*(k-2)],ROUNDING);
		     } 
		     else
		     {
			mpfr_mul(a[(k-1)+n*(k-2)],s,x,ROUNDING);
			mpfr_neg(a[(k-1)+n*(k-2)],a[(k-1)+n*(k-2)],ROUNDING);
		     }

		     mpfr_add(p,p,s,ROUNDING);
		     mpfr_div(x,p,s,ROUNDING);
		     mpfr_div(y,q,s,ROUNDING);
		     mpfr_div(z,r,s,ROUNDING);
		     mpfr_div(q,q,p,ROUNDING);
		     mpfr_div(r,r,p,ROUNDING);

		     for (j=k;j<=nn;j++) 
		     {
			mpfr_mul(p,q,a[k+n*(j-1)],ROUNDING);
			mpfr_add(p,p,a[(k-1)+n*(j-1)],ROUNDING);

			if (k != (nn-1)) 
			{
			   mpfr_mul(t1,r,a[(k+1)+n*(j-1)],ROUNDING);
			   mpfr_add(p,p,t1,ROUNDING);

			   mpfr_mul(t1,p,z,ROUNDING);
			   mpfr_sub(a[(k+1)+n*(j-1)],a[(k+1)+n*(j-1)],t1,ROUNDING);
			}
			mpfr_mul(t1,p,y,ROUNDING);
			mpfr_sub(a[k+n*(j-1)],a[k+n*(j-1)],t1,ROUNDING);

			mpfr_mul(t1,p,x,ROUNDING);
			mpfr_sub(a[(k-1)+n*(j-1)],a[(k-1)+n*(j-1)],t1,ROUNDING);
		     }
		     mmin = nn<k+3 ? nn : k+3;
		     for (i=l;i<=mmin;i++) 
		     {
			mpfr_mul(t1,y,a[(i-1)+n*k],ROUNDING);
			mpfr_mul(t2,x,a[(i-1)+n*(k-1)],ROUNDING);

			mpfr_add(p,t1,t2,ROUNDING);

			if (k != (nn-1)) 
			{
			   mpfr_mul(t1,z,a[(i-1)+n*(k+1)],ROUNDING);
			   mpfr_add(p,p,t1,ROUNDING);

			   mpfr_mul(t1,p,r,ROUNDING);
			   mpfr_sub(a[(i-1)+n*(k+1)],a[(i-1)+n*(k+1)],t1,ROUNDING);
			}
			mpfr_mul(t1,p,q,ROUNDING);
			mpfr_sub(a[(i-1)+n*k],a[(i-1)+n*k],t1,ROUNDING);

			mpfr_sub(a[(i-1)+n*(k-1)],a[(i-1)+n*(k-1)],p,ROUNDING);
		     }
		  }
	       }
	    }
	 }
      } while (l < nn-1);
   }

   mpfr_clears(t1,t2,z,y,x,w,v,u,t,s,r,q,p,anorm,NULL);
}


static void zrhqr(mpfr_t a[],int m,mpfr_t roots[])
{
   int j,k;
   mpfr_t *hess,xr,xi,t1,t2;

   hess=malloc(m*m*sizeof(mpfr_t));
   error(hess==NULL,1,"zrhqr [remes.c]","Unable to allocate auxiliary array");

   mpfr_inits(xr,xi,t1,t2,NULL);
   for (k=0;k<m*m;k++) 
      mpfr_init(hess[k]);

   for (k=0;k<m;k++) 
   {
      mpfr_div(hess[m*k],a[m-1-k],a[m],ROUNDING);
      mpfr_neg(hess[m*k],hess[m*k],ROUNDING);

      for (j=1;j<m;j++) 
	 mpfr_set_d(hess[j+m*k],0.0,ROUNDING);
      if (k != m-1) 
	 mpfr_set_d(hess[k+1+m*k],1.0,ROUNDING);
   }

   balanc(hess,m);
   hqr(hess,m,roots);

   for (j=1;j<m;j++) 
   {
      mpfr_set(xr,roots[2*j],ROUNDING);
      mpfr_set(xi,roots[2*j+1],ROUNDING);
      for (k=j-1;k>=0;k--) 
      {
	 mpfr_abs(t1,roots[2*j],ROUNDING);
	 mpfr_abs(t2,xr,ROUNDING);

	 if ( mpfr_greaterequal_p(t1,t2)) 
	    break;
	 mpfr_set(roots[2*(k+1)],roots[2*k],ROUNDING);
	 mpfr_set(roots[2*(k+1)+1],roots[2*k+1],ROUNDING);
      }
      mpfr_set(roots[2*(k+1)],xr,ROUNDING);
      mpfr_set(roots[2*(k+1)+1],xi,ROUNDING);
   }

   mpfr_clears(xr,xi,t1,t2,NULL);
   for (k=0;k<m*m;k++) 
      mpfr_clear(hess[k]);
   
   free(hess);
}



static void change_root_basis(remes_t *rem)
{
   int i;
   mpfr_t r;

   mpfr_init(r);

   for(i=0;i<rem->n;++i)
   {
      mpfr_mul(rem->rn[2*i],rem->rn[2*i],rem->bma,ROUNDING);
      mpfr_add(rem->rn[2*i],rem->rn[2*i],rem->bpa,ROUNDING);
      mpfr_sub(rem->rn[2*i],rem->rn[2*i],rem->mu2,ROUNDING);
      mpfr_neg(rem->rn[2*i],rem->rn[2*i],ROUNDING);

      mpfr_abs(r,rem->rn[2*i+1],ROUNDING);
      if(mpfr_cmp_d(r,1.0e-16)>0)
      {
	 mpfr_mul(rem->rn[2*i+1],rem->rn[2*i+1],rem->bma,ROUNDING);
	 mpfr_neg(rem->rn[2*i+1],rem->rn[2*i+1],ROUNDING);
      }

      mpfr_div(rem->cn[rem->n],rem->cn[rem->n],rem->bma,ROUNDING);
   }

   for(i=0;i<rem->d;++i)
   {
      mpfr_mul(rem->rd[2*i],rem->rd[2*i],rem->bma,ROUNDING);
      mpfr_add(rem->rd[2*i],rem->rd[2*i],rem->bpa,ROUNDING);
      mpfr_sub(rem->rd[2*i],rem->rd[2*i],rem->mu2,ROUNDING);
      mpfr_neg(rem->rd[2*i],rem->rd[2*i],ROUNDING);

      mpfr_abs(r,rem->rd[2*i+1],ROUNDING);
      if(mpfr_cmp_d(r,1.0e-16)>0)
      {
	 mpfr_mul(rem->rd[2*i+1],rem->rd[2*i+1],rem->bma,ROUNDING);
	 mpfr_neg(rem->rd[2*i+1],rem->rd[2*i+1],ROUNDING);
      }

      mpfr_div(rem->cd[rem->d],rem->cd[rem->d],rem->bma,ROUNDING);
   }

   for(i=0;i<rem->nd2;++i)
   {
      mpfr_mul(rem->x[i],rem->x[i],rem->bma,ROUNDING);
      mpfr_add(rem->x[i],rem->x[i],rem->bpa,ROUNDING);
      mpfr_sub(rem->x[i],rem->x[i],rem->mu2,ROUNDING);
   }

   mpfr_clear(r);
}


static void find_roots(remes_t *rem)
{
   zrhqr(rem->cn,rem->n,rem->rn);
   zrhqr(rem->cd,rem->d,rem->rd);
   change_root_basis(rem);
}


static void print_logfile_header(remes_t *rem)
{
   mpfr_t ra2,rb2;

   mpfr_inits(ra2,rb2,NULL);

   mpfr_mul(ra2,rem->ra,rem->ra,ROUNDING);
   mpfr_mul(rb2,rem->rb,rem->rb,ROUNDING);

   printf("\n\n");
   printf("MinMax rational approximation\n\n");
   printf("\t R^n_{ra rb}(x) = A * { ( x + a_1 )*( x + a_3 ) ... ( x + a_2n-1 ) }\n");
   printf("\t                    / { ( x + a_2 )*( x + a_4 ) ... ( x + a_2n   ) }\n\n");
   printf("of the function\n\n");
   printf("\t f(x)= ( x + mu^2 )^{%d/%d},\t mu= %.8e,\n\n",
	  rem->p,rem->q,
	  mpfr_get_d(rem->mu,ROUNDING)
      );
   printf("\t x in [ra^2,rb^2]= [%.8e,%.8e].\n\n",
	  mpfr_get_d(ra2,ROUNDING),
	  mpfr_get_d(rb2,ROUNDING)
      );
   printf("The program minimizes the error ");
   if(rem->relflag==1)
   {
      printf("\n\n");
      printf("\t delta= max_[ra^2,rb^2]{ |delta(x)| },\t delta(x)= 1 - R^n_{ra rb}(x)/f(x).\n\n");
   }
   else
   {
      printf("(the -abs option has been set)\n\n");
      printf("\t delta= max_[ra^2,rb^2]{ |delta(x)| },\t delta(x)= f(x) - R^n_{ra rb}(x).\n\n");
   }

   printf("The MinMax approximation of degree n is obtained by searching for the nodes,\n\n");
   printf("\t x_i,  i= 0,...,2n+1,\t e_i= delta(x_i).\n\n");
   printf("These are 2n+2 points such that\n\n");
   printf("\t x_{i+1} > x_i,\t e_i = (-1)^i e_0,\n\n");
   printf("and the approximation is such that\n\n"); 
   printf("\t delta <= |e_0|.\n\n");

   if(rem->nstart>0)
      printf("The iteration starts at n= %d and stops when\n\n",rem->nstart);
   else
      printf("The iteration starts at n= %d and stops when\n\n",1);
   printf("\t delta<= %.8e,\t n= n_max\n\n",mpfr_get_d(rem->goal,ROUNDING));
   printf("The arithmetic accuracy corresponds to floating point numbers with %d bits\n\n",rem->precision);
   printf("\n\n\n\n");

   fflush(stdout);

   mpfr_clears(ra2,rb2,NULL);
}


static void print_info(remes_t *rem)
{
   int i;
   mpfr_t A,re;

   mpfr_inits(A,re,NULL);

   for(i=0;i<80;++i)
      message("-");
   message("\n");

   message("n= %d error <= %.1e\n",rem->n,mpfr_get_d(rem->emax,ROUNDING));

   for(i=0;i<80;++i)
      message("-");
   message("\n\n");

   mpfr_div(A,rem->cn[rem->n],rem->cd[rem->d],ROUNDING);
   message("A\t %24.16e\n\n",mpfr_get_d(A,ROUNDING));
   for(i=0;i<rem->n;++i)
   {
      message("a_%d\t ( %24.16e, %24.16e )\n",2*i+1,mpfr_get_d(rem->rn[2*i],ROUNDING),mpfr_get_d(rem->rn[2*i+1],ROUNDING));
      message("a_%d\t ( %24.16e, %24.16e )\n",2*i+2,mpfr_get_d(rem->rd[2*i],ROUNDING),mpfr_get_d(rem->rd[2*i+1],ROUNDING));
   }
   message("\n\n");

   for(i=0;i<rem->nd2;++i)
   {
      rdeviation(rem,re,rem->x[i]);
      mpfr_sub(A,rem->e[i],re,ROUNDING);
      mpfr_abs(A,A,ROUNDING);

      error(mpfr_cmp_d(A,1.0e-13)>0,1,"remes [remes.c]","Precision lost in change of basis");      

      message("x_%d\t %24.16e    e_%d\t %24.16e\n",
	      i,mpfr_get_d(rem->x[i],ROUNDING),
	      i,mpfr_get_d(re,ROUNDING)
	 );
   }
   message("\n");


   for(i=0;i<80;++i)
      message("-");
   message("\n\n\n\n\n");
   
   fflush(stdout);

   mpfr_clears(A,re,NULL);
}


static void print_openQcdQed_parameters(remes_t *rem)
{
   int i,nflag,cflag;
   char ologname[NAME_SIZE];
   FILE *olog=NULL;
   mpfr_t A,re,nu,mu;

   mpfr_inits(A,re,nu,mu,NULL);

   nflag=0;
   cflag=0;
   for(i=0;i<rem->n;++i)
   {
      mpfr_abs(nu,rem->rn[2*i+1],ROUNDING);
      mpfr_abs(mu,rem->rd[2*i+1],ROUNDING);

      if(mpfr_cmp_d(nu,1.0e-16)>0)
	 cflag=1;

      if(mpfr_cmp_d(mu,1.0e-16)>0)
	 cflag=1;

      if(mpfr_cmp_d(rem->rn[2*i],0.0)<0)
	 nflag=1;

      if(mpfr_cmp_d(rem->rd[2*i],0.0)<0)
	 nflag=1;
   }

   sprintf(ologname,"%s/n%d.in",rem->path,rem->n);
   olog=fopen(ologname,"w");

   if(nflag==0 && cflag==0)
   {
      fprintf(olog,"power\t\t %d %d\n",rem->p,rem->q);
      
      fprintf(olog,"degree\t\t %d\n",rem->n);
      
      fprintf(olog,"range\t\t %.8e %.8e\n",
	      mpfr_get_d(rem->ra,ROUNDING),
	      mpfr_get_d(rem->rb,ROUNDING)
	 );

      fprintf(olog,"mu\t\t %.8e\n",mpfr_get_d(rem->mu,ROUNDING));
      
      rdeviation(rem,re,rem->x[0]);
      fprintf(olog,"delta\t\t %.16e\n",mpfr_get_d(re,ROUNDING));
      
      mpfr_div(A,rem->cn[rem->n],rem->cd[rem->d],ROUNDING);
      fprintf(olog,"A\t\t %.20e\n",mpfr_get_d(A,ROUNDING));
      for(i=0;i<rem->n;++i)
      {
	 mpfr_sqrt(nu,rem->rn[2*i],ROUNDING);
	 mpfr_sqrt(mu,rem->rd[2*i],ROUNDING);
	 fprintf(olog,"nu[%d]\t\t %.20e\n",i,mpfr_get_d(nu,ROUNDING));
	 fprintf(olog,"mu[%d]\t\t %.20e\n",i,mpfr_get_d(mu,ROUNDING));
      }
      
      for(i=0;i<rem->nd2;++i)
	 fprintf(olog,"x[%d]\t\t %.20e\n",i,mpfr_get_d(rem->x[i],ROUNDING));
   }
   else
   {
      fprintf(olog,"nflag= %2d  cflag= %2d\n",nflag,cflag);
      fprintf(olog,"some of the roots are negative or complex, the approximation cannot be used in openQcdQed!\n");
   }

   fflush(olog);
   fclose(olog);      

   mpfr_clears(A,re,nu,mu,NULL);
}


void remes(remes_t *rem)
{ 
   int i,done;
   double power;
   mpfr_t f;

   error(rem->p==0 || rem->q<0,1,"remes [remes.c]","Parameters p,q have not been set properly");
   error(mpfr_lessequal_p(rem->rb,rem->ra),1,"remes [remes.c]","Parameters ra,rb,mu have not been set properly");
   error(rem->relflag!=0 && rem->relflag!=1,1,"remes [remes.c]","Parameter relflag has not been set properly");
   error(rem->verbose!=0 && rem->verbose!=1,1,"remes [remes.c]","Parameter verbose has not been set properly");

   mpfr_init(f);

   power=(double)(rem->p);
   power/=(double)(rem->q);
   mpfr_set_d(rem->power,power,ROUNDING);

   mpfr_mul(rem->mu2,rem->mu,rem->mu,ROUNDING);

   mpfr_mul(rem->a,rem->ra,rem->ra,ROUNDING);
   mpfr_mul(rem->b,rem->rb,rem->rb,ROUNDING);

   mpfr_add(rem->a,rem->a,rem->mu2,ROUNDING);
   mpfr_add(rem->b,rem->b,rem->mu2,ROUNDING);

   error(mpfr_cmp_d(rem->a,0.0)==0,1,"remes [remes.c]","Parameters ra,rb,mu have not been set properly");


   mpfr_add(rem->bpa,rem->b,rem->a,ROUNDING);
   mpfr_sub(rem->bma,rem->b,rem->a,ROUNDING);
   mpfr_mul_d(rem->bpa,rem->bpa,0.5,ROUNDING);
   mpfr_mul_d(rem->bma,rem->bma,0.5,ROUNDING);

   mpfr_set_d(rem->sdx,REMESSDX,ROUNDING);

   mpfr_set_d(rem->eeps,REMESEPS,ROUNDING);
   mpfr_set_d(rem->eeps2,REMESEPS,ROUNDING);
   mpfr_set_d(rem->reeps,REMESREPS,ROUNDING);
   mpfr_set_d(rem->reeps2,REMESREPS,ROUNDING);

   rem->n=0;
   rem->d=0;
   rem->nd=rem->n+rem->d;
   rem->nd2=rem->nd+2;

   rem->nx=2;

   mpfr_set_d(rem->x[0],-1.0,ROUNDING);
   mpfr_set_d(rem->x[1],1.0,ROUNDING);

   function(rem,f,rem->x[0]);
   mpfr_set(rem->cn[0],f,ROUNDING);

   function(rem,f,rem->x[1]);
   mpfr_add(rem->cn[0],rem->cn[0],f,ROUNDING);
   mpfr_mul_d(rem->cn[0],rem->cn[0],0.5,ROUNDING);

   mpfr_set_d(rem->cd[0],1.0,ROUNDING);

   deviation(rem,rem->e[0],rem->x[0]);
   deviation(rem,rem->e[1],rem->x[1]);

   mpfr_set_d(rem->se[0],mpfr_sgn(rem->e[0]),ROUNDING);
   mpfr_set_d(rem->se[1],mpfr_sgn(rem->e[1]),ROUNDING);

   print_logfile_header(rem);

   for(i=1;i<rem->nstart;++i)
      starting_algorithm_cheby(rem);

   done=0;
   while(done==0)
   {
      starting_algorithm_cheby(rem);

      second_algorithm(rem);
      
      find_roots(rem);   

      if(rem->n>1)
	 set_nfit(rem);

      print_info(rem);   

      if(rem->verbose>0)
	 print_openQcdQed_parameters(rem);

      if ( mpfr_lessequal_p(rem->emax,rem->goal) )
	 done=1;
   }

   mpfr_clear(f);
}
