
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2017, 2019 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Analytic check of RWRTM reweighting factor.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "archive.h"
#include "mdflds.h"
#include "sflds.h"
#include "linalg.h"
#include "dirac.h"
#include "sap.h"
#include "dfl.h"
#include "forces.h"
#include "update.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static int my_rank;


#ifndef REWEIGHT_CHECK_PROGRAMS
#error This program must be compiled with the REWEIGHT_CHECK_PROGRAMS macro
#endif

static complex_dble dhat;


/*******************************************************************************
Dwhat(mu) = dhat.re + i ( dhat.im g0 + mu g5 )
Dwhat(mu)^dag = dhat.re - i ( dhat.im g0 + mu g5 ) = g5 Dwhat(-mu) g5
Dwhat(mu)^dag Dwhat(mu) = dhat.re^2 + ( dhat.im g0 + mu g5 )^2 =
                        = dhat.re^2 + dhat.im^2 + mu^2
*******************************************************************************/
void Dwhat_dble_rcp(double mu,spinor_dble *s,spinor_dble *r)
{
   spinor_dble *sm;

   sm=s+VOLUME/2;

   for (;s<sm;s++)
   {
      _vector_mul((*r).c1,dhat.re,(*s).c1);
      _vector_mul((*r).c2,dhat.re,(*s).c2);
      _vector_mul((*r).c3,dhat.re,(*s).c3);
      _vector_mul((*r).c4,dhat.re,(*s).c4);
      
      _vector_mulir_assign((*r).c1,-dhat.im,(*s).c3);
      _vector_mulir_assign((*r).c2,-dhat.im,(*s).c4);
      _vector_mulir_assign((*r).c3,-dhat.im,(*s).c1);
      _vector_mulir_assign((*r).c4,-dhat.im,(*s).c2);

      _vector_mulir_assign((*r).c1,+mu,(*s).c1);
      _vector_mulir_assign((*r).c2,+mu,(*s).c2);
      _vector_mulir_assign((*r).c3,-mu,(*s).c3);
      _vector_mulir_assign((*r).c4,-mu,(*s).c4);

      r+=1;
   }   
}


/*******************************************************************************
Multi-shift CG solver for the normal even-odd preconditioned Wilson-Dirac
equation (Dwhat^dag*Dwhat+mu^2)*psi=eta with a twisted-mass term.

(Dwhat^dag*Dwhat+mu^2)^{-1} =
  = [Dwhat(mu)^dag Dwhat(mu)]^{-1} =
  = [dhat.re^2 + dhat.im^2 + mu^2]^{-1}
*******************************************************************************/
void tmcgm_rcp(int nmx,double *res,int nmu,double *mu,
                spinor_dble *eta,spinor_dble **psi,int *status)
{
   int k;
   double sqn;
   
   sqn=dhat.re*dhat.re+dhat.im*dhat.im;
   
   for (k=0;k<nmu;k++)
   {
      set_sd2zero(VOLUME/2,psi[k]);
      mulr_spinor_add_dble(VOLUME/2,psi[k],eta,1.0/(sqn+mu[k]*mu[k]));
   }
   
   status[0]=1;
}


void tmcgeo_rcp(int nmx,double res,double mu,
                spinor_dble *eta,spinor_dble *psi,int *status)
{
   double sqn;
   
   sqn=dhat.re*dhat.re+dhat.im*dhat.im;
   
   set_sd2zero(VOLUME/2,psi);
   mulr_spinor_add_dble(VOLUME/2,psi,eta,1.0/(sqn+mu*mu));
   
   status[0]=1;
}


/*******************************************************************************
Since the twisted-mass flag is set, the sap_gcr program solves the equation
  (Dw+i*mu*gamma_5*1e)*psi=eta
This equation is equivalent to the set of two equations
  (Dwhat+i*mu*gamma_5*1e)*psi(e) = eta(e) - Deo*Doo^{-1}*eta(o)
  psi(o) = Doo^{-1}*(eta(o) - Doe*psi(e))

If eta(o)=0, as it is set it the rw* programs, then these equations become
  (Dwhat+i*mu*gamma_5*1e)*psi(e) = eta(e)
  psi(o) = - Doo^{-1}*Doe*psi(e)
i.e. if we discard psi(o) then the even-odd preconditioned Dirac equation is
solved.

Assume restriction to the even sites.

(Dwhat+i*mu*gamma_5)^{-1} = Dwhat(mu)^{-1}
  = Dwhat(mu)^dag [ Dwhat(mu) Dwhat(mu)^dag ]^{-1} =
  = [ dhat.re - i ( dhat.im g0 + mu g5 ) ] / [ dhat.re^2 + dhat.im^2 + mu^2 ]
*******************************************************************************/
double sap_gcr_rcp(int nkv,int nmx,double res,double mu,
                    spinor_dble *eta,spinor_dble *psi,int *status)
{
   spinor_dble *s,*sm,*r;
   double sqn;
   
   sqn=dhat.re*dhat.re+dhat.im*dhat.im+mu*mu;

   r=psi;
   s=eta;
   sm=s+VOLUME/2;

   for (;s<sm;s++)
   {
      _vector_mul((*r).c1,dhat.re/sqn,(*s).c1);
      _vector_mul((*r).c2,dhat.re/sqn,(*s).c2);
      _vector_mul((*r).c3,dhat.re/sqn,(*s).c3);
      _vector_mul((*r).c4,dhat.re/sqn,(*s).c4);
      
      _vector_mulir_assign((*r).c1,dhat.im/sqn,(*s).c3);
      _vector_mulir_assign((*r).c2,dhat.im/sqn,(*s).c4);
      _vector_mulir_assign((*r).c3,dhat.im/sqn,(*s).c1);
      _vector_mulir_assign((*r).c4,dhat.im/sqn,(*s).c2);

      _vector_mulir_assign((*r).c1,-mu/sqn,(*s).c1);
      _vector_mulir_assign((*r).c2,-mu/sqn,(*s).c2);
      _vector_mulir_assign((*r).c3,+mu/sqn,(*s).c3);
      _vector_mulir_assign((*r).c4,+mu/sqn,(*s).c4);

      r+=1;
   }   
   
   status[0]=1;
   
   return 0.0;
}


double dfl_sap_gcr2_rcp(int idfl,int nkv,int nmx,double res,double mu,
                               spinor_dble *eta,spinor_dble *psi,int *status)
{
   sap_gcr_rcp(nkv,nmx,res,mu,eta,psi,status);
   status[1]=1;
   status[2]=1;
   return 0.0;
}


static void read_rw_factors(void)
{
   int k,np;
   rw_parms_t rwp;
   rat_parms_t rp[2];

   read_rw_parms(0);
   rwp=rw_parms(0);

   error_root(rwp.rwfact!=RWRTM,1,"read_rw_factors [check2.c]",
              "This program works only for the RWRTM reweighting factor");

   rp[0]=rat_parms(rwp.irp[0]);
   rp[1]=rat_parms(rwp.irp[1]);

   if (rp[0].degree==0)
   {
      read_rat_parms(rwp.irp[0]);
      rp[0]=rat_parms(rwp.irp[0]);
   }
   if (rp[1].degree==0)
   {
      read_rat_parms(rwp.irp[1]);
      rp[1]=rat_parms(rwp.irp[1]);
   }

   error_root(rp[0].degree!=rp[1].degree,1,"read_rw_factors [ms1.c]",
              "Rational approximations in RWRTM reweighting factor "
              "must have the same order");

   error_root((rp[0].power[0]!=rp[1].power[0])||
              (rp[0].power[1]!=rp[1].power[1]),1,
              "read_rw_factors [ms1.c]",
              "Rational approximations in RWRTM reweighting factor "
              "must have the same power");

   error_root(rp[1].mu0!=0.0,1,"read_rw_factors [ms1.c]",
              "Rational approximations in RWRTM reweighting factor "
              "must have zero twisted mass");

   np=rp[0].degree;
   for (k=0;k<rwp.nfct;k++)
      np-=rwp.np[k];

   error_root(np!=0,1,"read_rw_factors [ms1.c]",
              "Invalid pole decomposition in reweighting factor (%d,%d,%d)",np,rp[0].degree,rwp.nfct);
}


static void read_solvers(int *isap, int *idfl)
{
   int bs[4]={4,4,4,4};
   int nfct,ifct,isp;
   rw_parms_t rwp;
   solver_parms_t sp;

   (*isap)=0;
   (*idfl)=0;

   rwp=rw_parms(0);
   nfct=rwp.nfct;

   for (ifct=0;ifct<nfct;ifct++)
   {
      isp=rwp.isp[ifct];
      sp=solver_parms(isp);

      if (sp.solver==SOLVERS)
      {
         read_solver_parms(isp);
         sp=solver_parms(isp);

         if (sp.solver==SAP_GCR)
            (*isap)=1;
         else if (sp.solver==DFL_SAP_GCR)
         {
            (*isap)=1;
            (*idfl)=1;
         }
      }
   }

   if (*isap)
      set_sap_parms(bs,1,4,5);
}


double calc_lnr(int np1,int np2,double sqn)
{
   int k;
   double q2,r;
   rw_parms_t rwp;
   rat_parms_t rp1;
   rat_parms_t rp2;
   
   rwp=rw_parms(0);
   rp1=rat_parms(rwp.irp[0]);
   rp2=rat_parms(rwp.irp[1]);
   
   q2=dhat.re*dhat.re+dhat.im*dhat.im;
   
   r=1.0;
   for (k=np1;k<=np2;k++)
   {
      r*=(q2+rp2.nu[k]*rp2.nu[k])/(q2+rp1.nu[k]*rp1.nu[k]);
      r*=(q2+rp1.mu[k]*rp1.mu[k])/(q2+rp2.mu[k]*rp2.mu[k]);
   }
   
   return sqn*(r-1.0);
}


int main(int argc,char *argv[])
{
   int level,seed;
   int isap=0,idfl=0,ifct;
   int nwsd;
   int isrc,np1,np2;
   int status[3];
   double lnr1,lnr2,d,dmax,sqn,r[2],pi;
   rw_parms_t rwp;
   rat_parms_t rp;
   FILE *flog=NULL,*fin=NULL;
   
   pi=4.0*atan(1.0);

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Analytic check of RWRTM reweighting factor\n");
      printf("------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
      printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d process block size\n\n",
             NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);
   }

   if (my_rank==0)
   {
      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);
   }

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

   set_flds_parms(1,0);

   read_rw_factors();
   read_solvers(&isap,&idfl);

   if (my_rank==0)
   {
      fclose(fin);

      printf("Random number generator:\n");
      printf("level = %d, seed = %d\n\n",level,seed);
   }

   print_rw_parms();
   print_rat_parms();
   print_solver_parms(&isap,&idfl);

   start_ranlux(level,seed);
   geometry();


   rwp=rw_parms(0);
   rp=rat_parms(rwp.irp[0]);

   nwsd=4;
   for (ifct=0;ifct<rwp.nfct;ifct++)
   {
      if(3+rwp.np[ifct]>nwsd)
         nwsd=3+rwp.np[ifct];
   }
   alloc_wsd(nwsd);
   
   random_ud();
   
   dmax=0.0;
   for (isrc=0;isrc<rwp.nsrc;isrc++)
   {
      if (my_rank==0)
         ranlxd(r,2);
      MPI_Bcast(r,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
      dhat.re=(rp.range[0]+r[0]*(rp.range[1]-rp.range[0]))*cos(2.0*pi*r[1]);
      dhat.im=(rp.range[0]+r[0]*(rp.range[1]-rp.range[0]))*sin(2.0*pi*r[1]);
      
      lnr1=0.0;
      lnr2=0.0;

      np1=0;
      for (ifct=0;ifct<rwp.nfct;ifct++)
      {
         np2=np1+rwp.np[ifct]-1;
         lnr1+=rwrtm(rwp.irp[0],rwp.irp[1],np1,np2,rwp.isp[ifct],&sqn,status);
         lnr2+=calc_lnr(np1,np2,sqn);

         np1=np2+1;
      }

      d=fabs((lnr1-lnr2)/lnr1);
      if (d>dmax)
         dmax=d;
      
      if (my_rank==0)
      {
         printf("isrc= %3d  |dhat|= %+.2e  lnr1= %+.3e  lnr2= %+.3e  rdev= %.2e \n",isrc,sqrt(dhat.re*dhat.re+dhat.im*dhat.im),lnr1,lnr2,d);
         
         fflush(flog);
      }
   }

   if (my_rank==0)
   {
      printf("\n========================================"
             "=========================================\n\n"
             "Maximum relative deviation (expected < 1e-14) = %.2e\n\n",dmax);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
