
/*******************************************************************************
*
* File rwrtm.c
*
* Copyright (C) 2017, 2019 Agostino Patella
*
* Based on openQCD-1.6/modules/update/rwrat.c
*      and openQCD-1.6/modules/update/rwtmeo.c
* Copyright (C) 2012-2014 Martin Luescher, Stefan Schaefer
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* RHMC twisted-mass reweighting factors.
*
* The externally accessible functions are
*
*   double rwrtm(int irp1,int irp2,int np1,int np2,int isp,
*                double *sqn,int *status)
*     Generates a random pseudo-fermion field with normal distribution,
*     assigns its square norm to sqn and returns -ln(r) (see the notes).
*
* Notes:
*
* The computation of the reweighting factor needed to correct for the fictitious
* twisted mass in the RHMC action is discussed in the notes "RHMC alogorithm in
* openQ*D" (doc/rhmc.pdf).
*
* R is a rational approximation for (Dwhat^dag*Dwhat+muhat^2)^a of order n, and
* R' is a rational approximation for (Dwhat^dag*Dwhat)^a of the same order n.
* Then the operator R^{-1} R' is a rational function of Dwhat^dag*Dwhat of order
* 2n, and is decomposed in factors P(np1,np2) constructed by taking
* 
*  P(np1,np2) = Sum{
*     (Dwhat^dag*Dwhat+mu[j]^2)/(Dwhat^dag*Dwhat+nu[j]^2)
*     *(Dwhat^dag*Dwhat+nu'[j]^2)/(Dwhat^dag*Dwhat+mu'[j]^2), j=np1..np2}
*
* mu[j] and nu[j] being the poles and zeros of R, and mu'[j] and nu'[j] being
* the poles and zeros of R'.
*
* The contribution of P(np1,np2) to the reweighting factor is given by
*
*  R(np1,np2) = det P(np1,np2)^{-1}
*
* For more details, refer to doc/rhmc.pdf.
*
* The function rwrtm return an unbiased stochastic estimate r of the reweighting
* factor R(np1,np2) given by
*
*  -ln(r) = (eta, [P(np1,np2)-1] eta)
*
* where eta is a pseudo-fermion field, generated randomly with distribution
* proportional to exp{-(eta,eta)}.
*
* The arguments of the program rwrtm() are:
*
*  irp1,irp2     Indices of the Zolotarev rational functions R and R'
*                respectively, as stored in the parameter data base.
*                It is checked that R and R' have the same order.
*
*  np1,np2       Integers defining the pole range to be used in P(np1,np2).
*                It is checked that 0<=np1<=np2<n, where n is the order of
*                both rational approxximations.
*
*  isp           Indices of the solvers to be used for the computation of the
*                action of P(np1,np2) on a given spinor field (the supported
*                solvers are MSCG, SAP_GCR and DFL_SAP_GCR).
*
*  sqn           Square norm of the generated random pseudo-fermion field.
*
*  status        Array of the average of the status variables returned by the
*                solvers. The array status must have as many elements as are
*                returned by the solver with index isp (1 for MSCG and SAP_GCR
*                and 3 for DFL_SAP_GCR). In the case of the DFL_SAP_GCR solver,
*                status[2] reports the number of subspace regenerations that
*                were required.
*
* It is taken for granted that the solver parameters have been set by
* set_solver_parms() [flags/sparms.c]. The bare quark mass is taken to be
* the one last set by set_sw_parms() [flags/lat_parms.c].
*
* The program requires a workspace of double-precision spinor fields, whose
* size depends on the chosen solver. The required workspace is 3+(np2-np1+1)
* for MSCG, and 4 for SAP_GCR or DFL_SAP_GCR. Note that these figures do not
* include the workspace for the solver programs (see forces/tmcgm.c,
* sap/sap_gcr.c and dfl/dfl_sap_gcr.c).
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes. Some debugging information
* is printed to stdout if the macro RWRAT_DBG is defined.
*
*******************************************************************************/

#define RWRTM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "sflds.h"
#include "sw_term.h"
#include "dirac.h"
#include "linalg.h"
#include "sap.h"
#include "dfl.h"
#include "forces.h"
#include "update.h"
#include "global.h"


static int nps=0;
static double *rs;


#ifdef REWEIGHT_CHECK_PROGRAMS

extern void Dwhat_dble_rcp(double mu,spinor_dble *s,spinor_dble *r);
extern void tmcgm_rcp(int nmx,double *res,int nmu,double *mu,
                      spinor_dble *eta,spinor_dble **psi,int *status);
extern double tmcgeo_rcp(int nmx,double res,double mu,
                     spinor_dble *eta,spinor_dble *psi,int *status);
extern double sap_gcr_rcp(int nkv,int nmx,double res,double mu,
                          spinor_dble *eta,spinor_dble *psi,int *status);
extern double dfl_sap_gcr2_rcp(int idfl,int nkv,int nmx,double res,double mu,
                               spinor_dble *eta,spinor_dble *psi,int *status);

#define Dwhat_dble Dwhat_dble_rcp
#define tmcgm tmcgm_rcp
#define tmcgeo tmcgeo_rcp
#define sap_gcr sap_gcr_rcp
#define dfl_sap_gcr2 dfl_sap_gcr2_rcp

#endif


static double set_eta(spinor_dble *eta)
{
   random_sd(VOLUME/2,eta,1.0);
   set_sd2zero(VOLUME/2,eta+(VOLUME/2));
   bnd_sd2zero(EVEN_PTS,eta);

   return norm_square_dble(VOLUME/2,1,eta);
}


static void set_res(int np,double res)
{
   int k;

   if (np>nps)
   {
      if (nps>0)
         free(rs);

      rs=malloc(np*sizeof(*rs));
      error(rs==NULL,1,"set_res [rwrtm.c]",
            "Unable to allocate auxiliary array");
      nps=np;
   }

   for (k=0;k<np;k++)
      rs[k]=res;
}


/*******************************************************************************
psi += sum_k ( rnu1(k)/(Dwhat^dag*Dwhat+nu1(k)^2)
               + rmu2(k)/(Dwhat^dag*Dwhat+mu2(k)^2) ) eta
*******************************************************************************/
static void apply_Rk(int np,int isp,double *nu1,double *rnu1,double *mu2,
                     double *rmu2,spinor_dble *eta,spinor_dble *psi,int *status)
{
   int k,l,stat[6];
   spinor_dble **rsd;
   solver_parms_t sp;
   sap_parms_t sap;

   sp=solver_parms(isp);
   
   if (sp.solver==MSCG)
   {
      rsd=reserve_wsd(np+1);

      set_res(np,sp.res);
      tmcgm(sp.nmx,rs,np,nu1,eta,rsd,status);

      error_root(status[0]<0,1,"apply_Rk [rwrtm.c]","MSCG solver failed "
                 "(isp=%d, status=%d)",isp,status[0]);

      status[1]=0;
      for (k=0;k<np;k++)
      {
         tmcgeo(sp.nmx,sp.res,mu2[k],rsd[k],rsd[np],stat);

         error_root(stat[0]<0,1,"apply_Rk [rwrtm.c]","CGNE solver failed "
                    "(k=%d, isp=%d, stat=%d)",k,isp,stat[0]);
      
         status[1]+=stat[0];

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[np],rnu1[k]*mu2[k]*mu2[k]+rmu2[k]*nu1[k]*nu1[k]);

         sw_term(ODD_PTS);
         Dwhat_dble(0.0,rsd[np],rsd[k]);
         mulg5_dble(VOLUME/2,rsd[k]);
         Dwhat_dble(0.0,rsd[k],rsd[np]);
         mulg5_dble(VOLUME/2,rsd[np]);

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[np],rnu1[k]+rmu2[k]);
      }

      release_wsd();
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(2);
      mulg5_dble(VOLUME/2,eta);
      set_sd2zero(VOLUME/2,eta+(VOLUME/2));
      status[0]=0;

      for (k=0;k<np;k++)
      {
         sap_gcr(sp.nkv,sp.nmx,sp.res,nu1[k],eta,rsd[0],stat);
         mulg5_dble(VOLUME/2,rsd[0]);
         set_sd2zero(VOLUME/2,rsd[0]+(VOLUME/2));
         sap_gcr(sp.nkv,sp.nmx,sp.res,-nu1[k],rsd[0],rsd[1],stat+1);

         error_root((stat[0]<0)||(stat[1]<0),1,"apply_Rk [rwrtm.c]",
                    "SAP_GCR solver failed (isp=%d, status=%d;%d)",
                    isp,stat[0],stat[1]);

         status[0]+=stat[0];
         status[0]+=stat[1];

         mulg5_dble(VOLUME/2,rsd[1]);
         sap_gcr(sp.nkv,sp.nmx,sp.res,mu2[k],rsd[1],rsd[0],stat);
         mulg5_dble(VOLUME/2,rsd[0]);
         set_sd2zero(VOLUME/2,rsd[0]+(VOLUME/2));
         sap_gcr(sp.nkv,sp.nmx,sp.res,-mu2[k],rsd[0],rsd[1],stat+1);

         error_root((stat[0]<0)||(stat[1]<0),1,"apply_Rk [rwrtm.c]",
                    "SAP_GCR solver failed (isp=%d, status=%d;%d)",
                    isp,stat[0],stat[1]);

         status[0]+=stat[0];
         status[0]+=stat[1];

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[1],rnu1[k]*mu2[k]*mu2[k]+rmu2[k]*nu1[k]*nu1[k]);

         sw_term(ODD_PTS);
         Dwhat_dble(0.0,rsd[1],rsd[0]);
         mulg5_dble(VOLUME/2,rsd[0]);
         Dwhat_dble(0.0,rsd[0],rsd[1]);
         mulg5_dble(VOLUME/2,rsd[1]);

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[1],rnu1[k]+rmu2[k]);
      }

      status[0]=(status[0]+2*np)/(4*np);
      mulg5_dble(VOLUME/2,eta);
      release_wsd();
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      rsd=reserve_wsd(2);
      mulg5_dble(VOLUME/2,eta);
      set_sd2zero(VOLUME/2,eta+(VOLUME/2));

      for (l=0;l<3;l++)
         status[l]=0;

      for (k=0;k<np;k++)
      {
         dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,nu1[k],eta,rsd[0],stat);
         mulg5_dble(VOLUME/2,rsd[0]);
         set_sd2zero(VOLUME/2,rsd[0]+(VOLUME/2));
         dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,-nu1[k],rsd[0],rsd[1],stat+3);

         error_root((stat[0]<0)||(stat[1]<0)||(stat[3]<0)||(stat[4]<0),1,
                    "apply_Rk [rwrtm.c]","DFL_SAP_GCR solver failed (isp=%d, "
                    "status=%d,%d,%d;%d,%d,%d)",isp,stat[0],stat[1],
                    stat[2],stat[3],stat[4],stat[5]);

         for (l=0;l<3;l++)
         {
            status[l]+=stat[l];
            status[l]+=stat[l+3];
         }

         status[2]+=(stat[2]!=0);
         status[2]+=(stat[5]!=0);

         mulg5_dble(VOLUME/2,rsd[1]);
         dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,mu2[k],rsd[1],rsd[0],stat);
         mulg5_dble(VOLUME/2,rsd[0]);
         set_sd2zero(VOLUME/2,rsd[0]+(VOLUME/2));
         dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,-mu2[k],rsd[0],rsd[1],stat+3);

         error_root((stat[0]<0)||(stat[1]<0)||(stat[3]<0)||(stat[4]<0),1,
                    "apply_Rk [rwrtm.c]","DFL_SAP_GCR solver failed (isp=%d, "
                    "status=%d,%d,%d;%d,%d,%d)",isp,stat[0],stat[1],
                    stat[2],stat[3],stat[4],stat[5]);

         for (l=0;l<3;l++)
         {
            status[l]+=stat[l];
            status[l]+=stat[l+3];
         }

         status[2]+=(stat[2]!=0);
         status[2]+=(stat[5]!=0);


         mulr_spinor_add_dble(VOLUME/2,psi,rsd[1],rnu1[k]*mu2[k]*mu2[k]+rmu2[k]*nu1[k]*nu1[k]);

         sw_term(ODD_PTS);
         Dwhat_dble(0.0,rsd[1],rsd[0]);
         mulg5_dble(VOLUME/2,rsd[0]);
         Dwhat_dble(0.0,rsd[0],rsd[1]);
         mulg5_dble(VOLUME/2,rsd[1]);

         mulr_spinor_add_dble(VOLUME/2,psi,rsd[1],rnu1[k]+rmu2[k]);
      }

      for (l=0;l<2;l++)
         status[l]=(status[l]+2*np)/(4*np);

      mulg5_dble(VOLUME/2,eta);
      release_wsd();
   }
   else
      error_root(1,1,"apply_Rk [rwrtm.c]","Unknown or unsupported solver");
}


double rwrtm(int irp1,int irp2,int np1,int np2,int isp,double *sqn,int *status)
{
   int j,k,n;
   double lnr,*nu1,*nu2,*mu1,*mu2,*rnu1,*rmu2,z1,z2;
   spinor_dble *eta,*psi,**wsd;
   tm_parms_t tm;
   rat_parms_t rp1,rp2;
   
   check_global_int("rwrtm_single [rwrtm.c]",5,irp1,irp2,np1,np2,isp);

   tm=tm_parms();
   if (tm.eoflg!=1)
      set_tm_parms(1);

   rp1=rat_parms(irp1);
   rp2=rat_parms(irp2);

   error_root((np1<0)||(np2<np1)||(np2>=rp1.degree),1,
              "rwrtm_single [rwrtm.c]",
              "Parameters np1 and np2 are out of range");


   wsd=reserve_wsd(2);
   eta=wsd[0];
   psi=wsd[1];
   (*sqn)=set_eta(eta);
   
   nu1=rp1.nu+np1;
   mu1=rp1.mu+np1;
   nu2=rp2.nu+np1;
   mu2=rp2.mu+np1;
   
   n=np2-np1+1;
   rnu1=malloc(2*n*sizeof(double));
   rmu2=rnu1+n;
   
   for (j=0;j<n;j++)
   {
      z1=nu1[j]*nu1[j];
      z2=mu2[j]*mu2[j];
      rnu1[j]=(mu1[j]*mu1[j]-z1)*(nu2[j]*nu2[j]-z1)/(mu2[j]*mu2[j]-z1);
      rmu2[j]=(mu1[j]*mu1[j]-z2)*(nu2[j]*nu2[j]-z2)/(nu1[j]*nu1[j]-z2);
      for (k=0;k<n;k++)
      {
         if (k!=j)
         {
            rnu1[j]*=(mu1[k]*mu1[k]-z1)*(nu2[k]*nu2[k]-z1)/
                      ((nu1[k]*nu1[k]-z1)*(mu2[k]*mu2[k]-z1));
            rmu2[j]*=(mu1[k]*mu1[k]-z2)*(nu2[k]*nu2[k]-z2)/
                      ((nu1[k]*nu1[k]-z2)*(mu2[k]*mu2[k]-z2));
         }
      }
   }

#ifdef RWRAT_DBG
   message("[rwrtm]: irp = (%d,%d), poles = %d ... %d\n",irp1,irp2,np1,np2);
   message("[rwrtm]: sqn = %.4e\n\n",*sqn);
   message("[rwrtm]: denominators           numerators             residues\n");
   for (j=0;j<n;j++)
   {
      message("[rwrtm]: nu1[%2d] = %.4e   mu1[%2d] = %.4e   rnu1[%2d] = %+.4e\n",np1+j,nu1[j],np1+j,mu1[j],np1+j,rnu1[j]);
      message("[rwrtm]: mu2[%2d] = %.4e   nu2[%2d] = %.4e   rmu2[%2d] = %+.4e\n",np1+j,mu2[j],np1+j,nu2[j],np1+j,rmu2[j]);
   }
#endif
   
   set_sd2zero(VOLUME,psi);
   apply_Rk(n,isp,nu1,rnu1,mu2,rmu2,eta,psi,status);
   
   lnr=spinor_prod_re_dble(VOLUME/2,1,eta,psi);

#ifdef RWRAT_DBG
   message("\n[rwrtm]: -ln(r) = %.6e\n\n",lnr);
#endif

   release_wsd();
   free(rnu1);

   return lnr;
}
