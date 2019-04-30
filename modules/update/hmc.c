
/*******************************************************************************
*
* File hmc.c
*
* Copyright (C) 2005, 2007, 2009-2013 Martin Luescher, Filippo Palombi,
*               2016                  Stefan Schaefer, Isabel Campos
*               2017                  Agostino Patella
*               2019                  Agostino Patella, Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* HMC simulation algorithm.
*
* The externally accessible functions are
*
*   void hmc_sanity_check(void)
*     Performs various checks on the chosen parameters for the HMC
*     algorithm and terminates with an error message if an inconsistency
*     is discovered.
*
*   void hmc_wsize(int *nwud,int *nwad,int *nws,int *nwsd,int *nwv,int *nwvd)
*     Determines the minimal sizes of the workspaces required for the
*     HMC algorithm based on the information in the parameter data base.
*     On exit the program returns the numbers of double-precision SU(3) gauge
*     (nwud), double-precision U(1) gauge (nwad), spinor (nwsd) and complex
*     vector (nwvd) fields as well as the numbers of single-precision spinor
*     (nws) and complex vector (nwv) fields that must be allocated.
*
*   int run_hmc(double *act0,double *act1)
*     Generates a random momentum field, integrates the MD equations and
*     applies the HMC acceptance step to the fields at the end of the MD
*     trajectory (see the notes).
*      The arrays act0 and act1 must have at least nact+1 elements, where
*     nact is the number of actions that take part in the HMC algorithm
*     (see flags/hmc_parms.c). On exit act0 and act1 contain the part of
*     the actions computed on the local lattice at the beginning and the
*     end of the MD evolution (see the notes).
*      The program returns 1 or 0 depending on whether the field generated
*     by the molecular-dynamics evolution was accepted or not. If it was
*     not accepted, the gauge field is restored to its initial value.
*
* Notes:
*
* The molecular-dynamics equations are integrated using the integrator
* specified by the list of elementary operations returned by mdsteps()
* (see update/mdsteps.c and update/mdint.c). The elements of the action
* arrays act0 and act1 are
*
*  actx[0]        Action of the momentum field,
*  actx[1+n]      Action number n,
*
* where the pseudo-fermion actions are counted from 0 in steps of 1, as
* they appear in the action array hmc.iact returned by hmc_parms().
*
* The boundary conditions are imposed as specified in the parameter data
* base (see flags/lat_parms.c). Accepted new gauge field configurations
* are renormalized to SU(3) on all active links.
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define HMC_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "lattice.h"
#include "utils.h"
#include "uflds.h"
#include "u1flds.h"
#include "mdflds.h"
#include "linalg.h"
#include "dfl.h"
#include "forces.h"
#include "update.h"
#include "global.h"

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static int nrs=0,*rs;
static su3_dble ubnd[3] ALIGNED16;
static double abnd[3] ALIGNED16;


static void init_rs(int nr)
{
   int k;

   if (nr>nrs)
   {
      if (nrs>0)
         free(rs);

      rs=malloc(nr*sizeof(*rs));
      error_root(rs==NULL,1,"init_rs [hmc.c]",
                 "Unable to allocate auxiliary array");
      nrs=nr;
   }

   for (k=0;k<nr;k++)
      rs[k]=0;
}


static int check_rat_actions(void)
{
   int k,l,j,ie;
   int nact,*iact,ir,nr,ifl,isw;
   hmc_parms_t hmc;
   action_parms_t ap;
   rat_parms_t rp;

   hmc=hmc_parms();
   nact=hmc.nact;
   iact=hmc.iact;
   ie=0;

   for (k=0;k<nact;k++)
   {
      ap=action_parms(iact[k]);

      if ((ap.action==ACF_RAT)||(ap.action==ACF_RAT_SDET))
      {
         ir=ap.irat[0];
         ifl=ap.ifl;
         rp=rat_parms(ir);
         nr=rp.degree;
         init_rs(nr);
         isw=0;

         for (l=0;l<nact;l++)
         {
            ap=action_parms(iact[l]);

            if ((ap.action==ACF_RAT)||(ap.action==ACF_RAT_SDET))
            {
               if ((ap.irat[0]==ir)&&(ap.ifl==ifl))
               {
                  if (ap.action==ACF_RAT_SDET)
                     isw+=1;

                  for (j=ap.irat[1];j<=ap.irat[2];j++)
                     rs[j]+=1;
               }
            }
         }

         for (l=0;l<nr;l++)
            ie|=(rs[l]!=isw);
      }
   }

   return ie;
}


static int match_force(action_parms_t ap,force_parms_t fp)
{
   int ie;

   ie=1;

   if (ap.action==ACG_SU3)
      ie&=(fp.force==FRG_SU3);
   else if (ap.action==ACG_U1)
      ie&=(fp.force==FRG_U1);
   else if (ap.action==ACF_TM1)
   {
      ie&=(fp.force==FRF_TM1);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.imu[0]==fp.imu[0]);
   }
   else if (ap.action==ACF_TM1_EO)
   {
      ie&=(fp.force==FRF_TM1_EO);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.imu[0]==fp.imu[0]);
   }
   else if (ap.action==ACF_TM1_EO_SDET)
   {
      ie&=(fp.force==FRF_TM1_EO_SDET);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.imu[0]==fp.imu[0]);
   }
   else if (ap.action==ACF_TM2)
   {
      ie&=(fp.force==FRF_TM2);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.imu[0]==fp.imu[0]);
      ie&=(ap.imu[1]==fp.imu[1]);
   }
   else if (ap.action==ACF_TM2_EO)
   {
      ie&=(fp.force==FRF_TM2_EO);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.imu[0]==fp.imu[0]);
      ie&=(ap.imu[1]==fp.imu[1]);
   }
   else if (ap.action==ACF_RAT)
   {
      ie&=(fp.force==FRF_RAT);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.irat[0]==fp.irat[0]);
      ie&=(ap.irat[1]==fp.irat[1]);
      ie&=(ap.irat[2]==fp.irat[2]);
   }
   else if (ap.action==ACF_RAT_SDET)
   {
      ie&=(fp.force==FRF_RAT_SDET);
      ie&=(ap.ipf==fp.ipf);
      ie&=(ap.ifl==fp.ifl);
      ie&=(ap.irat[0]==fp.irat[0]);
      ie&=(ap.irat[1]==fp.irat[1]);
      ie&=(ap.irat[2]==fp.irat[2]);
   }
   else
      ie=0;

   return ie;
}


void hmc_sanity_check(void)
{
   int my_rank;
   int nlv,nact,*iact,npf,nmu;
   int iepf,iemu,iem0,ierat,iacg3,iacg1;
   int ie,ic;
   int nfr,*ifr,i,j,k;
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   action_parms_t ap;
   force_parms_t fp;
   rat_parms_t rp;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      hmc=hmc_parms();
      nlv=hmc.nlv;
      nact=hmc.nact;
      iact=hmc.iact;
      npf=hmc.npf;
      nmu=hmc.nmu;

      error_root((nlv<1)||(nact<1),1,"hmc_sanity_check [hmc.c]",
                 "hmc.nlv or hmc.nact must be at least 1");

      iepf=0;
      iemu=0;
      iem0=0;
      ierat=0;
      iacg1=0;
      iacg3=0;

      for (i=0;i<nact;i++)
      {
         ap=action_parms(iact[i]);

         if (ap.action==ACG_SU3)
            iacg3+=1;
         else if (ap.action==ACG_U1)
            iacg1+=1;
         else if ((ap.action==ACF_TM1)||
                  (ap.action==ACF_TM1_EO)||
                  (ap.action==ACF_TM1_EO_SDET)||
                  (ap.action==ACF_TM2)||
                  (ap.action==ACF_TM2_EO))
         {
            iepf|=(ap.ipf<0);
            iepf|=(ap.ipf>=npf);
            iemu|=(ap.imu[0]<0);
            iemu|=(ap.imu[0]>=nmu);
            iem0|=(qlat_parms(ap.ifl).m0==DBL_MAX);

            if ((ap.action==ACF_TM2)||
                (ap.action==ACF_TM2_EO))
            {
               iemu|=(ap.imu[1]<0);
               iemu|=(ap.imu[1]>=nmu);
            }
         }
         else if ((ap.action==ACF_RAT)||
                  (ap.action==ACF_RAT_SDET))
         {
            iepf|=(ap.ipf<0);
            iepf|=(ap.ipf>=npf);
            iem0|=(qlat_parms(ap.ifl).m0==DBL_MAX);

            rp=rat_parms(ap.irat[0]);
            ierat|=(ap.irat[2]>=rp.degree);
         }
      }

      for (i=0;i<nlv;i++)
      {
         mdp=mdint_parms(i);
         nfr=mdp.nfr;
         ifr=mdp.ifr;

         for (j=0;j<nfr;j++)
         {
            fp=force_parms(ifr[j]);

            if ((fp.force==FRF_TM1)||
                (fp.force==FRF_TM1_EO)||
                (fp.force==FRF_TM1_EO_SDET)||
                (fp.force==FRF_TM2)||
                (fp.force==FRF_TM2_EO))
            {
               iepf|=(fp.ipf<0);
               iepf|=(fp.ipf>=npf);
               iemu|=(fp.imu[0]<0);
               iemu|=(fp.imu[0]>=nmu);
               iem0|=(qlat_parms(fp.ifl).m0==DBL_MAX);

               if ((fp.force==FRF_TM2)||
                   (fp.force==FRF_TM2_EO))
               {
                  iemu|=(fp.imu[1]<0);
                  iemu|=(fp.imu[1]>=nmu);
               }
            }
            else if ((fp.force==FRF_RAT)||
                     (fp.force==FRF_RAT_SDET))
            {
               iepf|=(fp.ipf<0);
               iepf|=(fp.ipf>=npf);
               iem0|=(qlat_parms(fp.ifl).m0==DBL_MAX);

               rp=rat_parms(fp.irat[0]);
               ierat|=(fp.irat[2]>=rp.degree);
            }
         }
      }

      error_root(iepf!=0,1,"hmc_sanity_check [hmc.c]",
                 "Some pseudo-fermion indices are out of range");
      error_root(iemu!=0,1,"hmc_sanity_check [hmc.c]",
                 "Some twisted-mass indices are out of range");
      error_root(iem0!=0,1,"hmc_sanity_check [hmc.c]",
                 "Some sea-quark mass indices are out of range");
      error_root(ierat!=0,1,"hmc_sanity_check [hmc.c]",
                 "Some rational functions are not correctly specified");
      error_root(((gauge()&1)!=0)&&iacg3==0,1,"hmc_sanity_check [hmc.c]",
                 "SU(3) gauge action is missing");
      error_root(((gauge()&1)==0)&&iacg3>0,1,"hmc_sanity_check [hmc.c]",
                 "SU(3) gauge action is requested but SU(3) gauge field is inactive");
      error_root(iacg3>1,1,"hmc_sanity_check [hmc.c]",
                 "SU(3) gauge action occurs several times");
      error_root(((gauge()&2)!=0)&&iacg1==0,1,"hmc_sanity_check [hmc.c]",
                 "U(1) gauge action is missing");
      error_root(((gauge()&2)==0)&&iacg1>0,1,"hmc_sanity_check [hmc.c]",
                 "U(1) gauge action is requested but U(1) gauge field is inactive");
      error_root(iacg1>1,1,"hmc_sanity_check [hmc.c]",
                 "U(1) gauge action occurs several times");

      ie=0;

      for (k=0;k<nact;k++)
      {
         ap=action_parms(iact[k]);
         ic=0;

         for (i=0;i<nlv;i++)
         {
            mdp=mdint_parms(i);
            nfr=mdp.nfr;
            ifr=mdp.ifr;

            for (j=0;j<nfr;j++)
            {
               fp=force_parms(ifr[j]);
               ic+=match_force(ap,fp);
            }
         }

         ie|=(ic!=1);
      }

      for (i=0;i<nlv;i++)
      {
         mdp=mdint_parms(i);
         nfr=mdp.nfr;
         ifr=mdp.ifr;

         for (j=0;j<nfr;j++)
         {
            fp=force_parms(ifr[j]);
            ic=0;

            for (k=0;k<nact;k++)
            {
               ap=action_parms(iact[k]);
               ic+=match_force(ap,fp);
            }

            ie|=(ic!=1);
         }
      }

      error_root(ie!=0,1,"hmc_sanity_check [hmc.c]",
                 "Specified actions and forces do not match");

      ie=check_rat_actions();
      error_root(ie!=0,1,"hmc_sanity_check [hmc.c]",
                 "Inconsistent rational function actions");
   }
}


static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpp;

   dp=dfl_parms();
   dpp=dfl_pro_parms();

   MAX(*nws,dp.Ns+2);
   MAX(*nwv,2*dpp.nkv+2);
   MAX(*nwvd,4);
}


static void solver_wsize(int isp,int nsd,int np,
                         int *nws,int *nwsd,int *nwv,int *nwvd)
{
   solver_parms_t sp;

   sp=solver_parms(isp);

   if (sp.solver==CGNE)
   {
      MAX(*nws,5);
      MAX(*nwsd,nsd+5);
   }
   else if (sp.solver==MSCG)
   {
      if (np>1)
      {
         MAX(*nwsd,nsd+np+3);
      }
      else
      {
         MAX(*nwsd,nsd+5);
      }
   }
   else if (sp.solver==SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+1);
      MAX(*nwsd,nsd+2);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+2);
      MAX(*nwsd,nsd+3);
      dfl_wsize(nws,nwv,nwvd);
   }
}


void hmc_wsize(int *nwud,int *nwad,int *nws,int *nwsd,int *nwv,int *nwvd)
{
   int nlv,nact,*iact;
   int nfr,*ifr,nsd,np,i,j;
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   action_parms_t ap;
   force_parms_t fp;
   solver_parms_t sp;

   (*nwud)=(gauge()&1);
   (*nwad)=(gauge()&2)/2;
   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;

   hmc=hmc_parms();
   nlv=hmc.nlv;
   nact=hmc.nact;
   iact=hmc.iact;

   for (i=0;i<nact;i++)
   {
      ap=action_parms(iact[i]);

      if ((ap.action==ACF_TM1)||
          (ap.action==ACF_TM1_EO)||
          (ap.action==ACF_TM1_EO_SDET)||
          (ap.action==ACF_TM2)||
          (ap.action==ACF_TM2_EO))
      {
         nsd=1;
         solver_wsize(ap.isp[0],nsd,0,nws,nwsd,nwv,nwvd);

         if ((ap.action==ACF_TM2)||
             (ap.action==ACF_TM2_EO))
            solver_wsize(ap.isp[1],nsd,0,nws,nwsd,nwv,nwvd);
      }
      else if ((ap.action==ACF_RAT)||
               (ap.action==ACF_RAT_SDET))
      {
         np=ap.irat[2]-ap.irat[1]+1;
         sp=solver_parms(ap.isp[0]);

         if (sp.solver==MSCG)
            nsd=np;
         else
            nsd=2;

         solver_wsize(ap.isp[0],nsd,np,nws,nwsd,nwv,nwvd);
      }
   }

   for (i=0;i<nlv;i++)
   {
      mdp=mdint_parms(i);
      nfr=mdp.nfr;
      ifr=mdp.ifr;

      for (j=0;j<nfr;j++)
      {
         fp=force_parms(ifr[j]);

         if ((fp.force==FRF_TM1)||
             (fp.force==FRF_TM1_EO)||
             (fp.force==FRF_TM1_EO_SDET)||
             (fp.force==FRF_TM2)||
             (fp.force==FRF_TM2_EO))
         {
            sp=solver_parms(fp.isp[0]);

            if (fp.icr[0]==0)
               nsd=2;
            else if (sp.solver==CGNE)
               nsd=3;
            else
               nsd=4;

            solver_wsize(fp.isp[0],nsd,0,nws,nwsd,nwv,nwvd);
         }
         else if ((fp.force==FRF_RAT)||
                  (fp.force==FRF_RAT_SDET))
         {
            np=fp.irat[2]-fp.irat[1]+1;
            sp=solver_parms(fp.isp[0]);

            if (sp.solver==MSCG)
               nsd=np;
            else
               nsd=3;

            solver_wsize(fp.isp[0],nsd,np,nws,nwsd,nwv,nwvd);
         }
      }
   }
}


static void chk_mode_regen(int isp,int *status)
{
   solver_parms_t sp;

   sp=solver_parms(isp);

   if ((sp.solver==DFL_SAP_GCR)&&(status[2]>0))
      add2counter("modes",2+3*sp.idfl,status+2);
}


static void start_hmc(double *act0,su3_dble *uold,double *aold)
{
   int i,nact,*iact,idfl;
   int status[3];
   double *mu;
   su3_dble *udb;
   double *adb;
   dfl_parms_t dfl;
   hmc_parms_t hmc;
   action_parms_t ap;
   dirac_parms_t dp;
   dflst_t dfl_status;

   hmc=hmc_parms();
   nact=hmc.nact;
   iact=hmc.iact;
   mu=hmc.mu;

   clear_counters();

   act0[0]=0.0;
   if((gauge()&1)!=0)
   {
      udb=udfld();
      cm3x3_assign(4*VOLUME,udb,uold);
      if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
         cm3x3_assign(3,udb+4*VOLUME+7*(BNDRY/4),ubnd);
      random_su3mom();
      act0[0]+=su3momentum_action(0);
   }
   if((gauge()&2)!=0)
   {
      adb=adfld();
      assign_dvec2dvec(4*VOLUME,adb,aold);
      if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
         assign_dvec2dvec(3,adb+4*VOLUME+7*(BNDRY/4),abnd);
      random_u1mom();
      act0[0]+=u1momentum_action(0);
      if (hmc.facc)
         u1mom_Delta_no0(2,mdflds()->u1mom,mdflds()->u1mom);
   }

   dfl=dfl_parms();

   if (dfl.Ns)
   {
      idfl=0;
      while(1)
      {
         dfl_status=dfl_gen_parms(idfl).status;
         if(dfl_status==DFL_OUTOFRANGE) break;
         if(dfl_status==DFL_DEF)
         {
            dfl_modes2(idfl,status);
            error_root((status[1]<0)||((status[1]==0)&&(status[0]<0)),1,
                       "start_hmc [hmc.c]","Generation of deflation subspace %d "
                       "failed (status = %d;%d)",idfl,status[0],status[1]);

            if (status[1]==0)
               add2counter("modes",3*idfl,status);
            else
               add2counter("modes",2+3*idfl,status+1);
         }
         idfl++;
      }
   }

   for (i=0;i<nact;i++)
   {
      ap=action_parms(iact[i]);

      if (ap.action==ACG_SU3)
         act0[i+1]=action0(0);
      else if (ap.action==ACG_U1)
         act0[i+1]=action6(0);
      else
      {
         dp=qlat_parms(ap.ifl);
         set_dirac_parms1(&dp);

         if (ap.action==ACF_TM1)
            act0[i+1]=setpf1(mu[ap.imu[0]],ap.ipf,0);
         else if (ap.action==ACF_TM1_EO)
            act0[i+1]=setpf4(mu[ap.imu[0]],ap.ipf,0,0);
         else if (ap.action==ACF_TM1_EO_SDET)
            act0[i+1]=setpf4(mu[ap.imu[0]],ap.ipf,1,0);
         else if (ap.action==ACF_TM2)
         {
            status[2]=0;
            act0[i+1]=setpf2(mu[ap.imu[0]],mu[ap.imu[1]],ap.ipf,ap.isp[1],
                           0,status);
            chk_mode_regen(ap.isp[1],status);
            add2counter("field",ap.ipf,status);
         }
         else if (ap.action==ACF_TM2_EO)
         {
            status[2]=0;
            act0[i+1]=setpf5(mu[ap.imu[0]],mu[ap.imu[1]],ap.ipf,ap.isp[1],
                           0,status);
            chk_mode_regen(ap.isp[1],status);
            add2counter("field",ap.ipf,status);
         }
         else if (ap.action==ACF_RAT)
         {
            status[2]=0;
            act0[i+1]=setpf3(ap.irat,ap.ipf,0,ap.isp[0],0,status);
            chk_mode_regen(ap.isp[0],status);
            add2counter("field",ap.ipf,status);
         }
         else if (ap.action==ACF_RAT_SDET)
         {
            status[2]=0;
            act0[i+1]=setpf3(ap.irat,ap.ipf,1,ap.isp[0],0,status);
            chk_mode_regen(ap.isp[0],status);
            add2counter("field",ap.ipf,status);
         }
         else
            error_root(1,1,"start_hmc [hmc.c]","Unknown action");
      }
   }
}


static void end_hmc(double *act1)
{
   int i,nact,*iact;
   int status[3];
   double *mu;
   hmc_parms_t hmc;
   action_parms_t ap;
   dirac_parms_t dp;

   hmc=hmc_parms();

   act1[0]=0.0;
   if((gauge()&1)!=0)
      act1[0]+=su3momentum_action(0);
   if((gauge()&2)!=0)
   {
      if (hmc.facc)
         u1mom_Delta_no0(3,mdflds()->u1mom,mdflds()->u1mom);
      act1[0]+=u1momentum_action(0);
   }
   nact=hmc.nact;
   iact=hmc.iact;
   mu=hmc.mu;

   for (i=0;i<nact;i++)
   {
      ap=action_parms(iact[i]);

      if (ap.action==ACG_SU3)
         act1[i+1]=action0(0);
      else if (ap.action==ACG_U1)
         act1[i+1]=action6(0);
      else
      {
         dp=qlat_parms(ap.ifl);
         set_dirac_parms1(&dp);
         status[2]=0;

         if (ap.action==ACF_TM1)
            act1[i+1]=action1(mu[ap.imu[0]],ap.ipf,ap.isp[0],0,status);
         else if (ap.action==ACF_TM1_EO)
            act1[i+1]=action4(mu[ap.imu[0]],ap.ipf,0,ap.isp[0],0,status);
         else if (ap.action==ACF_TM1_EO_SDET)
            act1[i+1]=action4(mu[ap.imu[0]],ap.ipf,1,ap.isp[0],0,status);
         else if (ap.action==ACF_TM2)
            act1[i+1]=action2(mu[ap.imu[0]],mu[ap.imu[1]],ap.ipf,ap.isp[0],
                            0,status);
         else if (ap.action==ACF_TM2_EO)
            act1[i+1]=action5(mu[ap.imu[0]],mu[ap.imu[1]],ap.ipf,ap.isp[0],
                            0,status);
         else if (ap.action==ACF_RAT)
            act1[i+1]=action3(ap.irat,ap.ipf,0,ap.isp[0],0,status);
         else if (ap.action==ACF_RAT_SDET)
            act1[i+1]=action3(ap.irat,ap.ipf,1,ap.isp[0],0,status);

         chk_mode_regen(ap.isp[0],status);
         add2counter("action",iact[i],status);
      }
   }
}


static int accept_hmc(double *act0,double *act1,su3_dble *uold,double *aold)
{
   int my_rank,nact,iac,i;
   double da,r;
   su3_dble *udb;
   double *adb;
   hmc_parms_t hmc;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   hmc=hmc_parms();
   nact=hmc.nact;
   iac=0;
   da=0.0;

   for (i=0;i<=nact;i++)
      da+=(act1[i]-act0[i]);

   if (NPROC>1)
   {
      r=da;
      MPI_Reduce(&r,&da,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   }

   if (my_rank==0)
   {
      ranlxd(&r,1);

      if (da<=0.0)
         iac=1;
      else if (r<=exp(-da))
         iac=1;
   }

   if (NPROC>1)
      MPI_Bcast(&iac,1,MPI_INT,0,MPI_COMM_WORLD);

   if (iac==0)
   {
      if ((gauge()&1)!=0)
      {
         udb=udfld();
         cm3x3_assign(4*VOLUME,uold,udb);
         if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
            cm3x3_assign(3,ubnd,udb+4*VOLUME+7*(BNDRY/4));
         set_flags(UPDATED_UD);
      }
      if ((gauge()&2)!=0)
      {
         adb=adfld();
         assign_dvec2dvec(4*VOLUME,aold,adb);
         if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
            assign_dvec2dvec(3,abnd,adb+4*VOLUME+7*(BNDRY/4));
         set_flags(UPDATED_AD);
      }
   }
   else
   {
      if ((gauge()&1)!=0)
         renormalize_ud();
      if ((gauge()&2)!=0)
         renormalize_ad();
   }

   return iac;
}


int run_hmc(double *act0,double *act1)
{
   int iac;
   su3_dble **usv,*uold=NULL;
   double **asv,*aold=NULL;

   if ((gauge()&1)!=0)
   {
      usv=reserve_wud(1);
      uold=usv[0];
   }
   if ((gauge()&2)!=0)
   {
      asv=reserve_wad(1);
      aold=asv[0];
   }

   start_hmc(act0,uold,aold);
   run_mdint();
   end_hmc(act1);
   iac=accept_hmc(act0,act1,uold,aold);

   if ((gauge()&1)!=0) release_wud();
   if ((gauge()&2)!=0) release_wad();

   return iac;
}
