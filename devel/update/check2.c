
/*******************************************************************************
*
* File check2.c
*gauge
* Copyright (C) 2005, 2007, 2009-2013, Martin Luescher, Filippo Palombi,
*               2016                   Stefan Schaefer,
*               2017                   Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reversibility of the MD evolution and stability of boundary conditions.
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
#include "linalg.h"
#include "uflds.h"
#include "u1flds.h"
#include "mdflds.h"
#include "archive.h"
#include "forces.h"
#include "dfl.h"
#include "update.h"
#include "global.h"

static int my_rank;


static void read_flds_bc_lat_parms(void)
{
   int gg,nfl,ifl,bc,sf,cs,type,qhat;
   double phi[2],phi_prime[2];
   double beta,c0,cG,cG_prime;
   double alpha,lambda,invqel;
   double kappa,su3csw,u1csw,cF,cF_prime,th1,th2,th3;
   char line[NAME_SIZE];

   if (my_rank==0)
   {
      find_section("Gauge group");
      read_line("gauge","%s",line);

      gg=0;
      if (strcmp(line,"SU(3)")==0)
         gg=1;
      else if (strcmp(line,"U(1)")==0)
         gg=2;
      else if (strcmp(line,"SU(3)xU(1)")==0)
         gg=3;
      else
         error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                    "Unknown gauge group %s",line);

      find_section("Quark action");
      read_line("nfl","%d",&nfl);
   }
   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nfl,1,MPI_INT,0,MPI_COMM_WORLD);

   set_flds_parms(gg,nfl);

   if (my_rank==0)
   {
      find_section("Boundary conditions");
      read_line("type","%s",&line);
      bc=4;
      if ((strcmp(line,"open")==0)||(strcmp(line,"0")==0))
         bc=0;
      else if ((strcmp(line,"SF")==0)||(strcmp(line,"1")==0))
         bc=1;
      else if ((strcmp(line,"open-SF")==0)||(strcmp(line,"2")==0))
         bc=2;
      else if ((strcmp(line,"periodic")==0)||(strcmp(line,"3")==0))
         bc=3;
      else
         error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                    "Unknown time boundary condition type %s",line);
      
      read_line("cstar","%d",&cs);

      sf=0;
      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
      if ((cs==0)&&(bc==1)&&((gg&1)!=0))
         read_dprms("phi",2,phi);
      if ((bc==1)||(bc==2))
      {
         read_line("SFtype","%s",&line);
         if ((strcmp(line,"orbifold")==0)||
             (strcmp(line,"openQCD-1.4")==0)||(strcmp(line,"0")==0))
            sf=0;
         else if ((strcmp(line,"AFW-typeB")==0)||
                  (strcmp(line,"openQCD-1.2")==0)||(strcmp(line,"1")==0))
            sf=1;
         else
            error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                       "Unknown SF type %s",line);
         
         if ((cs==0)&&((gg&1)!=0))
            read_dprms("phi'",2,phi_prime);
      }
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&sf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_bc_parms(bc,sf,cs,phi,phi_prime);
   
   if ((gg&1)!=0)
   {
      if (my_rank==0)
      {
         find_section("SU(3) action");
         read_line("beta","%lf",&beta);
         read_line("c0","%lf",&c0);

         cG=1.0;
         cG_prime=1.0;
         if (bc!=3)
            read_line("cG","%lf",&cG);
         if (bc==2)
            read_line("cG'","%lf",&cG_prime);
      }
      
      MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
      set_su3lat_parms(beta,c0,cG,cG_prime);
   }
   
   if ((gg&2)!=0)
   {
      if (my_rank==0)
      {
         find_section("U(1) action");
         read_line("type","%s",line);

         type=0;
         if ((strcmp(line,"compact")==0)||(strcmp(line,"0")==0))
            type=0;
         else if ((strcmp(line,"non-compact")==0)||(strcmp(line,"1")==0))
            type=1;
         else
            error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                       "Unknown U(1) action type %s",line);
      
         read_line("alpha","%lf",&alpha);
         read_line("invqel","%lf",&invqel);

         lambda=c0=cG=cG_prime=0.0;
         if(type==0)
         {
            read_line("c0","%lf",&c0);
            if (bc!=3) read_line("cG","%lf",&cG);
            if (bc==2) read_line("cG'","%lf",&cG_prime);
         }
         else
         {
            read_line("lambda","%lf",&lambda);
         }
      }
      
      MPI_Bcast(&type,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&invqel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&lambda,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
      set_u1lat_parms(type,alpha,invqel,lambda,c0,cG,cG_prime);
   }
   
   for (ifl=0;ifl<nfl;ifl++)
   {
      qhat=0;
      su3csw=u1csw=0.0;
      cF=cF_prime=0.0;
      th1=th2=th3=0.0;
      
      sprintf(line,"Flavour %d",ifl);
      if (my_rank==0)
      {
         find_section(line);
         read_line("kappa","%lf",&kappa);
         if (gg==1)
         {
            read_line("csw","%lf",&su3csw);
         }
         else if (gg==2)
         {
            read_line("qhat","%d",&qhat);
            read_line("csw","%lf",&u1csw);
         }
         else if (gg==3)
         {
            read_line("qhat","%d",&qhat);         
            read_line("su3csw","%lf",&su3csw);
            read_line("u1csw","%lf",&u1csw);
         }
         if (bc!=3) read_line("cF","%lf",&cF);
         if (bc==2) read_line("cF'","%lf",&cF_prime);
         if (cs==0) read_line("theta","%lf %lf %lf",&th1,&th2,&th3);
      }
      
      MPI_Bcast(&qhat,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&su3csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&u1csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&th1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&th2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&th3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      set_qlat_parms(ifl,kappa,qhat,su3csw,u1csw,cF,cF_prime,th1,th2,th3);
   }
}


static void read_hmc_parms(void)
{
   int nact,*iact;
   int npf,nmu,nlv;
   double tau,*mu;

   if (my_rank==0)
   {
      find_section("HMC parameters");
      nact=count_tokens("actions");
      read_line("npf","%d",&npf);
      nmu=count_tokens("mu");
      read_line("nlv","%d",&nlv);
      read_line("tau","%lf",&tau);
   }

   MPI_Bcast(&nact,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&npf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmu,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nlv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&tau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   if (nact>0)
   {
      iact=malloc(nact*sizeof(*iact));
      error(iact==NULL,1,"read_hmc_parms [check2.c]",
            "Unable to allocate temporary array");
      if (my_rank==0)
         read_iprms("actions",nact,iact);
      MPI_Bcast(iact,nact,MPI_INT,0,MPI_COMM_WORLD);
   }
   else
      iact=NULL;

   if (nmu>0)
   {
      mu=malloc(nmu*sizeof(*mu));
      error(mu==NULL,1,"read_hmc_parms [check2.c]",
            "Unable to allocate temporary array");
      if (my_rank==0)
         read_dprms("mu",nmu,mu);
      MPI_Bcast(mu,nmu,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
      mu=NULL;

   set_hmc_parms(nact,iact,npf,nmu,mu,nlv,tau);

   if (nact>0)
      free(iact);
   if (nmu>0)
      free(mu);
}


static void read_integrator(void)
{
   int nlv,i,j,k,l;
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   force_parms_t fp;
   rat_parms_t rp;

   hmc=hmc_parms();
   nlv=hmc.nlv;

   for (i=0;i<nlv;i++)
   {
      read_mdint_parms(i);
      mdp=mdint_parms(i);

      for (j=0;j<mdp.nfr;j++)
      {
         k=mdp.ifr[j];
         fp=force_parms(k);

         if (fp.force==FORCES)
            read_force_parms2(k);

         fp=force_parms(k);

         if ((fp.force==FRF_RAT)||(fp.force==FRF_RAT_SDET))
         {
            l=fp.irat[0];
            rp=rat_parms(l);

            if (rp.degree==0)
               read_rat_parms(l);
         }
      }
   }
}


static void read_actions(void)
{
   int i,k,l,nact,*iact;
   hmc_parms_t hmc;
   action_parms_t ap;
   rat_parms_t rp;

   hmc=hmc_parms();
   nact=hmc.nact;
   iact=hmc.iact;

   for (i=0;i<nact;i++)
   {
      k=iact[i];
      ap=action_parms(k);

      if (ap.action==ACTIONS)
         read_action_parms(k);

      ap=action_parms(k);

      if ((ap.action==ACF_RAT)||(ap.action==ACF_RAT_SDET))
      {
         l=ap.irat[0];
         rp=rat_parms(l);

         if (rp.degree==0)
            read_rat_parms(l);
      }
   }
}


static void read_sap_parms(void)
{
   int bs[4];

   if (my_rank==0)
   {
      find_section("SAP");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   set_sap_parms(bs,1,4,5);
}


static void read_dfl_parms(void)
{
   int bs[4],Ns;
   int ninv,nmr,ncy,nkv,nmx,nsm,qhat;
   double kappa,mu,su3csw,u1csw,cF,cF_prime,th1,th2,th3,res,dtau;

   if (my_rank==0)
   {
      find_section("Deflation subspace");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
      read_line("Ns","%d",&Ns);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_parms(bs,Ns);

   qhat=0;
   su3csw=u1csw=cF=cF_prime=th1=th2=th3=0.0;
   if (my_rank==0)
   {
      find_section("Deflation subspace generation");
      read_line("kappa","%lf",&kappa);
      read_line("mu","%lf",&mu);
      if((gauge()&2)!=0) read_line("qhat","%d",&qhat);
      if((gauge()&1)!=0) read_line("su3csw","%lf",&su3csw);
      if((gauge()&2)!=0) read_line("u1csw","%lf",&u1csw);
      if (bc_type()!=3) read_line("cF","%lf",&cF);
      if (bc_type()==2) read_line("cF'","%lf",&cF_prime);
      if(bc_cstar()==0) read_line("theta","%lf %lf %lf",&th1,&th2,&th3);
      read_line("ninv","%d",&ninv);
      read_line("nmr","%d",&nmr);
      read_line("ncy","%d",&ncy);
   }

   MPI_Bcast(&qhat,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&su3csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&th1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&th2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&th3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_gen_parms(kappa,mu,qhat,su3csw,u1csw,cF,cF_prime,th1,th2,th3,ninv,nmr,ncy);

   if (my_rank==0)
   {
      find_section("Deflation projection");
      read_line("nkv","%d",&nkv);
      read_line("nmx","%d",&nmx);
      read_line("res","%lf",&res);
   }

   MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   set_dfl_pro_parms(nkv,nmx,res);

   if (my_rank==0)
   {
      find_section("Deflation update scheme");
      read_line("dtau","%lf",&dtau);
      read_line("nsm","%d",&nsm);
   }

   MPI_Bcast(&dtau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nsm,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_upd_parms(dtau,nsm);
}


static void read_solvers(void)
{
   int nact,*iact,nlv,nsp;
   int nfr,*ifr;
   int isap,idfl,i,j,k;
   hmc_parms_t hmc;
   mdint_parms_t mdp;
   action_parms_t ap;
   force_parms_t fp;
   solver_parms_t sp;

   hmc=hmc_parms();
   nact=hmc.nact;
   iact=hmc.iact;
   nlv=hmc.nlv;
   isap=0;
   idfl=0;

   for (i=0;i<nact;i++)
   {
      ap=action_parms(iact[i]);

      if ((ap.action==ACF_TM1)||
          (ap.action==ACF_TM1_EO)||
          (ap.action==ACF_TM1_EO_SDET)||
          (ap.action==ACF_TM2)||
          (ap.action==ACF_TM2_EO)||
          (ap.action==ACF_RAT)||
          (ap.action==ACF_RAT_SDET))
      {
         if ((ap.action==ACF_TM2)||(ap.action==ACF_TM2_EO))
            nsp=2;
         else
            nsp=1;

         for (k=0;k<nsp;k++)
         {
            j=ap.isp[k];
            sp=solver_parms(j);

            if (sp.solver==SOLVERS)
            {
               read_solver_parms(j);
               sp=solver_parms(j);

               if (sp.solver==SAP_GCR)
                  isap=1;
               else if (sp.solver==DFL_SAP_GCR)
               {
                  isap=1;
                  idfl=1;
               }
            }
         }
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
             (fp.force==FRF_TM2_EO)||
             (fp.force==FRF_RAT)||
             (fp.force==FRF_RAT_SDET))
         {
            k=fp.isp[0];
            sp=solver_parms(k);

            if (sp.solver==SOLVERS)
            {
               read_solver_parms(k);
               sp=solver_parms(k);

               if (sp.solver==SAP_GCR)
                  isap=1;
               else if (sp.solver==DFL_SAP_GCR)
               {
                  isap=1;
                  idfl=1;
               }
            }
         }
      }
   }

   if (isap)
      read_sap_parms();

   if (idfl)
      read_dfl_parms();
}


static void chk_mode_regen(int isp,int *status)
{
   solver_parms_t sp;

   sp=solver_parms(isp);

   if ((sp.solver==DFL_SAP_GCR)&&(status[2]>0))
      add2counter("modes",2,status+2);
}


static void start_hmc(double *act0,su3_dble *uold,double *aold)
{
   int i,nact,*iact;
   int status[3];
   double *mu;
   su3_dble *udb;
   double *adb;
   dfl_parms_t dfl;
   hmc_parms_t hmc;
   action_parms_t ap;
   dirac_parms_t dp;

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
      random_su3mom();
      act0[0]+=su3momentum_action(0);
   }
   if((gauge()&2)!=0)
   {
      adb=adfld();
      assign_dvec2dvec(4*VOLUME,adb,aold);
      random_u1mom();
      act0[0]+=u1momentum_action(0);
   }

   dfl=dfl_parms();

   if (dfl.Ns)
   {
      dfl_modes2(status);
      error_root((status[1]<0)||((status[1]==0)&&(status[0]<0)),1,
                 "start_hmc [check2.c]","Deflation subspace generation "
                 "failed (status = %d;%d)",status[0],status[1]);

      if (status[1]==0)
         add2counter("modes",0,status);
      else
         add2counter("modes",2,status+1);
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
         dp=qlat_parms(ap.im0);
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
            error_root(1,1,"start_hmc [check2.c]","Unknown action");
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

   act1[0]=0.0;
   if((gauge()&1)!=0)
      act1[0]+=su3momentum_action(0);
   if((gauge()&2)!=0)
      act1[0]+=u1momentum_action(0);

   hmc=hmc_parms();
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
         dp=qlat_parms(ap.im0);
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


static void flip_mom(void)
{
   int status;
   su3_alg_dble *su3mom,*su3momm;
   double *u1mom,*u1momm;
   mdflds_t *mdfs;
   dfl_parms_t dfl;
   dfl_upd_parms_t dup;

   mdfs=mdflds();
   su3mom=(*mdfs).su3mom;
   su3momm=su3mom+4*VOLUME;

   for (;su3mom<su3momm;su3mom++)
   {
      (*su3mom).c1=-(*su3mom).c1;
      (*su3mom).c2=-(*su3mom).c2;
      (*su3mom).c3=-(*su3mom).c3;
      (*su3mom).c4=-(*su3mom).c4;
      (*su3mom).c5=-(*su3mom).c5;
      (*su3mom).c6=-(*su3mom).c6;
      (*su3mom).c7=-(*su3mom).c7;
      (*su3mom).c8=-(*su3mom).c8;
   }
   
   u1mom=(*mdfs).u1mom;
   u1momm=u1mom+4*VOLUME;

   for (;u1mom<u1momm;u1mom++)
   {
      (*u1mom)=-(*u1mom);
   }

   dfl=dfl_parms();

   if (dfl.Ns)
   {
      dup=dfl_upd_parms();
      dfl_update(dup.nsm,&status);
      error_root(status<0,1,"flip_mom [check2.c]",
                 "Deflation subspace update failed (status = %d)",status);
      add2counter("modes",1,&status);
   }
}


static double cmp_ud(su3_dble *u,su3_dble *v)
{
   int i;
   double r[18],dev,dmax;

   r[ 0]=(*u).c11.re-(*v).c11.re;
   r[ 1]=(*u).c11.im-(*v).c11.im;
   r[ 2]=(*u).c12.re-(*v).c12.re;
   r[ 3]=(*u).c12.im-(*v).c12.im;
   r[ 4]=(*u).c13.re-(*v).c13.re;
   r[ 5]=(*u).c13.im-(*v).c13.im;

   r[ 6]=(*u).c21.re-(*v).c21.re;
   r[ 7]=(*u).c21.im-(*v).c21.im;
   r[ 8]=(*u).c22.re-(*v).c22.re;
   r[ 9]=(*u).c22.im-(*v).c22.im;
   r[10]=(*u).c23.re-(*v).c23.re;
   r[11]=(*u).c23.im-(*v).c23.im;

   r[12]=(*u).c31.re-(*v).c31.re;
   r[13]=(*u).c31.im-(*v).c31.im;
   r[14]=(*u).c32.re-(*v).c32.re;
   r[15]=(*u).c32.im-(*v).c32.im;
   r[16]=(*u).c33.re-(*v).c33.re;
   r[17]=(*u).c33.im-(*v).c33.im;

   dmax=0.0;

   for (i=0;i<18;i+=2)
   {
      dev=r[i]*r[i]+r[i+1]*r[i+1];
      if (dev>dmax)
         dmax=dev;
   }

   return dmax;
}


static double max_dev_ud(su3_dble *v)
{
   double d,dmax;
   su3_dble *u,*um;

   u=udfld();
   um=u+4*VOLUME;
   dmax=0.0;

   for (;u<um;u++)
   {
      d=cmp_ud(u,v);

      if (d>dmax)
         dmax=d;

      v+=1;
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return sqrt(dmax);
}


static double max_dev_ad(double *v)
{
   double d,dmax;
   double *u,*um;

   u=adfld();
   um=u+4*VOLUME;
   dmax=0.0;

   for (;u<um;u++)
   {
      d=fabs((*u)-(*v));

      if (d>dmax)
         dmax=d;

      v+=1;
   }

   if (NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static void bcs_dev(double *dev)
{
   su3_dble *udb,**usv;
   double *adb,**asv;
   double d[4];
   
   if((gauge()&1)!=0)
   {
      udb=udfld();
      usv=reserve_wud(2);
      cm3x3_assign(4*VOLUME,udb,usv[1]);
      cm3x3_assign(4*VOLUME,udb,usv[0]);
      set_bc();
      d[0]=max_dev_ud(usv[0]);
      cm3x3_assign(4*VOLUME,udb,usv[0]);
      orbi_cpy_ud();
      d[1]=max_dev_ud(usv[0]);
      cm3x3_assign(4*VOLUME,usv[1],udb);
      release_wud();
   }
   if((gauge()&2)!=0)
   {
      adb=adfld();
      asv=reserve_wad(2);
      assign_dvec2dvec(4*VOLUME,adb,asv[1]);
      assign_dvec2dvec(4*VOLUME,adb,asv[0]);
      set_ad_bc();
      d[2]=max_dev_ad(asv[0]);
      assign_dvec2dvec(4*VOLUME,adb,asv[0]);
      orbi_cpy_ad();
      d[3]=max_dev_ad(asv[0]);
      assign_dvec2dvec(4*VOLUME,asv[1],adb);
      release_wad();
   }
   MPI_Allreduce(d,dev,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int first,last,step;
   int nc,nsize,icnfg,nact,i;
   int isap,idfl;
   int nwud,nwad,nws,nwsd,nwv,nwvd;
   double *act0,*act1,*act2;
   double sm0[2],sm1[2],dud[2],dH;
   double dudmin[2],dudmax[2],dudavg[2],dHmin,dHmax,dHavg;
   double bcd[4],bcdmax[4];
   su3_dble **usv;
   double **asv;
   hmc_parms_t hmc;
   char cnfg_dir[NAME_SIZE],cnfg_file[NAME_SIZE];
   char nbase[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;
   int cnfg_type;
   const char gg_str[3][256]={"SU(3)","U(1)","SU(3)xU(1)"};

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);

      printf("\n");
      printf("Reversibility of the MD evolution and stability of boundary conditions\n");
      printf("----------------------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      find_section("Configurations");
      read_line("cnfg_dir","%s",cnfg_dir);
      read_line("name","%s",nbase);
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);
   }

   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   read_flds_bc_lat_parms();
   read_hmc_parms();
   read_actions();
   read_integrator();
   read_solvers();

   if (my_rank==0)
      fclose(fin);

   hmc_wsize(&nwud,&nwad,&nws,&nwsd,&nwv,&nwvd);
   alloc_wud(3);
   alloc_wad(3);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);
   usv=reserve_wud(1);
   asv=reserve_wad(1);

   hmc=hmc_parms();
   nact=hmc.nact;
   act0=malloc(3*(nact+1)*sizeof(*act0));
   act1=act0+nact+1;
   act2=act1+nact+1;
   error(act0==NULL,1,"main [check2.c]","Unable to allocate action arrays");

   print_flds_bc_lat_parms();
   print_hmc_parms();
   print_action_parms();
   print_rat_parms();
   print_mdint_parms();
   print_force_parms2();
   print_solver_parms(&isap,&idfl);
   if (isap)
      print_sap_parms(0);
   if (idfl)
      print_dfl_parms(1);

   if (my_rank==0)
   {
      printf("Configurations %sn%d -> %sn%d in steps of %d\n\n",
             nbase,first,nbase,last,step);
      fflush(flog);
   }

   start_ranlux(0,1234);
   geometry();

   error_root(((last-first)%step)!=0,1,"main [check2.c]",
              "last-first is not a multiple of step");
   check_dir_root(cnfg_dir);

   nsize=name_size("%s/%sn%d",cnfg_dir,nbase,last);
   error_root(nsize>=NAME_SIZE,1,"main [check2.c]",
              "Configuration file name is too long");

   hmc_sanity_check();
   set_mdsteps();
   setup_counters();
   setup_chrono();

   dud[0]=dud[1]=0.0;
   dudmin[0]=dudmin[1]=0.0;
   dudmax[0]=dudmax[1]=0.0;
   dudavg[0]=dudavg[1]=0.0;
   dHmin=0.0;
   dHmax=0.0;
   dHavg=0.0;
   bcd[0]=bcd[1]=bcd[2]=bcd[3]=0.0;
   bcdmax[0]=bcdmax[1]=bcdmax[2]=bcdmax[3]=0.0;
   

   for (icnfg=first;icnfg<=last;icnfg+=step)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      cnfg_type=import_cnfg(cnfg_file);
      if (((gauge()&1)!=0)&&((cnfg_type&1)==0)) random_ud();
      if (((gauge()&2)!=0)&&((cnfg_type&2)==0)) random_ad();

      if (my_rank==0)
      {
         printf("Configuration type %s no %d\n",gg_str[cnfg_type-1],icnfg);
         fflush(flog);
      }
      
      bcs_dev(bcd);
      if (my_rank==0)
      {
         if((gauge()&1)!=0)
            printf("U time & C* bc dev = %.1e %.1e\n",bcdmax[0],bcdmax[1]);
         if((gauge()&2)!=0)
            printf("A time & C* bc dev = %.1e %.1e\n",bcdmax[2],bcdmax[3]);
         fflush(flog);
      }
      if(bcd[0]>bcdmax[0]) bcdmax[0]=bcd[0];
      if(bcd[1]>bcdmax[1]) bcdmax[1]=bcd[1];
      if(bcd[2]>bcdmax[2]) bcdmax[2]=bcd[2];
      if(bcd[3]>bcdmax[3]) bcdmax[3]=bcd[3];

      start_hmc(act0,usv[0],asv[0]);
      if((gauge()&1)!=0) dud[0]=max_dev_ud(usv[0]);
      if((gauge()&2)!=0) dud[1]=max_dev_ad(asv[0]);
      run_mdint();
      end_hmc(act1);

      sm0[0]=0.0;
      sm0[1]=0.0;

      for (i=0;i<=nact;i++)
      {
         sm0[0]+=act0[i];
         sm0[1]+=(act1[i]-act0[i]);
      }

      MPI_Reduce(sm0,sm1,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(sm1,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         printf("start_hmc:\n");
         if((gauge()&1)!=0) printf("max|U_ij-U'_ij| = %.1e\n",dud[0]);
         if((gauge()&2)!=0) printf("max|A_ij-A'_ij| = %.1e\n",dud[1]);
         printf("run_mdint:\n");
         printf("H = %.6e\n",sm1[0]);
         printf("dH = %.2e\n",sm1[1]);
         fflush(flog);
      }
      
      bcs_dev(bcd);
      if (my_rank==0)
      {
         if((gauge()&1)!=0)
            printf("U time & C* bc dev = %.1e %.1e\n",bcdmax[0],bcdmax[1]);
         if((gauge()&2)!=0)
            printf("A time & C* bc dev = %.1e %.1e\n",bcdmax[2],bcdmax[3]);
         fflush(flog);
      }
      if(bcd[0]>bcdmax[0]) bcdmax[0]=bcd[0];
      if(bcd[1]>bcdmax[1]) bcdmax[1]=bcd[1];
      if(bcd[2]>bcdmax[2]) bcdmax[2]=bcd[2];
      if(bcd[3]>bcdmax[3]) bcdmax[3]=bcd[3];

      print_all_avgstat();

      flip_mom();
      run_mdint();
      end_hmc(act2);

      sm0[0]=0.0;
      sm0[1]=0.0;

      for (i=0;i<=nact;i++)
      {
         sm0[0]+=act2[i];
         sm0[1]+=(act2[i]-act0[i]);
      }

      MPI_Reduce(sm0,sm1,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(sm1,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      dH=fabs(sm1[1]);
      if((gauge()&1)!=0) dud[0]=max_dev_ud(usv[0]);
      if((gauge()&2)!=0) dud[1]=max_dev_ad(asv[0]);

      if (my_rank==0)
      {
         printf("Flip momenta and run_mdint:\n");
         printf("H  = %.6e\n",sm1[0]);
         printf("|dH| = % .2e\n",dH);
         if((gauge()&1)!=0) printf("max|U_ij-U'_ij| = %.2e\n",dud[0]);
         if((gauge()&2)!=0) printf("max|A_ij-A'_ij| = %.2e\n",dud[1]);
         fflush(flog);
      }
      
      bcs_dev(bcd);
      if (my_rank==0)
      {
         if((gauge()&1)!=0)
            printf("U time & C* bc dev = %.1e %.1e\n",bcdmax[0],bcdmax[1]);
         if((gauge()&2)!=0)
            printf("A time & C* bc dev = %.1e %.1e\n",bcdmax[2],bcdmax[3]);
         fflush(flog);
      }
      if(bcd[0]>bcdmax[0]) bcdmax[0]=bcd[0];
      if(bcd[1]>bcdmax[1]) bcdmax[1]=bcd[1];
      if(bcd[2]>bcdmax[2]) bcdmax[2]=bcd[2];
      if(bcd[3]>bcdmax[3]) bcdmax[3]=bcd[3];

      if (icnfg==first)
      {
         dudmin[0]=dud[0];
         dudmax[0]=dud[0];
         dudavg[0]=dud[0];
         
         dudmin[1]=dud[1];
         dudmax[1]=dud[1];
         dudavg[1]=dud[1];

         dHmin=dH;
         dHmax=dH;
         dHavg=dH;
      }
      else
      {
         if (dud[0]<dudmin[0])
            dudmin[0]=dud[0];
         if (dud[0]>dudmax[0])
            dudmax[0]=dud[0];
         dudavg[0]+=dud[0];

         if (dud[1]<dudmin[1])
            dudmin[1]=dud[1];
         if (dud[1]>dudmax[1])
            dudmax[1]=dud[1];
         dudavg[1]+=dud[1];

         if (dH<dHmin)
            dHmin=dH;
         if (dH>dHmax)
            dHmax=dH;
         dHavg+=dH;
      }
   }

   if (my_rank==0)
   {
      nc=(last-first)/step+1;

      printf("\nTest summary\n");
      printf("------------\n\n");

      printf("Considered %d configurations in the range %d -> %d\n\n",
             nc,first,last);

      printf("Check boundary conditions, maximal deviations are quoted\n\n");
      if((gauge()&1)!=0)
         printf("U time & C* bc dev = %.1e %.1e\n",bcdmax[0],bcdmax[1]);
      if((gauge()&2)!=0)
         printf("A time & C* bc dev = %.1e %.1e\n",bcdmax[2],bcdmax[3]);

      printf("\nThe three figures quoted in each case are the minimal,\n");
      printf("maximal and average values\n\n");

      if((gauge()&1)!=0) printf("max|U_ij-U'_ij| = %.2e, %.2e, %.2e\n",
             dudmin[0],dudmax[0],dudavg[0]/(double)(nc));
      if((gauge()&2)!=0) printf("max|A_ij-A'_ij| = %.2e, %.2e, %.2e\n",
             dudmin[1],dudmax[1],dudavg[1]/(double)(nc));
      printf("|dH|            = %.2e, %.2e, %.2e\n\n",
             dHmin,dHmax,dHavg/(double)(nc));

      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
