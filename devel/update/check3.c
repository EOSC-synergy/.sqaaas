
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2007, 2009-2013, Martin Luescher, Filippo Palombi,
*               2016                   Stefan Schaefer,
*               2017, 2019             Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Conservation of the Hamilton function by the MD evolution.
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
#include "u1flds.h"
#include "mdflds.h"
#include "linalg.h"
#include "archive.h"
#include "forces.h"
#include "dfl.h"
#include "update.h"
#include "global.h"

static int my_rank;


static void read_flds_bc_lat_parms(void)
{
   int gg,nfl;
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

   read_bc_parms();
   read_glat_parms();
   read_qlat_parms();
}


static void read_hmc_parms(void)
{
   int nact,*iact;
   int npf,nmu,nlv,facc;
   double tau,*mu;

   facc=0;
   if (my_rank==0)
   {
      find_section("HMC parameters");
      if ((gauge()&2)!=0)
         read_line("facc","%d",&facc);
      nact=count_tokens("actions");
      read_line("npf","%d",&npf);
      nmu=count_tokens("mu");
      read_line("nlv","%d",&nlv);
      read_line("tau","%lf",&tau);
   }

   MPI_Bcast(&facc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nact,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&npf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmu,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nlv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&tau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   if (nact>0)
   {
      iact=malloc(nact*sizeof(*iact));
      error(iact==NULL,1,"read_hmc_parms [check3.c]",
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
      error(mu==NULL,1,"read_hmc_parms [check3.c]",
            "Unable to allocate temporary array");
      if (my_rank==0)
         read_dprms("mu",nmu,mu);
      MPI_Bcast(mu,nmu,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
      mu=NULL;

   set_hmc_parms(nact,iact,npf,nmu,mu,nlv,tau,facc);

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
                  
                  if (dfl_gen_parms(sp.idfl).status!=DFL_DEF)
                     read_dfl_parms(sp.idfl);
               }
            }
         }
      }
   }

   if (isap)
      read_sap_parms();

   if (idfl)
      read_dfl_parms(-2);
}


static void chk_mode_regen(int isp,int *status)
{
   solver_parms_t sp;

   sp=solver_parms(isp);

   if ((sp.solver==DFL_SAP_GCR)&&(status[2]>0))
      add2counter("modes",2,status+2);
}


static void start_hmc(double *act0,su3_dble *uold,double *aold,su3_alg_dble *su3mold,double *u1mold)
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
      random_su3mom();
      assign_alg2alg(4*VOLUME,mdflds()->su3mom,su3mold);
      act0[0]+=su3momentum_action(0);
   }
   if((gauge()&2)!=0)
   {
      adb=adfld();
      assign_dvec2dvec(4*VOLUME,adb,aold);
      random_u1mom();
      act0[0]+=u1momentum_action(0);
      if (hmc.facc)
         u1mom_Delta_no0(2,mdflds()->u1mom,mdflds()->u1mom);
      assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,u1mold);
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
            dfl_modes(idfl,status);
            error_root(status[0]<0,1,"start_hmc [check3.c]",
                       "Generation of eflation subspace %d failed (status = %d)",
                       idfl,status[0]);
            add2counter("modes",0,status);
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
            error_root(1,1,"start_hmc [check3.c]","Unknown action");
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
      if (hmc.facc)
         u1mom_Delta_no0(2,mdflds()->u1mom,mdflds()->u1mom);
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


static void restart_hmc(su3_dble *uold,double *aold,su3_alg_dble *su3mold,double *u1mold)
{
   int status,idfl;
   su3_dble *udb;
   double *adb;
   dfl_parms_t dfl;
   dflst_t dfl_status;

   clear_counters();

   if((gauge()&1)!=0)
   {
      udb=udfld();
      cm3x3_assign(4*VOLUME,uold,udb);
      assign_alg2alg(4*VOLUME,su3mold,mdflds()->su3mom);
      set_flags(UPDATED_UD);
   }
   if((gauge()&2)!=0)
   {
      adb=adfld();
      assign_dvec2dvec(4*VOLUME,aold,adb);
      assign_dvec2dvec(4*VOLUME,u1mold,mdflds()->u1mom);
      set_flags(UPDATED_AD);
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
            dfl_modes(idfl,&status);
            error_root(status<0,1,"start_hmc [check3.c]",
                       "Generation of eflation subspace %d failed (status = %d)",
                       idfl,status);
            add2counter("modes",0,&status);
         }
         idfl++;
      }
   }
}


int main(int argc,char *argv[])
{
   int first,last,step;
   int nsize,icnfg,nact,i,j;
   int isap,idfl;
   int nwud,nwad,nws,nwsd,nwv,nwvd;
   double *act0,*act1,tau[4];
   double sm0[3],sm1[3],dH[4];
   su3_dble **usv;
   su3_alg_dble **f3sv;
   double **asv;
   double **f1sv;
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
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Conservation of the Hamilton function by the MD evolution\n");
      printf("---------------------------------------------------------\n\n");

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
   alloc_wud(1);
   alloc_wad(1);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);
   alloc_wf3d(1);
   alloc_wf1d(1);
   usv=reserve_wud(1);
   asv=reserve_wad(1);
   f3sv=reserve_wf3d(1);
   f1sv=reserve_wf1d(1);

   hmc=hmc_parms();
   nact=hmc.nact;
   act0=malloc(2*(nact+1)*sizeof(*act0));
   act1=act0+nact+1;
   error(act0==NULL,1,"main [check3.c]","Unable to allocate action arrays");
   tau[0]=hmc.tau;

   for (i=1;i<4;i++)
      tau[i]=tau[i-1]/pow(4.0,1.0/3.0);

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

   error_root(((last-first)%step)!=0,1,"main [check3.c]",
              "last-first is not a multiple of step");
   check_dir_root(cnfg_dir);

   nsize=name_size("%s/%sn%d",cnfg_dir,nbase,last);
   error_root(nsize>=NAME_SIZE,1,"main [check3.c]",
              "Configuration file name is too long");

   hmc_sanity_check();
   setup_counters();
   setup_chrono();

   for (icnfg=first;icnfg<=last;icnfg+=step)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      cnfg_type=import_cnfg(cnfg_file);
      error_root(((~gauge())&&cnfg_type)==0,1,"main [check3.c]",
                 "Imported inactive gauge field");
      if (((gauge()&1)!=0)&&((cnfg_type&1)==0)) random_ud();
      if (((gauge()&2)!=0)&&((cnfg_type&2)==0)) random_ad();

      if (my_rank==0)
      {
         printf("Configuration type %s no %d\n",gg_str[cnfg_type-1],icnfg);
         fflush(flog);
      }

      for (i=0;i<4;i++)
      {
         set_hmc_parms(hmc.nact,hmc.iact,hmc.npf,
                       hmc.nmu,hmc.mu,hmc.nlv,tau[i],hmc.facc);
         set_mdsteps();
         
         if (i==0)
            start_hmc(act0,usv[0],asv[0],f3sv[0],f1sv[0]);
         else
            restart_hmc(usv[0],asv[0],f3sv[0],f1sv[0]);
         
         run_mdint();
         end_hmc(act1);

         sm0[0]=0.0;
         sm0[1]=0.0;
         sm0[2]=0.0;

         for (j=0;j<=nact;j++)
         {
            sm0[0]+=act0[j];
            sm0[1]+=act1[j];
            sm0[2]+=(act1[j]-act0[j]);
         }

         MPI_Reduce(sm0,sm1,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         MPI_Bcast(sm1,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
         dH[i]=fabs(sm1[2]);

         if (my_rank==0)
         {
            if (i==0)
            {
               printf("start_hmc:\n");
               printf("H = %.6e\n",sm1[0]);
               fflush(flog);
            }

            printf("run_md:\n");
            printf("tau = %.3f\n",tau[i]);
            printf("H = %.6e, |dH| = %.2e\n",sm1[1],dH[i]);
            fflush(flog);
         }

         print_all_avgstat();
      }

      if (my_rank==0)
      {
         printf("\n");
         printf("tau = %.2e, |dH| = %.2e\n",tau[0],dH[0]);

         for (i=1;i<4;i++)
         {
            printf("tau = %.2e, |dH| = %.2e, |dH[i]|/|dH[i-1]| = %.2e\n",
                   tau[i],dH[i],dH[i]/dH[i-1]);
         }

         printf("\n");
         printf("(From one tau to the next, the scale factor s is 4^(1/3),\n"
                "i.e. s^{-3,-4,-5} = {%.2e,%.2e,%.2e})\n\n",
                pow(4.0,-1.0),pow(4.0,-4.0/3.0),pow(4.0,-5.0/3.0));
         fflush(flog);
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
