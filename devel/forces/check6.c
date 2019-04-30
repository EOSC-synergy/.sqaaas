
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*               2016 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of force1() and action1().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dfl.h"
#include "forces.h"
#include "global.h"
#include "gflds_utils.h"
#include "mdflds_utils.h"

#define N0 (NPROC0*L0)


static void unconstrained_random_su3mom(void)
{
   random_alg(4*VOLUME,mdflds()->su3mom);
   bnd_su3mom2zero();
}


static void unconstrained_random_u1mom(void)
{
   random_dvec(4*VOLUME,mdflds()->u1mom);
   bnd_u1mom2zero();
}


static double su3dSdt(double mu,int ipf,int isp,int *status)
{
   mdflds_t *mdfs;

   mdfs=mdflds();
   set_su3frc2zero();
   force1(mu,ipf,isp,0,1.2345,status);
   check_bnd_su3frc((*mdfs).su3frc);

   return scalar_prod_alg(4*VOLUME,0,(*mdfs).su3mom,(*mdfs).su3frc);
}


static double u1dSdt(double mu,int ipf,int isp,int *status)
{
   mdflds_t *mdfs;
   
   mdfs=mdflds();
   set_u1frc2zero();
   force1(mu,ipf,isp,0,1.2345,status);
   check_bnd_u1frc((*mdfs).u1frc);

   return scalar_prod_dvec(4*VOLUME,0,(*mdfs).u1mom,(*mdfs).u1frc);
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,isp,status[6],mnkv;
   int bs[4],Ns,nmx,nkv,nmr,ncy,ninv,qhat;
   int isap,idfl;
   double chi[2],chi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   double kappa,mu,res;
   double eps,act0,act1,dact,dsdt;
   double dev_act[2],dev_frc,sig_loss,rdmy;
   solver_parms_t sp;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);
      fin=freopen("check6.in","r",stdin);

      printf("\n");
      printf("Check of force1() and action1()\n");
      printf("-------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-cs <cstar>]");
   }
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   chi[0]=0.123;
   chi[1]=-0.534;
   chi_prime[0]=0.912;
   chi_prime[1]=0.078;
   set_bc_parms(bc,0,cs,chi,chi_prime);
   print_bc_parms();

   set_su3lat_parms(3.50,0.95,0.82,1.32);
   set_u1lat_parms(0,1.5,1.2,0.0,0.482,0.87,0.57);
   print_lat_parms();

   if (my_rank==0)
   {
      find_section("SAP");
      read_iprms("bs",4,bs);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   set_sap_parms(bs,1,4,5);

   if (my_rank==0)
   {
      find_section("Deflation subspace");
      read_iprms("bs",4,bs);
      read_line("Ns","%d",&Ns);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_parms(bs,Ns);

   if (my_rank==0)
   {
      find_section("Deflation subspace generation");
      read_line("qhat","%d",&qhat);
      read_line("kappa","%lf",&kappa);
      read_line("mu","%lf",&mu);
      read_line("ninv","%d",&ninv);
      read_line("nmr","%d",&nmr);
      read_line("ncy","%d",&ncy);
   }

   MPI_Bcast(&qhat,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_gen_parms(kappa,mu,qhat,0.95,0.8,1.301,0.789,0.38,-1.25,0.54,ninv,nmr,ncy);

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

   set_hmc_parms(0,NULL,1,0,NULL,1,1.0);
   mnkv=0;

   for (isp=0;isp<3;isp++)
   {
      read_solver_parms(isp);
      sp=solver_parms(isp);

      if (sp.nkv>mnkv)
         mnkv=sp.nkv;
   }

   if (my_rank==0)
      fclose(fin);

   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   if ((gauge()&1)!=0) su3csw=0.95;
   if ((gauge()&2)!=0) u1csw=0.8;
   if (bc_type()!=3)
   {
      cF[0]=1.301;
      cF[1]=0.789;
   }
   if (bc_cstar()==0)
   {
      theta[0]=0.35;
      theta[1]=-1.25;
      theta[2]=0.78;
   }
   set_dirac_parms9(-3,-0.0123,su3csw,u1csw,cF[0],cF[1],
                    theta[0],theta[1],theta[2]);
   
   print_dirac_parms();
   print_solver_parms(&isap,&idfl);
   print_sap_parms(1);
   print_dfl_parms(0);

   start_ranlux(0,1245);
   geometry();

   mnkv=2*mnkv+2;
   if (mnkv<(Ns+2))
      mnkv=Ns+2;
   if (mnkv<5)
      mnkv=5;

   alloc_ws(mnkv);
   alloc_wsd(6);
   alloc_wv(2*nkv+2);
   alloc_wvd(4);
   
   if (my_rank==0)
      printf("++++++ SU(3) forces\n\n");

   for (isp=0;isp<3;isp++)
   {
      if (isp==0)
      {
         mu=1.0;
         eps=1.0e-4;
      }
      else if (isp==1)
      {
         mu=0.1;
         eps=2.0e-4;
      }
      else
      {
         mu=0.01;
         eps=3.0e-4;
      }

      random_gflds();
      unconstrained_random_su3mom();

      if (isp==2)
      {
         dfl_modes(status);
         error_root(status[0]<0,1,"main [check6.c]",
                    "dfl_modes failed");
      }

      status[0]=0;
      status[1]=0;

      act0=setpf1(mu,0,0);
      act1=action1(mu,0,isp,0,status);
      error_root((status[0]<0)||(status[1]<0),1,"main [check6.c]",
                 "action1 failed (mu = %.2e, isp=%d)",mu,isp);

      rdmy=fabs(act1-act0);
      MPI_Reduce(&rdmy,dev_act,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      rdmy=act1-act0;
      MPI_Reduce(&rdmy,dev_act+1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(dev_act,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      dsdt=su3dSdt(mu,0,isp,status);

      if (my_rank==0)
      {
         printf("Solver number %d, mu = %.2e\n",isp,mu);

         if (isp==0)
            printf("Status = %d\n",status[0]);
         else if (isp==1)
            printf("Status = %d,%d\n",status[0],status[1]);
         else
            printf("Status = (%d,%d,%d),(%d,%d,%d)\n",
                   status[0],status[1],status[2],status[3],
                   status[4],status[5]);

         printf("Absolute action difference |setpf1-action1| = %.1e,",
                fabs(dev_act[1]));
         printf(" %.1e (local)\n",dev_act[0]);
         fflush(flog);
      }

      rot_ud(eps);
      act0=2.0*action1(mu,0,isp,0,status)/3.0;
      rot_ud(-eps);

      rot_ud(-eps);
      act1=2.0*action1(mu,0,isp,0,status)/3.0;
      rot_ud(eps);

      rot_ud(2.0*eps);
      act0-=action1(mu,0,isp,0,status)/12.0;
      rot_ud(-2.0*eps);

      rot_ud(-2.0*eps);
      act1-=action1(mu,0,isp,0,status)/12.0;
      rot_ud(2.0*eps);

      dact=1.2345*(act0-act1)/eps;
      dev_frc=dsdt-dact;
      sig_loss=-log10(fabs(1.0-act0/act1));

      rdmy=dsdt;
      MPI_Reduce(&rdmy,&dsdt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dsdt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=dev_frc;
      MPI_Reduce(&rdmy,&dev_frc,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dev_frc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=sig_loss;
      MPI_Reduce(&rdmy,&sig_loss,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&sig_loss,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         printf("Relative deviation of dS/dt = %.2e ",fabs(dev_frc/dsdt));
         printf("[significance loss = %d digits]\n\n",(int)(sig_loss));
         fflush(flog);
      }
   }
   
   if (my_rank==0)
      printf("\n++++++ U(1) forces\n\n");

   for (isp=0;isp<3;isp++)
   {
      if (isp==0)
      {
         mu=1.0;
         eps=1.0e-4;
      }
      else if (isp==1)
      {
         mu=0.1;
         eps=2.0e-4;
      }
      else
      {
         mu=0.01;
         eps=3.0e-4;
      }

      random_gflds();
      unconstrained_random_u1mom();

      if (isp==2)
      {
         dfl_modes(status);
         error_root(status[0]<0,1,"main [check6.c]",
                    "dfl_modes failed");
      }

      status[0]=0;
      status[1]=0;

      act0=setpf1(mu,0,0);
      act1=action1(mu,0,isp,0,status);
      error_root((status[0]<0)||(status[1]<0),1,"main [check6.c]",
                 "action1 failed (mu = %.2e, isp=%d)",mu,isp);

      rdmy=fabs(act1-act0);
      MPI_Reduce(&rdmy,dev_act,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      rdmy=act1-act0;
      MPI_Reduce(&rdmy,dev_act+1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(dev_act,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      dsdt=u1dSdt(mu,0,isp,status);

      if (my_rank==0)
      {
         printf("Solver number %d, mu = %.2e\n",isp,mu);

         if (isp==0)
            printf("Status = %d\n",status[0]);
         else if (isp==1)
            printf("Status = %d,%d\n",status[0],status[1]);
         else
            printf("Status = (%d,%d,%d),(%d,%d,%d)\n",
                   status[0],status[1],status[2],status[3],
                   status[4],status[5]);

         printf("Absolute action difference |setpf1-action1| = %.1e,",
                fabs(dev_act[1]));
         printf(" %.1e (local)\n",dev_act[0]);
         fflush(flog);
      }

      rot_ad(eps);
      act0=2.0*action1(mu,0,isp,0,status)/3.0;
      rot_ad(-eps);

      rot_ad(-eps);
      act1=2.0*action1(mu,0,isp,0,status)/3.0;
      rot_ad(eps);

      rot_ad(2.0*eps);
      act0-=action1(mu,0,isp,0,status)/12.0;
      rot_ad(-2.0*eps);

      rot_ad(-2.0*eps);
      act1-=action1(mu,0,isp,0,status)/12.0;
      rot_ad(2.0*eps);

      dact=1.2345*(act0-act1)/eps;
      dev_frc=dsdt-dact;
      sig_loss=-log10(fabs(1.0-act0/act1));

      rdmy=dsdt;
      MPI_Reduce(&rdmy,&dsdt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dsdt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=dev_frc;
      MPI_Reduce(&rdmy,&dev_frc,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dev_frc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=sig_loss;
      MPI_Reduce(&rdmy,&sig_loss,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&sig_loss,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      if (my_rank==0)
      {
         printf("Relative deviation of dS/dt = %.2e ",fabs(dev_frc/dsdt));
         printf("[significance loss = %d digits]\n\n",(int)(sig_loss));
         fflush(flog);
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
