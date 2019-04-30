
/*******************************************************************************
*
* File mdint.c
*
* Copyright (C) 2011-2013 Stefan Schaefer, Martin Luescher, John Bulava
*               2017      Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Integration of the molecular-dynamics equations.
*
* The externally accessible functions are
*
*   void run_mdint(void)
*     Integrates the molecular-dynamics equations using the current
*     integrator (see the notes).
*
* Notes:
*
* The integrator used is the one defined by the array of elementary operations
* returned by mdsteps() (see update/mdsteps.c). It is assumed that the fields
* and the integrator have been properly initialized.
*
* In the course of the integration, the solver iteration numbers are added
* to the appropriate counters provided by the module update/counters.c.
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes.
*
* Some debugging information is printed to stdout if the macro MDINT_DBG is
* defined. The norm of the forces printed is the norm per active link.
*
*******************************************************************************/

#define MDINT_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "mdflds.h"
#include "su3fcts.h"
#include "linalg.h"
#include "dfl.h"
#include "forces.h"
#include "update.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int nsm;
static double rtau,dtau;


static void chk_mode_regen(int isp,int *status)
{
   int i,is;
   solver_parms_t sp;

   sp=solver_parms(isp);

   if (sp.solver==DFL_SAP_GCR)
   {
      is=status[2];

      for (i=2;i<4;i++)
         status[i]=status[i+1];

      status[4]=is;

      if (status[4]>0)
         add2counter("modes",2,status+4);
      if (status[5]>0)
         add2counter("modes",2,status+5);
   }
}


static void update_su3mom(int isym)
{
   int bc,cs,ix,t,ifc;
   int mirror, tag;
   su3_alg_dble *mom,*frc;
   mdflds_t *mdfs;
   MPI_Status stat;

   bc=bc_type();
   cs=bc_cstar();
   mdfs=mdflds();
   
   if ((cs!=0)&&(isym!=0))
   {
      mom=(*mdfs).su3mom;
      for (ix=0;ix<4*VOLUME;ix++)
      {
         _su3_alg_mul_assign((*mom),0.5);
         mom+=1;
      }
   }

   if ((cs!=0)&&(isym==0))
   {
      mom=(*mdfs).su3mom;
      frc=(*mdfs).su3frc;
      for (ix=(VOLUME/2);ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (t==0)
         {
            _su3_alg_mul_add_assign((*mom),-2.0,(*frc));
            mom+=1;
            frc+=1;

            if (bc!=0)
            {
               _su3_alg_mul_add_assign((*mom),-2.0,(*frc));
            }

            mom+=1;
            frc+=1;

            for (ifc=2;ifc<8;ifc++)
            {
               if (bc!=1)
               {
                  _su3_alg_mul_add_assign((*mom),-2.0,(*frc));
               }

               mom+=1;
               frc+=1;
            }
         }
         else if (t==(N0-1))
         {
            if (bc!=0)
            {
               _su3_alg_mul_add_assign((*mom),-2.0,(*frc));
            }

            mom+=1;
            frc+=1;

            for (ifc=1;ifc<8;ifc++)
            {
               _su3_alg_mul_add_assign((*mom),-2.0,(*frc));
               mom+=1;
               frc+=1;
            }
         }
         else
         {
            for (ifc=0;ifc<8;ifc++)
            {
               _su3_alg_mul_add_assign((*mom),-2.0,(*frc));
               mom+=1;
               frc+=1;
            }
         }
      }
   }
   else
   {
      mom=(*mdfs).su3mom;
      frc=(*mdfs).su3frc;
      for (ix=(VOLUME/2);ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (t==0)
         {
            _su3_alg_sub_assign((*mom),(*frc));
            mom+=1;
            frc+=1;

            if (bc!=0)
            {
               _su3_alg_sub_assign((*mom),(*frc));
            }

            mom+=1;
            frc+=1;

            for (ifc=2;ifc<8;ifc++)
            {
               if (bc!=1)
               {
                  _su3_alg_sub_assign((*mom),(*frc));
               }

               mom+=1;
               frc+=1;
            }
         }
         else if (t==(N0-1))
         {
            if (bc!=0)
            {
               _su3_alg_sub_assign((*mom),(*frc));
            }

            mom+=1;
            frc+=1;

            for (ifc=1;ifc<8;ifc++)
            {
               _su3_alg_sub_assign((*mom),(*frc));
               mom+=1;
               frc+=1;
            }
         }
         else
         {
            for (ifc=0;ifc<8;ifc++)
            {
               _su3_alg_sub_assign((*mom),(*frc));
               mom+=1;
               frc+=1;
            }
         }
      }
   }

   if ((cs!=0)&&(isym!=0))
   {
      mirror=get_mirror_rank();
      tag=mpi_tag();

      mom=mdflds()->su3mom;
      frc=mdflds()->su3frc;

      MPI_Sendrecv(mom,8*4*VOLUME,MPI_DOUBLE,mirror,tag,
                   frc,8*4*VOLUME,MPI_DOUBLE,mirror,tag,
                   MPI_COMM_WORLD,&stat);

      for (ix=0;ix<4*VOLUME;ix++)
      {
         (*mom).c1-=(*frc).c1;
         (*mom).c2-=(*frc).c2;
         (*mom).c3+=(*frc).c3;
         (*mom).c4-=(*frc).c4;
         (*mom).c5+=(*frc).c5;
         (*mom).c6-=(*frc).c6;
         (*mom).c7+=(*frc).c7;
         (*mom).c8-=(*frc).c8;
         mom+=1;
         frc+=1;
      }
   }
}


static void update_u1mom(int isym)
{
   int bc,cs,ix,t,ifc;
   int mirror, tag;
   double *mom,*frc;
   mdflds_t *mdfs;
   MPI_Status stat;

   bc=bc_type();
   cs=bc_cstar();
   mdfs=mdflds();
   
   if ((cs!=0)&&(isym!=0))
   {
      mom=(*mdfs).u1mom;
      for (ix=0;ix<4*VOLUME;ix++)
      {
         (*mom)*=0.5;
         mom+=1;
      }
   }

   if ((cs!=0)&&(isym==0))
   {
      mom=(*mdfs).u1mom;
      frc=(*mdfs).u1frc;
      for (ix=(VOLUME/2);ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (t==0)
         {
            (*mom)-=2.0*(*frc);
            mom+=1;
            frc+=1;

            if (bc!=0)
            {
               (*mom)-=2.0*(*frc);
            }

            mom+=1;
            frc+=1;

            for (ifc=2;ifc<8;ifc++)
            {
               if (bc!=1)
               {
                  (*mom)-=2.0*(*frc);
               }

               mom+=1;
               frc+=1;
            }
         }
         else if (t==(N0-1))
         {
            if (bc!=0)
            {
               (*mom)-=2.0*(*frc);
            }

            mom+=1;
            frc+=1;

            for (ifc=1;ifc<8;ifc++)
            {
               (*mom)-=2.0*(*frc);
               mom+=1;
               frc+=1;
            }
         }
         else
         {
            for (ifc=0;ifc<8;ifc++)
            {
               (*mom)-=2.0*(*frc);
               mom+=1;
               frc+=1;
            }
         }
      }
   }
   else
   {
      mom=(*mdfs).u1mom;
      frc=(*mdfs).u1frc;
      for (ix=(VOLUME/2);ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (t==0)
         {
            (*mom)-=(*frc);
            mom+=1;
            frc+=1;

            if (bc!=0)
            {
               (*mom)-=(*frc);
            }

            mom+=1;
            frc+=1;

            for (ifc=2;ifc<8;ifc++)
            {
               if (bc!=1)
               {
                  (*mom)-=(*frc);
               }

               mom+=1;
               frc+=1;
            }
         }
         else if (t==(N0-1))
         {
            if (bc!=0)
            {
               (*mom)-=(*frc);
            }

            mom+=1;
            frc+=1;

            for (ifc=1;ifc<8;ifc++)
            {
               (*mom)-=(*frc);
               mom+=1;
               frc+=1;
            }
         }
         else
         {
            for (ifc=0;ifc<8;ifc++)
            {
               (*mom)-=(*frc);
               mom+=1;
               frc+=1;
            }
         }
      }
   }

   if ((cs!=0)&&(isym!=0))
   {
      mirror=get_mirror_rank();
      tag=mpi_tag();

      mom=mdflds()->u1mom;
      frc=mdflds()->u1frc;

      MPI_Sendrecv(mom,4*VOLUME,MPI_DOUBLE,mirror,tag,
                   frc,4*VOLUME,MPI_DOUBLE,mirror,tag,
                   MPI_COMM_WORLD,&stat);

      for (ix=0;ix<4*VOLUME;ix++)
      {
         (*mom)-=(*frc);
         mom+=1;
         frc+=1;
      }
   }
}


static void update_ud(double eps)
{
   int bc,ix,t,ifc;
   su3_dble *u;
   su3_alg_dble *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).su3mom;
   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         if (bc!=0)
            expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static void update_ad(double eps)
{
   int bc,ix,t,ifc;
   double *a;
   double *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).u1mom;
   a=adfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         (*a)+=eps*(*mom);
         a+=1;
         mom+=1;

         if (bc!=0)
            (*a)+=eps*(*mom);
         a+=1;
         mom+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               (*a)+=eps*(*mom);
            a+=1;
            mom+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            (*a)+=eps*(*mom);
         a+=1;
         mom+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            (*a)+=eps*(*mom);
            a+=1;
            mom+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            (*a)+=eps*(*mom);
            a+=1;
            mom+=1;
         }
      }
   }

   set_flags(UPDATED_AD);
}


static void start_dfl_upd(void)
{
   dfl_upd_parms_t dup;

   dup=dfl_upd_parms();
   dtau=dup.dtau;
   nsm=dup.nsm;
   rtau=0.0;
}


static void dfl_upd(int isp)
{
   int status[2];
   solver_parms_t sp;

   if ((nsm>0)&&(rtau>dtau))
   {
      sp=solver_parms(isp);

      if (sp.solver==DFL_SAP_GCR)
      {
         dfl_update2(nsm,status);
         error_root((status[1]<0)||((status[1]==0)&&(status[0]<0)),1,
                    "dfl_upd [mdint.c]","Deflation subspace update "
                    "failed (status = %d;%d)",status[0],status[1]);

         if (status[1]==0)
            add2counter("modes",1,status);
         else
            add2counter("modes",2,status+1);

         rtau=0.0;
      }
   }
}

#ifdef MDINT_DBG

void run_mdint(void)
{
   int my_rank,nop;
   int iop,status[6],isym;
   double *mu,eps,nlk,nrm;
   mdflds_t *mdfs;
   mdstep_t *s,*sm;
   hmc_parms_t hmc;
   force_parms_t fp;
   dirac_parms_t dp;
   double wt1, wt2;
   double u1mdt,su3mdt,dt;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   mdfs=mdflds();
   hmc=hmc_parms();
   mu=hmc.mu;
   u1mdt=su3mdt=0.0;
   reset_chrono();
   start_dfl_upd();

   nlk=(double)(4*N0*N1)*(double)(N2*N3);
   if (bc_type()==0)
      nlk-=(double)(N1)*(double)(N2*N3);
   else if (bc_type()==1)
      nlk-=(double)(3*N1)*(double)(N2*N3);

   s=mdsteps(&nop);
   sm=s+nop;

   for (;s<sm;s++)
   {
      iop=(*s).iop;
      eps=(*s).eps;

      if (iop>=0)
      {
         fp=force_parms(iop);

         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();

         if (fp.force==FRG_SU3)
         {
            isym=0;
            force0(eps);
            if ((gauge()&2)!=0) set_u1frc2zero();
         }
         else if (fp.force==FRG_U1)
         {
            isym=0;
            force6(eps);
            if ((gauge()&1)!=0) set_su3frc2zero();
         }
         else
         {
            isym=1;
            dfl_upd(fp.isp[0]);
            dp=qlat_parms(fp.im0);
            set_dirac_parms1(&dp);
            set_su3frc2zero();
            set_u1frc2zero();
            status[2]=0;
            status[5]=0;

            if (fp.force==FRF_TM1)
               force1(mu[fp.imu[0]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO)
               force4(mu[fp.imu[0]],fp.ipf,0,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO_SDET)
               force4(mu[fp.imu[0]],fp.ipf,1,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2)
               force2(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2_EO)
               force5(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_RAT)
               force3(fp.irat,fp.ipf,0,fp.isp[0],
                      eps,status);
            else if (fp.force==FRF_RAT_SDET)
               force3(fp.irat,fp.ipf,1,fp.isp[0],
                      eps,status);

            chk_mode_regen(fp.isp[0],status);
            add2counter("force",iop,status);
         }

         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();

         if((gauge()&1)!=0)
         {
            update_su3mom(isym);
            nrm=norm_square_alg(4*VOLUME,1,(*mdfs).su3frc);
            nrm=sqrt(nrm/nlk);

            if (my_rank==0)
            {
               if (fp.force==FRG_SU3)
                  printf("SU(3) Force FRG_SU3:          ");
               else if (fp.force==FRF_TM1)
                  printf("SU(3) Force FRF_TM1:          ");
               else if (fp.force==FRF_TM1_EO)
                  printf("SU(3) Force FRF_TM1_EO:       ");
               else if (fp.force==FRF_TM1_EO_SDET)
                  printf("SU(3) Force FRF_TM1_EO_SDET:  ");
               else if (fp.force==FRF_TM2)
                  printf("SU(3) Force FRF_TM2:          ");
               else if (fp.force==FRF_TM2_EO)
                  printf("SU(3) Force FRF_TM2_EO:       ");
               else if (fp.force==FRF_RAT)
                  printf("SU(3) Force FRF_RAT:          ");
               else if (fp.force==FRF_RAT_SDET)
                  printf("SU(3) Force FRF_RAT_SDET:     ");
               else if (fp.force==FRG_U1)
                  printf("SU(3) Force FRG_U1:         ");

               printf("nrm = %.2e, eps = % .2e, nrm*|eps| = %.2e, "
                      "time = %.2e sec\n",nrm/fabs(eps),eps,nrm,wt2-wt1);
            }
         }

         if ((gauge()&2)!=0)
         {
            update_u1mom(isym);
            nrm=norm_square_dvec(4*VOLUME,1,(*mdfs).u1frc);
            nrm=sqrt(nrm/nlk);

            if (my_rank==0)
            {
               if (fp.force==FRG_SU3)
                  printf("U(1) Force FRG_SU3:          ");
               else if (fp.force==FRF_TM1)
                  printf("U(1) Force FRF_TM1:          ");
               else if (fp.force==FRF_TM1_EO)
                  printf("U(1) Force FRF_TM1_EO:       ");
               else if (fp.force==FRF_TM1_EO_SDET)
                  printf("U(1) Force FRF_TM1_EO_SDET:  ");
               else if (fp.force==FRF_TM2)
                  printf("U(1) Force FRF_TM2:          ");
               else if (fp.force==FRF_TM2_EO)
                  printf("U(1) Force FRF_TM2_EO:       ");
               else if (fp.force==FRF_RAT)
                  printf("U(1) Force FRF_RAT:          ");
               else if (fp.force==FRF_RAT_SDET)
                  printf("U(1) Force FRF_RAT_SDET:     ");
               else if (fp.force==FRG_U1)
                  printf("U(1) Force FRG_U1:         ");

               printf("nrm = %.2e, eps = % .2e, nrm*|eps| = %.2e, "
                      "time = %.2e sec\n",nrm/fabs(eps),eps,nrm,wt2-wt1);
            }
         }
      }
      else if (iop==SU3UPDATE)
      {
         update_ud(eps);
         su3mdt+=eps;
         dt=su3mdt-mdtime();
         step_mdtime(dt);
         rtau+=dt;
      }
      else if (iop==U1UPDATE)
      {
         update_ad(eps);
         u1mdt+=eps;
         dt=u1mdt-mdtime();
         step_mdtime(dt);
         rtau+=dt;
      }
   }
}

#else

void run_mdint(void)
{
   int nop;
   int iop,status[6];
   int isym1,isym3;
   double *mu,eps;
   mdstep_t *s,*sm;
   hmc_parms_t hmc;
   force_parms_t fp;
   dirac_parms_t dp;
   double u1mdt,su3mdt,dt;

   hmc=hmc_parms();
   mu=hmc.mu;
   u1mdt=su3mdt=0.0;
   reset_chrono();
   start_dfl_upd();

   s=mdsteps(&nop);
   sm=s+nop;

   isym3=isym1=0;
   for (;s<sm;s++)
   {
      iop=(*s).iop;
      eps=(*s).eps;

      if (iop>=0)
      {
         fp=force_parms(iop);

         if (fp.force==FRG_SU3)
         {
            isym3=0;
            force0(eps);
         }
         else if (fp.force==FRG_U1)
         {
            isym1=0;
            force6(eps);
         }
         else
         {
            isym3=isym1=1;
            dfl_upd(fp.isp[0]);
            dp=qlat_parms(fp.im0);
            set_dirac_parms1(&dp);
            status[2]=0;
            status[5]=0;

            if (fp.force==FRF_TM1)
               force1(mu[fp.imu[0]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO)
               force4(mu[fp.imu[0]],fp.ipf,0,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM1_EO_SDET)
               force4(mu[fp.imu[0]],fp.ipf,1,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2)
               force2(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_TM2_EO)
               force5(mu[fp.imu[0]],mu[fp.imu[1]],fp.ipf,fp.isp[0],fp.icr[0],
                      eps,status);
            else if (fp.force==FRF_RAT)
               force3(fp.irat,fp.ipf,0,fp.isp[0],
                      eps,status);
            else if (fp.force==FRF_RAT_SDET)
               force3(fp.irat,fp.ipf,1,fp.isp[0],
                      eps,status);

            chk_mode_regen(fp.isp[0],status);
            add2counter("force",iop,status);
         }
      }
      else if (iop==SU3UPDATE)
      {
         update_su3mom(isym3);
         update_ud(eps);
         su3mdt+=eps;
         dt=su3mdt-mdtime();
         step_mdtime(dt);
         rtau+=dt;
      }
      else if (iop==U1UPDATE)
      {
         update_u1mom(isym1);
         update_ad(eps);
         u1mdt+=eps;
         dt=u1mdt-mdtime();
         step_mdtime(dt);
         rtau+=dt;
      }
      else
      {
         if((gauge()&1)!=0) update_su3mom(isym3);
         if((gauge()&2)!=0) update_u1mom(isym1);
      }
   }
}

#endif
