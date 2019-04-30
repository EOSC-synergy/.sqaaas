
/*******************************************************************************
*
* File wflow_u1.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/modules/wflow/wflow.c
* Copyright (C) 2009-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Integration of the Wilson flow for U(1) gauge field.
*
* The externally accessible functions are
*
*   void fwd_u1_euler(int n,double eps)
*     Applies n forward Euler integration steps, with step size eps, to the
*     current U(1) gauge field.
*
*   void fwd_u1_rk2(int n,double eps)
*     Applies n forward 2nd-order Runge-Kutta integration steps, with step
*     size eps, to the current U(1) gauge field.
*
*   void fwd_u1_rk3(int n,double eps)
*     Applies n forward 3rd-order Runge-Kutta integration steps, with step
*     size eps, to the current U(1) gauge field.
*
*******************************************************************************/

#define WFLOW_C

#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1flds.h"
#include "mdflds.h"
#include "linalg.h"
#include "forces.h"
#include "wflow.h"
#include "global.h"

#define N0 (NPROC0*L0)


static void update_ad(double eps,double *frc)
{
   int bc,ix,t,ifc;
   double *a;

   bc=bc_type();
   a=adfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         (*a)+=eps*(*frc);
         frc+=1;
         a+=1;

         if (bc!=0)
            (*a)+=eps*(*frc);
         frc+=1;
         a+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if ((bc==0)||(bc==2))
               (*a)+=2.0*eps*(*frc);
            else if (bc==3)
               (*a)+=eps*(*frc);
            frc+=1;
            a+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            (*a)+=eps*(*frc);
         frc+=1;
         a+=1;

         (*a)+=eps*(*frc);
         frc+=1;
         a+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc==0)
               (*a)+=2.0*eps*(*frc);
            else
               (*a)+=eps*(*frc);
            frc+=1;
            a+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            (*a)+=eps*(*frc);
            frc+=1;
            a+=1;
         }
      }
   }

   set_flags(UPDATED_AD);
}


static void update_fro1(double c,double *frc,double *fro)
{
   double *frm;

   frm=frc+4*VOLUME;

   for (;frc<frm;frc++)
   {
      (*fro)-=c*(*frc);

      fro+=1;
   }
}


static void update_fro2(double c,double *frc,double *fro)
{
   double *frm;

   frm=frc+4*VOLUME;

   for (;frc<frm;frc++)
   {
      (*fro)=(*frc)+c*(*fro);

      fro+=1;
   }
}


void fwd_u1_euler(int n,double eps)
{
   int iprms[1],k;
   double dprms[1];
   double *frc;
   mdflds_t *mdfs;

   if (NPROC>1)
   {
      iprms[0]=n;
      dprms[0]=eps;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=n)||(dprms[0]!=eps),1,
            "fwd_u1_euler [wflow_u1.c]","Parameters are not global");
   }

   if (n>0)
   {
      mdfs=mdflds();
      frc=(*mdfs).u1frc;

      for (k=0;k<n;k++)
      {
         plaq_u1frc();
         update_ad(-eps,frc);
      }
   }
}


void fwd_u1_rk2(int n,double eps)
{
   int iprms[1],k;
   double dprms[1];
   double *frc,*fro,**fsv;
   mdflds_t *mdfs;

   if (NPROC>1)
   {
      iprms[0]=n;
      dprms[0]=eps;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=n)||(dprms[0]!=eps),1,
            "fwd_u1_rk2 [wflow_u1.c]","Parameters are not global");
   }

   if (n>0)
   {
      mdfs=mdflds();
      frc=(*mdfs).u1frc;
      fsv=reserve_wf1d(1);
      fro=fsv[0];

      for (k=0;k<n;k++)
      {
         plaq_u1frc();
         assign_dvec2dvec(4*VOLUME,frc,fro);
         update_ad(-0.5*eps,frc);

         plaq_u1frc();
         update_fro2(-0.5,frc,fro);
         update_ad(-eps,fro);
      }

      release_wf1d();
   }
}


void fwd_u1_rk3(int n,double eps)
{
   int iprms[1],k;
   double dprms[1];
   double *frc,*fro,**fsv;
   mdflds_t *mdfs;

   if (NPROC>1)
   {
      iprms[0]=n;
      dprms[0]=eps;

      MPI_Bcast(iprms,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      error((iprms[0]!=n)||(dprms[0]!=eps),1,
            "fwd_u1_rk3 [wflow_u1.c]","Parameters are not global");
   }

   if (n>0)
   {
      mdfs=mdflds();
      frc=(*mdfs).u1frc;
      fsv=reserve_wf1d(1);
      fro=fsv[0];

      for (k=0;k<n;k++)
      {
         plaq_u1frc();
         assign_dvec2dvec(4*VOLUME,frc,fro);
         update_ad(-0.25*eps,frc);

         plaq_u1frc();
         update_fro1(32.0/17.0,frc,fro);
         update_ad((17.0/36.0)*eps,fro);

         plaq_u1frc();
         update_fro2(17.0/27.0,frc,fro);
         update_ad(-0.75*eps,fro);
      }

      release_wf1d();
   }
}
