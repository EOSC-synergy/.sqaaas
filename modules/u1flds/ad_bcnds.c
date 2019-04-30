
/*******************************************************************************
*
* File ad_bcnds.c
*
* Copyright (C) 2016 Marina Marinkovic
*
* Based on openQCD-1.6/modules/lattice/bcnds.c
* Copyright (C) 2005, 2010-2014 Martin Luescher, John Bulava
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs related to the boundary conditions in the time direction for the 
* non-compact U(1) gauge field.
*
* The externally accessible functions are
*
*   void set_ad_bc(void)
*     Sets the non-compact U(1) link variables at time 0 and T to the
*     values required by the chosen boundary conditions (see the notes).
*
*   int check_ad_bc(double tol)
*     Returns 1 if the non-compact U(1) gauge field has the proper boundary
*     values. Otherwise the program returns 0. The parameter tol>=0.0 sets an
*     upper bound on the tolerated difference of the boundary values of the
*     gauge field from the expected ones in the case of SF and open-SF
*     boundary conditions.
*
*
* Notes:
*
* The time extent T of the lattice is
*
*  NPROC0*L0-1      for open boundary conditions,
*
*  NPROC0*L0        for SF, open-SF and periodic boundary conditions.
*
* Note that in the latter cases the points at time T are not in the local
* lattice.
*
* The action performed by set_bc() is the following:
*
*  Open bc:         Set all link variables A(x,0) at time T to zero.
*
*  SF bc:           Set all link variables A(x,k) for k=1,2,3 at time 0 and T
*                   to the values specified in the parameter database.
*
*  Open-SF bc:      Set all link variables A(x,k) for k=1,2,3 at time T to zero
*                   to the values specified in the parameter database.
*
*  Periodic bc:     No action is performed.
*
* The programs in this module act globally and must be called simultaneously
* on all MPI processes.
*
*******************************************************************************/

#define AD_BCNDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "flags.h"
#include "mpi.h"
#include "random.h"
#include "u1flds.h"
#include "lattice.h"
#include "global.h"


#define N0 (NPROC0*L0)

static int init0=0,nlks,*lks;
static int init1=0,npts,*pts;
static int init2=0;
static double abnd[2][3];


static void alloc_lks(void)
{
   int ix,t,*lk;

   error(iup[0][0]==0,1,"alloc_lks [ad_bcnds.c]","Geometry arrays are not set");

   if ((cpr[0]==0)||(cpr[0]==(NPROC0-1)))
   {
      if (NPROC0>1)
         nlks=(L1*L2*L3)/2;
      else
         nlks=L1*L2*L3;

      lks=malloc(nlks*sizeof(*lks));

      if (lks!=NULL)
      {
         lk=lks;

         for (ix=(VOLUME/2);ix<VOLUME;ix++)
         {
            t=global_time(ix);

            if (t==0)
            {
               (*lk)=8*(ix-(VOLUME/2))+1;
               lk+=1;
            }
            else if (t==(N0-1))
            {
               (*lk)=8*(ix-(VOLUME/2));
               lk+=1;
            }
         }
      }
   }
   else
   {
      nlks=0;
      lks=NULL;
   }

   error((nlks>0)&&(lks==NULL),1,"alloc_lks [ad_bcnds.c]",
         "Unable to allocate index array");
   init0=1;
}


static void alloc_pts(void)
{
   int bc,ix,t,*pt;

   error(iup[0][0]==0,1,"alloc_pts [ad_bcnds.c]","Geometry arrays are not set");
   bc=bc_type();

   if (((cpr[0]==0)&&(bc!=3))||((cpr[0]==(NPROC0-1))&&(bc==0)))
   {
      if ((NPROC0==1)&&(bc==0))
         npts=2*L1*L2*L3;
      else
         npts=L1*L2*L3;

      pts=malloc(npts*sizeof(*pts));

      if (pts!=NULL)
      {
         pt=pts;

         for (ix=0;ix<VOLUME;ix++)
         {
            t=global_time(ix);

            if ((t==0)||((t==(N0-1))&&(bc==0)))
            {
               (*pt)=ix;
               pt+=1;
            }
         }
      }
   }
   else
   {
      npts=0;
      pts=NULL;
   }

   error((npts>0)&&(pts==NULL),1,"alloc_pts [ad_bcnds.c]",
         "Unable to allocate index array");
   init1=1;
}


static int ad_is_equal(double tol,double *u,double *v)
{
   return (fabs((*u)-(*v))<=tol);
}


static int ad_check_zero(int bc)
{
   int it,ix,t;
   double *u;

   it=1;
   u=adfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((bc==0)&&(t==0))
      {
         u+=1;
         it&=((*u)==0.0);
         u+=1;
      }
      else if ((bc==0)&&(t==(N0-1)))
      {
         it&=((*u)==0.0);
         u+=2;
      }
      else
      {
         u+=2;
      }
      u+=6;
   }

   return it;
}


static void set_abnd(void)
{
   int i,k;
   double s[3];
   bc_parms_t bcp;

   bcp=bc_parms();
   s[0]=(double)(NPROC1*L1);
   s[1]=(double)(NPROC2*L2);
   s[2]=(double)(NPROC3*L3);

   for (i=0;i<2;i++)
   {
      for (k=0;k<3;k++)
      {
         abnd[i][k]=bcp.phi1[i]/s[k];
         abnd[i][k]=bcp.phi1[i]/s[k];
      }
   }

   init2=1;
}


static void open_bc(void)
{
   int *lk,*lkm;
   double *ub;

   if (init0==0)
      alloc_lks();

   ub=adfld();
   lk=lks;
   lkm=lk+nlks;

   for (;lk<lkm;lk++)
      ub[*lk]=0.0;

   set_flags(UPDATED_AD);
}


static void SF_bc(void)
{
   int k,*pt,*ptm;
   double *ub,*u;

   if (init1==0)
      alloc_pts();
   if (init2==0)
      set_abnd();

   ub=adfld();

   if (cpr[0]==0)
   {
      pt=pts+(npts/2);
      ptm=pts+npts;

      for (;pt<ptm;pt++)
      {
         u=ub+8*(pt[0]-(VOLUME/2));

         for (k=0;k<3;k++)
         {
            u[2+2*k]=abnd[0][k];
            u[3+2*k]=abnd[0][k];
         }
      }
   }

   if (cpr[0]==(NPROC0-1))
   {
      u=ub+4*VOLUME+7*(BNDRY/4);

      for (k=0;k<3;k++)
         u[k]=abnd[1][k];
   }

   set_flags(UPDATED_AD);
}


static void openSF_bc(void)
{
   int k;
   double *ub,*u;

   if (init2==0)
      set_abnd();

   ub=adfld();

   if (cpr[0]==(NPROC0-1))
   {
      u=ub+4*VOLUME+7*(BNDRY/4);

      for (k=0;k<3;k++)
         u[k]=abnd[1][k];
   }

   set_flags(UPDATED_AD);
}


void set_ad_bc(void)
{
   int bc,it;

   bc=bc_type();

   if (bc==0)
      open_bc();
   else if (bc==1)
      SF_bc();
   else if (bc==2)
      openSF_bc();

   it=ad_check_zero(bc);
   error(it!=1,1,"set_bc [ad_bcnds.c]",
         "Non-compact U(1) link variables vanish on an incorrect set of links");
}


static int check_SF(double tol)
{
   int it,k,*pt,*ptm;
   double *ub,*u;

   if (init1==0)
      alloc_pts();
   if (init2==0)
      set_abnd();

   it=1;
   ub=adfld();

   if (cpr[0]==0)
   {
      pt=pts+(npts/2);
      ptm=pts+npts;

      for (;pt<ptm;pt++)
      {
         u=ub+8*(pt[0]-(VOLUME/2));

         for (k=0;k<3;k++)
         {
            it&=ad_is_equal(tol,u+2+2*k,abnd[0]+k);
            it&=ad_is_equal(tol,u+3+2*k,abnd[0]+k);
         }
      }
   }

   if (cpr[0]==(NPROC0-1))
   {
      u=ub+4*VOLUME+7*(BNDRY/4);

      for (k=0;k<3;k++)
         it&=ad_is_equal(tol,u+k,abnd[1]+k);
   }

   return it;
}


static int check_openSF(double tol)
{
   int it,k;
   double *ub,*u;

   if (init2==0)
      set_abnd();

   it=1;
   ub=adfld();

   if (cpr[0]==(NPROC0-1))
   {
      u=ub+4*VOLUME+7*(BNDRY/4);

      for (k=0;k<3;k++)
         it&=ad_is_equal(tol,u+k,abnd[1]+k);
   }

   return it;
}

int check_ad_bc(double tol)
{
   int bc,it,is;
   double dprms[1];

   if (NPROC>1)
   {
      dprms[0]=tol;
      MPI_Bcast(dprms,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      error(dprms[0]!=tol,1,"ad_check_bc [ad_bcnds.c]","Parameter is not global");
   }

   bc=bc_type();
   it=ad_check_zero(bc);

   if (bc==1)
      it&=check_SF(tol);
   else if (bc==2)
      it&=check_openSF(tol);

   if (NPROC>1)
   {
      is=it;
      MPI_Allreduce(&is,&it,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
   }

   return it;
}
