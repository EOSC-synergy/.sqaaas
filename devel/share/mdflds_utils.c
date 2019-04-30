
/*******************************************************************************
*
* File mdflds_utils.c
*
* Copyright (C) 2016 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization, translation and guage transformations of the gauge fields.
*
*******************************************************************************/

#define MDFLDS_UTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "flags.h"
#include "uflds.h"
#include "global.h"
#include "u1flds.h"
#include "mdflds.h"
#include "mdflds_utils.h"

#define FILENAME "mdflds_utils.c"

#define N0 (NPROC0*L0)


void rot_ud(double eps)
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
         mom+=1;
         u+=1;

         if (bc!=0)
            expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,mom,u);
         mom+=1;
         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            mom+=1;
            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


void rot_ad(double eps)
{
   int bc,ix,t,ifc;
   double *u;
   double *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).u1mom;
   u=adfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         (*u)+=eps*(*mom);
         mom+=1;
         u+=1;

         if (bc!=0)
            (*u)+=eps*(*mom);
         mom+=1;
         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               (*u)+=eps*(*mom);
            mom+=1;
            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            (*u)+=eps*(*mom);
         mom+=1;
         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            (*u)+=eps*(*mom);
            mom+=1;
            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            (*u)+=eps*(*mom);
            mom+=1;
            u+=1;
         }
      }
   }

   set_flags(UPDATED_AD);
}


static int is_su3frc_zero(su3_alg_dble *f)
{
   int ie;

   ie=1;
   ie&=((*f).c1==0.0);
   ie&=((*f).c2==0.0);
   ie&=((*f).c3==0.0);
   ie&=((*f).c4==0.0);
   ie&=((*f).c5==0.0);
   ie&=((*f).c6==0.0);
   ie&=((*f).c7==0.0);
   ie&=((*f).c8==0.0);

   return ie;
}


int check_bnd_su3frc(su3_alg_dble *frc)
{
   int bc,ix,t,ifc,ie;

   bc=bc_type();
   ie=0;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t==0)&&(bc==0))
      {
         ie|=is_su3frc_zero(frc);
         frc+=1;

         ie|=(is_su3frc_zero(frc)^0x1);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=is_su3frc_zero(frc);
            frc+=1;
         }
      }
      else if ((t==0)&&(bc==1))
      {
         ie|=is_su3frc_zero(frc);
         frc+=1;

         ie|=is_su3frc_zero(frc);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=(is_su3frc_zero(frc)^0x1);
            frc+=1;
         }
      }
      else if ((t==(N0-1))&&(bc==0))
      {
         ie|=(is_su3frc_zero(frc)^0x1);
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            ie|=is_su3frc_zero(frc);
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            ie|=is_su3frc_zero(frc);
            frc+=1;
         }
      }
   }

   return ie;
}


int check_zero_su3frc(su3_alg_dble *frc)
{
   su3_alg_dble *fm;
   int ie;

   ie=0;
   fm=frc+4*VOLUME;
   for (;frc<fm;frc++)
      ie|=(is_su3frc_zero(frc)^0x1);
   
   return ie;
}




int check_bnd_u1frc(double *frc)
{
   int bc,ix,t,ifc,ie;

   bc=bc_type();
   ie=0;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t==0)&&(bc==0))
      {
         ie|=((*frc)==0.0);
         frc+=1;

         ie|=((*frc)!=0.0);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=((*frc)==0.0);
            frc+=1;
         }
      }
      else if ((t==0)&&(bc==1))
      {
         ie|=((*frc)==0.0);
         frc+=1;

         ie|=((*frc)==0.0);
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            ie|=((*frc)!=0.0);
            frc+=1;
         }
      }
      else if ((t==(N0-1))&&(bc==0))
      {
         ie|=((*frc)!=0.0);
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            ie|=((*frc)==0.0);
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            ie|=((*frc)==0.0);
            frc+=1;
         }
      }
   }

   return ie;
}


int check_zero_u1frc(double *frc)
{
   double *fm;
   int ie;

   ie=0;
   fm=frc+4*VOLUME;
   for (;frc<fm;frc++)
      ie|=((*frc)!=0.0);

   return ie;
}
