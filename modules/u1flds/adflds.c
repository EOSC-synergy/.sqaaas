
/*******************************************************************************
*
* File adflds.c
*
* Copyright (C) 2016 Marina Marinkovic
*
* Based on openQCD-1.6/modules/udflds/udflds.c
* Copyright (C) 2006, 2010-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and initialization of the global non-compact U(1) gauge field.
*
* The externally accessible functions are
*
*   double *adfld(void)
*     Returns the base address of the non-compact U(1) gauge field. If it is not 
*     already allocated, the field is allocated and initialized to zero. Then
*     the boundary conditions are set according to the data base by calling
*     set_ad_bc() [ad_bcnds.c].
*
*   void random_ad(void)
*     Initializes the active non-compact U(1) link variables to random numbers,
*     uniformly distributed in [-pi,pi). The static link variables are
*     left untouched.
*
*   void renormalize_ad(void)
*     If the compact formulation is chosen, it adds a multiple of 2*pi to the
*     non-compact U(1) gauge field, in such a way that it takes values in the
*     [-pi,pi) interval. C* boundary conditions are imposed by means of the
*     orbifold constraint.
*
* Notes:
*
* The double-precision field can only be allocated after the geometry arrays
* are set up. All programs in this module act globally and must be called on
* all MPI processes simultaneously.
*
*******************************************************************************/

#define ADFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1flds.h"
#include "global.h"

#define N0 (NPROC0*L0)

static double *adb=NULL;


static void alloc_ad(void)
{
   int bc;
   size_t n;
   double *a,*am;


   error(iup[0][0]==0,1,"alloc_ad [adflds.c]","Geometry arrays are not set");

   bc=bc_type();
   n=4*VOLUME+7*(BNDRY/4);

   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
      n+=3;
   
   adb=amalloc(n*sizeof(*adb),ALIGN);
   error(adb==NULL,1,"alloc_ad [adflds.c]",
         "Unable to allocate memory space for the gauge field");

   a=adb;
   am=adb+n;
   for (;a<am;a++)
      (*a)=0.0; 

   set_flags(UPDATED_AD);
   set_ad_bc(); 
}


double *adfld(void)
{
   if (adb==NULL)
   {
      alloc_ad();
   }

   return adb;
}


void random_ad(void)
{
   double *a,*am;
   double twopi;

   a=adfld();
   twopi=8.0*atan(1.0);

   ranlxd(a,4*VOLUME);
   am=a+4*VOLUME;
   for (;a<am;a++)
   {
      (*a)=((*a)-0.5)*twopi;
   }

   set_flags(UPDATED_AD);
   set_ad_bc();
   orbi_cpy_ad();
}


void renormalize_ad(void)
{
   int bc,ix,t,ifc;
   double *ad;
   double twopi;
   
   if(u1lat_parms().type==0)
   {
      twopi=8.0*atan(1.0);

      bc=bc_type();
      ad=adfld();

      for (ix=(VOLUME/2);ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if (t==0)
         {
            (*ad)-=twopi*(floor((*ad)/twopi+0.5));
            ad+=1;

            if (bc!=0)
               (*ad)-=twopi*(floor((*ad)/twopi+0.5));
            ad+=1;

            for (ifc=2;ifc<8;ifc++)
            {
               if (bc!=1)
                  (*ad)-=twopi*(floor((*ad)/twopi+0.5));
               ad+=1;
            }
         }
         else if (t==(N0-1))
         {
            if (bc!=0)
               (*ad)-=twopi*(floor((*ad)/twopi+0.5));
            ad+=1;

            for (ifc=1;ifc<8;ifc++)
            {
               (*ad)-=twopi*(floor((*ad)/twopi+0.5));
               ad+=1;
            }
         }
         else
         {
            for (ifc=0;ifc<8;ifc++)
            {
               (*ad)-=twopi*(floor((*ad)/twopi+0.5));
               ad+=1;
            }
         }
      }
      set_flags(UPDATED_AD);
   }
   
   orbi_cpy_ad();
}
