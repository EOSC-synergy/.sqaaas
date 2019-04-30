
/*******************************************************************************
*
* File u1flds.c
*
* Copyright (C) 2015,2016 Marina Marinkovic, Nazario Tantalo, Agostino Patella
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and calculation of the global compact U(1) gauge field.
*
* The externally accessible functions are
*
*   complex_dble *u1dfld(u1lat_t lat)
*     Returns the base address of the compact U(1) gauge field. If it is not
*     allocated yet, the field is allocated and calculated by exponentiating
*     the non-compact U(1) field. If lat==LOC the compact field is calculated
*     only on the local lattice, if lat==EXT it is calculated also on the
*     communication buffers. The compact U(1) link variable are set to be zero
*     on the inactive links in case of open boundary conditions.
*
* Notes:
*
* The double-precision compact U(1) field can only be allocated after the 
* geometry arrays are set up. All programs in this module act globally and must 
* be called on all MPI processes simultaneously.
*
*******************************************************************************/

#define U1FLDS_C

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

static complex_dble *u1db=NULL;
static size_t size=0;
static int nlks=0,*lks=NULL;



static void alloc_lks(void)
{
   int ix,t,*lk,ip[4];
      
   error(iup[0][0]==0,1,"alloc_lks [u1flds.c]","Geometry arrays are not set");
   
   uidx();
   
   if (cpr[0]==(NPROC0-1))
   {
      nlks=L1*L2*L3;
      if (NPROC1>1)
      {
         nlks+=L2*L3;
      }
      if (NPROC2>1)
      {
         nlks+=L1*L3;
      }
      if (NPROC3>1)
      {
         nlks+=L2*L1;
      }
      lks=malloc(nlks*sizeof(*lks));
      if (lks!=NULL)
      {
         lk=lks;
         for (ix=0;ix<VOLUME;ix++)
         {
            t=global_time(ix);
            if (t==(N0-1))
            {
               plaq_uidx(0,ix,ip);
               (*lk)=ip[0];
               lk+=1;
               if (iup[ix][1]>=VOLUME)
               {
                  (*lk)=ip[3];
                  lk+=1;
               }
               if (iup[ix][2]>=VOLUME)
               {
                  plaq_uidx(1,ix,ip);
                  (*lk)=ip[3];
                  lk+=1;
               }
               if (iup[ix][3]>=VOLUME)
               {
                  plaq_uidx(2,ix,ip);
                  (*lk)=ip[3];
                  lk+=1;
               }
            }
         }
      }
   }
   else if ((cpr[0]==0)&&(NPROC0>1))
   {
      nlks=(L1*L2*L3)/2;
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
         }
      }
   }
   else
   {
      nlks=0;
      lks=NULL;
   }

   error((nlks>0)&&(lks==NULL),1,"alloc_lks [u1flds.c]",
         "Unable to allocate index array");
}



static void expXu1(double *X,complex_dble *u)
{
   (*u).re=cos(*X);
   (*u).im=sin(*X);
}


static void alloc_u1d(void)
{
   int bc;

   error_root(sizeof(complex_dble)!=(2*sizeof(double)),1,"alloc_u1d [u1flds.c]",
              "The complex_dble structures are not properly packed");

   error(iup[0][0]==0,1,"alloc_u1d [u1flds.c]","Geometry arrays are not set");

   bc=bc_type(); 
   size=4*VOLUME+7*(BNDRY/4);

   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
      size+=3;
   
   u1db=amalloc(size*sizeof(*u1db),ALIGN);
   error(u1db==NULL,1,"alloc_u1d [u1flds.c]",
         "Unable to allocate memory space for the gauge field");

   if(bc==0)
      alloc_lks();
}


complex_dble *u1dfld(u1lat_t lat)
{
   complex_dble *u1,*u1m;
   double *a;
   size_t n;
   int *lk,*lkm;

   a=adfld();
   if (u1db==NULL)
   {
      alloc_u1d();
   }   

   if ((lat==LOC)&&(query_flags(U1D_UP2DATE)==1))
      return u1db;

   if ((lat==EXT)&&(query_flags(U1D_UP2DATE)==1)&&(query_flags(ADBUF_UP2DATE)==1)&&(query_flags(U1DBUF_UP2DATE)==1))
      return u1db;

   n=4*VOLUME;
   if (lat==EXT)
   {
      if (query_flags(ADBUF_UP2DATE)!=1)
         copy_bnd_ad();
  	   n=size;
	}

   u1m=u1db+n;
   for (u1=u1db;u1<u1m;u1++)
   {
   	expXu1(a,u1);
      a++;
   }

   if (nlks>0)
   {
      lk=lks;
      lkm=lk+nlks;
      for (;lk<lkm;lk++)
      {
         u1db[*lk].re=0.0;
         u1db[*lk].im=0.0;
      }
   }
   
   set_flags(COMPUTED_U1D);
   if (lat==EXT) set_flags(COMPUTED_BND_U1D);

   error_root((query_flags(U1D_UP2DATE)!=1)||
              ((lat==EXT)&&(query_flags(ADBUF_UP2DATE)!=1))||
              ((lat==EXT)&&(query_flags(U1DBUF_UP2DATE)!=1)),1,
              "u1fld [u1flds.c]",
              "Flags debug: Something is not working");

   return u1db;
}
