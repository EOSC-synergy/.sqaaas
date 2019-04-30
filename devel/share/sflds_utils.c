
/*******************************************************************************
*
* File sflds_utils.c
*
* Copyright (C) 2016 Agostino Patella
*
* Based on openQCD-1.4/devel/dirac/check1.c
* Copyright (C) 2005, 2011-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge transformations for the spinor fields.
*
*******************************************************************************/

#define SFLDS_UTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "global.h"
#include "gflds_utils.h"
#include "sflds_utils.h"

#define FILENAME "sflds_utils.c"



void transform_sd(spinor_dble *pk,spinor_dble *pl)
{
   int ix,i,q;
   su3_dble g3x;
   complex_dble zx,sc,*rc;
   spinor_dble r,s;
   su3_dble *g3;
   double *g1;

   if((gauge()&1)!=0)
   {
      g3=g3tr();
      
      for (ix=0;ix<VOLUME;ix++)
      {
         s=pk[ix];
         g3x=g3[ix];

         _su3_multiply(r.c1,g3x,s.c1);
         _su3_multiply(r.c2,g3x,s.c2);
         _su3_multiply(r.c3,g3x,s.c3);
         _su3_multiply(r.c4,g3x,s.c4);

         pl[ix]=r;
      }
   }
   else
   {
      for (ix=0;ix<VOLUME;ix++)
      {
         pl[ix]=pk[ix];
      }
   }
   
   if((gauge()&2)!=0)
   {
      g1=g1tr();
      q=dirac_parms().qhat;
      
      for (ix=0;ix<VOLUME;ix++)
      {
         zx.re=cos(q*g1[ix]);
         zx.im=sin(q*g1[ix]);
         
         #define cmpl_mul(a,b,c) \
            (a).re=(b).re*(c).re-(b).im*(c).im; \
            (a).im=(b).re*(c).im+(b).im*(c).re
         
         for(i=0;i<12;i++)
         {
            rc=(complex_dble*)(pl+ix)+i;
            sc=(*rc);
            cmpl_mul((*rc),zx,sc);
         }
         
         #undef cmpl_mul
      }
   }
}
