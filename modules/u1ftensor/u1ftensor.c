
/*******************************************************************************
*
* File u1ftensor.c
*
* Copyright (C) 2016, 2017 Nazario Tantalo
*               2017 Agostino Patella
*
* Based on openQCD-1.4/modules/tcharge/ftensor.c
* Copyright (C) 2010-2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the symmetric U(1) field tensor.
*
* The externally accessible function is
*
*   double **u1ftensor(void)
*     Computes the symmetric field tensor of the global double-precision
*     U(1) gauge field and returns the pointers ft[0],..,ft[5] to the field
*     components with the Lorentz indices (0 1),(0 2),(0 3),(2 3),(3 1),
*     (1 2). The arrays are automatically allocated if needed. Along the
*     boundaries of the lattice (if any), the program sets the field to
*     zero.
*
* Notes:
*
* At all points x in the interior of the lattice, the (mu nu)-component of
* the field tensor is defined by
*
*  Ahat_{mu nu}(x) = (1/4*qel)*Im{ z_{mu nu}(x)    + z_{mu nu}(x-mu) + 
*                                  z_{mu nu}(x-nu) + z_{mu nu}(x-mu-nu) },
*
* where
*
*  z_{mu nu}(x) = z(x,mu)*z(x+mu,nu)*z(x+nu,mu)^star*z(x,nu)^star
*
* denotes the plaquette loop at x in the (mu nu)-plane and qel is the elementary
* charge. Elsewhere the elements of the field arrays are set to zero. The 
* interior points are  those at global time x0 in the range
*
*  0<x0<NPROC0*L0-1        (open bc),
*
*  0<x0<NPROC0*L0          (SF and open-SF bc),
*
*  0<=x0<NPROC0*L0         (periodic bc).
*
* The program in this module performs global operations and must be called
* simultaneously on all MPI processes.
*
*******************************************************************************/

#define U1FTENSOR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1flds.h"
#include "linalg.h"
#include "u1ftensor.h"
#include "global.h"

#define N0 (NPROC0*L0)

static double **fts=NULL,**ft;
static ftidx_t *idx;


static void alloc_u1fts(void)
{
   int n,nbf;
   double **pp,*p;

   idx=ftidx();
   nbf=0;

   for (n=0;n<6;n++)
      nbf+=idx[n].nft[0]+idx[n].nft[1];

   pp=malloc(12*sizeof(*pp));
   p=amalloc((6*VOLUME+nbf)*sizeof(*p),ALIGN);
   error((pp==NULL)||(p==NULL),1,"alloc_u1fts [u1ftensor.c]",
	 "Unable to allocate field tensor arrays");

   fts=pp;
   ft=pp+6;

   for (n=0;n<6;n++)
   {
      (*pp)=p;
      pp+=1;
      p+=VOLUME+idx[n].nft[0]+idx[n].nft[1];
   }
}


static void build_u1fts(void)
{
   int bc,n,ix,t,ip[4],ipf[4],tmx;
   double *ftn,zmunu,iqel;
   complex_dble *u1b,w1,w2;
   u1lat_parms_t lat;

   lat=u1lat_parms();
   iqel=lat.invqel;

   bc=bc_type();
   u1b=u1dfld(EXT);

   for (n=0;n<6;n++)
   {
      ftn=fts[n];
      set_dvec2zero(VOLUME+idx[n].nft[0]+idx[n].nft[1],ftn);
      tmx=N0;
      if (bc==0)
         tmx-=1;
      if (n<3)
         tmx-=1;

      for (ix=0;ix<VOLUME;ix++)
      {
         t=global_time(ix);
         plaq_uidx(n,ix,ip);
         plaq_ftidx(n,ix,ipf);

         u1xu1(u1b+ip[0],u1b+ip[1],&w1);
         u1xu1(u1b+ip[2],u1b+ip[3],&w2);

         zmunu= 0.25*iqel*(w1.im*w2.re-w1.re*w2.im);

         if (((t>0)&&(t<tmx))||(bc==3))
         {
            ftn[ipf[0]]+=zmunu;
            ftn[ipf[1]]+=zmunu;
            ftn[ipf[2]]+=zmunu;
            ftn[ipf[3]]+=zmunu;
         }
         else if ((t==0)&&(n<3))
         {
            ftn[ipf[1]]+=zmunu;
            ftn[ipf[3]]+=zmunu;
         }
         else if ((t==tmx)&&(n<3))
         {
            ftn[ipf[0]]+=zmunu;
            ftn[ipf[2]]+=zmunu;
         }
      }

      add_bnd_u1ft(n,ftn);
   }
}

double **u1ftensor(void)
{
   int n;

   if (query_flags(U1FTS_UP2DATE)!=1)
   {
      if (fts==NULL)
         alloc_u1fts();

      build_u1fts();
      set_flags(COMPUTED_U1FTS);
   }

   for (n=0;n<6;n++)
      ft[n]=fts[n];

   error_root(query_flags(U1FTS_UP2DATE)!=1,1,
              "u1ftensor [u1ftensor.c]",
              "Flags debug: Something is not working");

   return ft;
}

#undef N0
