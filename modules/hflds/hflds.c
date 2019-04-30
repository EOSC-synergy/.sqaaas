
/*******************************************************************************
*
* File hflds.c
*
* Copyright (C) 2016,2017 Agostino Patella, Isabel Campos
*
* Based on openQCD-1.6/modules/uflds/uflds.c
*          openQCD-1.4/modules/tcharge/ftensor.c
* Copyright (C) 2005, 2006, 2009-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation and calculation of the global U(3) gauge fields.
*
* The externally accessible functions are
*
*   su3_dble *hdfld(void)
*     Computes the double-precision U(3) gauge field and returns its base
*     address. If it is not already allocated, the field is allocated first.
*     The SU(3) and U(1) gauge fields are assumed to have already the correct
*     boundary values.
*
*   su3 *hfld(void)
*     Computes the single-precision U(3) gauge field and returns its base
*     address. If it is not already allocated, the field is allocated first.
*     The SU(3) and U(1) gauge fields are assumed to have already the correct
*     boundary values.
*
* Notes:
*
* The U(3) gauge field is defined accordingly to
* 
*  W(x,mu) = [z(x,mu)]^q U(x,mu) exp(i sum_k theta[k]/Nk) s(x0,mu)
*  with  W = U(3) gauge field
*        U = SU(3) gauge field
*        z = U(1) complex gauge field
*        q = electric charge
*        theta = Angles specifying the phase-periodic boundary conditions for
*                the quark fields
*        s(x0) = (-1)^{delta(x0,N0-1) delta(mu,0)} for periodic b.c.s in time
*                1                                 otherwise
*
* The relevant parameters must be set with set_dirac_parms1() or
* set_dirac_parms9() (see flags/lat_parms.c).
*
*******************************************************************************/

#define HFLDS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lattice.h"
#include "su3.h"
#include "su3fcts.h"
#include "uflds.h"
#include "u1flds.h"
#include "flags.h"
#include "hflds.h"
#include "global.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static su3 *hb=NULL;
static su3_dble *hdb=NULL;



static void alloc_hd(void)
{
   size_t n;

   error_root(sizeof(su3_dble)!=(18*sizeof(double)),1,"alloc_hd [hflds.c]",
              "The su3_dble structures are not properly packed");

   error(iup[0][0]==0,1,"alloc_hd [hflds.c]","Geometry arrays are not set");

   n=4*VOLUME;
   hdb=amalloc(n*sizeof(*hdb),ALIGN);
   error(hdb==NULL,1,"alloc_hd [hflds.c]",
         "Unable to allocate memory space for the gauge field");
}


static void alloc_h(void)
{
   size_t n;

   error_root(sizeof(su3)!=(18*sizeof(float)),1,"alloc_h [hflds.c]",
              "The su3 structures are not properly packed");

   n=4*VOLUME;
   hb=amalloc(n*sizeof(*hb),ALIGN);
   error(hb==NULL,1,"alloc_h [hflds.c]",
         "Unable to allocate memory space for the gauge field");
}


static void mult_hd_phase(void)
{
   int k;
   double p;
   su3_dble *hd,*hm;
   dirac_parms_t dp;
   static complex_dble phase[3] ALIGNED16;
   
   dp=dirac_parms();
   
   if((dp.theta[0]!=0.0)||(dp.theta[1]!=0.0)||(dp.theta[2]!=0.0))
   {
      p=dp.theta[0]/(double)(N1);
      phase[0].re=cos(p);
      phase[0].im=sin(p);

      p=dp.theta[1]/(double)(N2);
      phase[1].re=cos(p);
      phase[1].im=sin(p);

      p=dp.theta[2]/(double)(N3);
      phase[2].re=cos(p);
      phase[2].im=sin(p);

      hd=hdb;
      hm=hdb+4*VOLUME;

      for (;hd<hm;)
      {
         hd+=2;

         for (k=0;k<3;k++)
         {
            cm3x3_mulc(phase+k,hd,hd);
            hd+=1;
            cm3x3_mulc(phase+k,hd,hd);
            hd+=1;
         }
      }
   }
}


static void chs_hd0(void)
{
   int nlks,*lks,*lkm;
   su3_dble *vd;

   lks=bnd_lks(&nlks);
   lkm=lks+nlks;

   for (;lks<lkm;lks++)
   {
      vd=hdb+(*lks);

      (*vd).c11.re=-(*vd).c11.re;
      (*vd).c11.im=-(*vd).c11.im;
      (*vd).c12.re=-(*vd).c12.re;
      (*vd).c12.im=-(*vd).c12.im;
      (*vd).c13.re=-(*vd).c13.re;
      (*vd).c13.im=-(*vd).c13.im;

      (*vd).c21.re=-(*vd).c21.re;
      (*vd).c21.im=-(*vd).c21.im;
      (*vd).c22.re=-(*vd).c22.re;
      (*vd).c22.im=-(*vd).c22.im;
      (*vd).c23.re=-(*vd).c23.re;
      (*vd).c23.im=-(*vd).c23.im;

      (*vd).c31.re=-(*vd).c31.re;
      (*vd).c31.im=-(*vd).c31.im;
      (*vd).c32.re=-(*vd).c32.re;
      (*vd).c32.im=-(*vd).c32.im;
      (*vd).c33.re=-(*vd).c33.re;
      (*vd).c33.im=-(*vd).c33.im;
   }
}


static complex_dble cpow(complex_dble base, int e)
{
   double tmp;
   complex_dble res;
   
   res.re=1.0;
   res.im=0.0;

   if(e==0) return res;

   if(e<0) 
   {
      base.im=-base.im;
      e=-e;
   }

   while(e)
   {
      if(e&1) 
      {
         /* res *= base; */
         tmp=res.re*base.re-res.im*base.im;
         res.im=res.re*base.im+res.im*base.re;
         res.re=tmp;
      }

      e>>=1;
      /* base *= base; */
      tmp=base.re*base.re-base.im*base.im;
      base.im=2.*base.re*base.im;
      base.re=tmp;
   }

   return res;
}


static void build_hd(void)
{
   su3_dble *hd,*hdm,*ud;
   complex_dble *u1d;
   complex_dble z; 
   int q,gg;
   
   q=dirac_parms().qhat;
   gg=gauge();

   hdm=hdb+4*VOLUME;
   if((gg==3)&&(q!=0))
   {
      ud=udfld();
      u1d=u1dfld(LOC);
      for(hd=hdb;hd<hdm;hd++) 
      {
         z=cpow(*u1d,q);
         cm3x3_mulc(&z,ud,hd);

         ud++;
         u1d++;
      }
   }
   else if((gg&1)!=0)
   {
      ud=udfld();
      for(hd=hdb;hd<hdm;hd++) {
         (*hd)=(*ud);
         ud++;
      }
   }
   else
   {
      u1d=u1dfld(LOC);
      for(hd=hdb;hd<hdm;hd++) 
      {
         z=cpow(*u1d,q);
         cm3x3_zero(1,hd);
         (*hd).c11=z;
         (*hd).c22=z;
         (*hd).c33=z;
         u1d++;
      }
   }

   mult_hd_phase();

   if (bc_type()==3)
      chs_hd0();
}


su3_dble *hdfld(void)
{

   if (query_flags(HD_UP2DATE)!=1)
   {
      if (hdb==NULL)
         alloc_hd();

      build_hd();
      set_flags(COMPUTED_HD);
   }
   
   error_root(query_flags(HD_UP2DATE)!=1,1,
              "hdfld [hflds.c]",
              "Flags debug: Something is not working");

   return hdb;
}


static void build_h(void)
{
   su3 *h,*hm;
   su3_dble *hd;
   
   hm=hb+4*VOLUME;
   hd=hdfld();

   for (h=hb;h<hm;h++)
   {
      (*h).c11.re=(float)((*hd).c11.re);
      (*h).c11.im=(float)((*hd).c11.im);
      (*h).c12.re=(float)((*hd).c12.re);
      (*h).c12.im=(float)((*hd).c12.im);
      (*h).c13.re=(float)((*hd).c13.re);
      (*h).c13.im=(float)((*hd).c13.im);

      (*h).c21.re=(float)((*hd).c21.re);
      (*h).c21.im=(float)((*hd).c21.im);
      (*h).c22.re=(float)((*hd).c22.re);
      (*h).c22.im=(float)((*hd).c22.im);
      (*h).c23.re=(float)((*hd).c23.re);
      (*h).c23.im=(float)((*hd).c23.im);

      (*h).c31.re=(float)((*hd).c31.re);
      (*h).c31.im=(float)((*hd).c31.im);
      (*h).c32.re=(float)((*hd).c32.re);
      (*h).c32.im=(float)((*hd).c32.im);
      (*h).c33.re=(float)((*hd).c33.re);
      (*h).c33.im=(float)((*hd).c33.im);
      
      hd+=1;
   }
}


su3 *hfld(void)
{
   if (query_flags(H_UP2DATE)!=1)
   {
      if (hb==NULL)
         alloc_h();

      build_h();
      set_flags(ASSIGNED_HD2H);
   }

   error_root(query_flags(H_UP2DATE)!=1,1,
              "hfld [hflds.c]",
              "Flags debug: Something is not working");

   return hb;
}
