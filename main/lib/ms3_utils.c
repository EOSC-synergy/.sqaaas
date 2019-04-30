
/*******************************************************************************
*
* File ms3_utils.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/main/qcd1.c
*      and openQCD-1.6/main/ym1.c
* Copyright (C) 2011-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#include <stdio.h>
#include "flags.h"
#include "utils.h"
#include "global.h"
#include "uflds.h"
#include "wflow.h"
#include "tcharge.h"
#include "su3fcts.h"
#include "main_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static struct
{
   int nt;
   double **Wsl,**Ysl,**Qsl;
} data3;

static double *Wact,*Yact,*Qtop;

static wflow_parms_t wfp;


static void init_ms3dat(void)
{
   int nn,tmax;
   int in;
   double **pp,*p;
   static int init=0;
   
   if(init) return;

   error_root((gauge()&1)==0,1,"init_ms3dat [ms3_utils.c]",
              "SU(3) gauge group is not active");
   
   wfp=wflow_parms();

   nn=wfp.nn;
   tmax=wfp.tmax;

   pp=amalloc(3*(nn+1)*sizeof(*pp),3);
   p=amalloc(3*(nn+1)*(tmax+1)*sizeof(*p),4);

   error_root((pp==NULL)||(p==NULL),1,"init_ms3dat [ms3_utils.c]",
              "Unable to allocate data arrays");

   data3.Wsl=pp;
   data3.Ysl=pp+nn+1;
   data3.Qsl=pp+2*(nn+1);

   for (in=0;in<(3*(nn+1));in++)
   {
      *pp=p;
      pp+=1;
      p+=tmax;
   }

   Wact=p;
   p+=nn+1;
   Yact=p;
   p+=nn+1;
   Qtop=p;
   
   init=1;
}


void write_ms3dat(FILE *fdat)
{
   int iw,nn,tmax;
   int in;

   init_ms3dat();

   iw=write_little_int(0,fdat,1,data3.nt);

   nn=wfp.nn;
   tmax=wfp.tmax;

   for (in=0;in<=nn;in++)
      iw+=write_little_dblearray(0,fdat,tmax,data3.Wsl[in]);

   for (in=0;in<=nn;in++)
      iw+=write_little_dblearray(0,fdat,tmax,data3.Ysl[in]);

   for (in=0;in<=nn;in++)
      iw+=write_little_dblearray(0,fdat,tmax,data3.Qsl[in]);

   error_root(iw!=(1+3*(nn+1)*tmax),1,"write_ms3dat [ms3_utils.c]",
              "Incorrect write count");
}


int read_ms3dat(FILE *fdat,int *nt)
{
   int ir,nn,tmax;
   int in;

   init_ms3dat();

   ir=read_little_int(0,fdat,1,&(data3.nt));
   if (ir!=1) return 0;

   nn=wfp.nn;
   tmax=wfp.tmax;

   for (in=0;in<=nn;in++)
      ir+=read_little_dblearray(0,fdat,tmax,data3.Wsl[in]);

   for (in=0;in<=nn;in++)
      ir+=read_little_dblearray(0,fdat,tmax,data3.Ysl[in]);

   for (in=0;in<=nn;in++)
      ir+=read_little_dblearray(0,fdat,tmax,data3.Qsl[in]);

   error_root(ir!=(1+3*(nn+1)*tmax),1,"read_ms3dat [ms3_utils.c]",
              "Read error or incomplete data record");

   (*nt)=data3.nt;

   return 1;
}


void write_ms3dat_head(FILE *fdat)
{
   int iw;

   init_ms3dat();

   iw=write_little_int(0,fdat,3,wfp.dn,wfp.nn,wfp.tmax);
   iw+=write_little_dble(0,fdat,1,wfp.eps);

   error_root(iw!=4,1,"write_ms3dat_head [ms3_utils.c]",
              "Incorrect write count");
}


void check_ms3dat_head(FILE *fdat)
{
   init_ms3dat();
   
   check_little_int("check_ms3dat_head [ms3_utils.c]",fdat,3,wfp.dn,wfp.nn,wfp.tmax);
   check_little_dble("check_ms3dat_head [ms3_utils.c]",fdat,1,wfp.eps);
}


void set_ms3dat(int nt)
{
   int flint,in,dn,nn,x0;
   double eps;
   su3_dble **usv;

   init_ms3dat();

   usv=reserve_wud(1);
   cm3x3_assign(4*VOLUME,udfld(),usv[0]);

   data3.nt=nt;
   dn=wfp.dn;
   nn=wfp.nn;
   eps=wfp.eps;
   flint=wfp.flint;

   for (in=0;in<nn;in++)
   {
      Wact[in]=plaq_action_slices(data3.Wsl[in]);
      Yact[in]=ym_action_slices(data3.Ysl[in]);
      Qtop[in]=tcharge_slices(data3.Qsl[in]);
      if (bc_cstar()!=0)
      {
         for(x0=0;x0<N0;x0++)
         {
            data3.Wsl[in][x0]*=0.5;
            data3.Ysl[in][x0]*=0.5;
            data3.Qsl[in][x0]*=0.5;
         }
         Wact[in]*=0.5;
         Yact[in]*=0.5;
         Qtop[in]*=0.5;
      }

      if (flint==0)
         fwd_su3_euler(dn,eps);
      else if (flint==1)
         fwd_su3_rk2(dn,eps);
      else
         fwd_su3_rk3(dn,eps);
   }

   Wact[in]=plaq_action_slices(data3.Wsl[in]);
   Yact[in]=ym_action_slices(data3.Ysl[in]);
   Qtop[in]=tcharge_slices(data3.Qsl[in]);
   if (bc_cstar()!=0)
   {
      for(x0=0;x0<N0;x0++)
      {
         data3.Wsl[in][x0]*=0.5;
         data3.Ysl[in][x0]*=0.5;
         data3.Qsl[in][x0]*=0.5;
      }
      Wact[in]*=0.5;
      Yact[in]*=0.5;
      Qtop[in]*=0.5;
   }

   cm3x3_assign(4*VOLUME,usv[0],udfld());
   set_flags(UPDATED_UD);
   release_wud();
}


void print_ms3dat(void)
{
   int in,dn,nn,din;
   double eps;

   init_ms3dat();

   dn=wfp.dn;
   nn=wfp.nn;
   eps=wfp.eps;

   din=nn/10;
   if (din<1)
      din=1;

   printf("SU(3) WF observables:\n\n");

   for (in=0;in<=nn;in+=din)
      printf("n = %3d, t = %.2e, Wact = %.6e, Yact = %.6e, Q = % .2e\n",
             in*dn,eps*(double)(in*dn),Wact[in],Yact[in],Qtop[in]);

   printf("\n");
}
