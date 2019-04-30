
/*******************************************************************************
*
* File ms5_utils.c
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
#include "u1flds.h"
#include "wflow.h"
#include "u1ftensor.h"
#include "linalg.h"
#include "main_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static struct
{
   int nt;
   double **Usl,**Msl,**Fsl[6];
} data5;

static double *Uact,*Mact,*Flux[6];

static wflow_parms_t wfp;


static void init_ms5dat(void)
{
   int nn,tmax;
   int in;
   double **pp,*p;
   static int init=0;
   
   if(init) return;

   error_root((gauge()&2)==0,1,"init_ms5dat [ms5_utils.c]",
              "U(1) gauge group is not active");

   wfp=wflow_parms();

   nn=wfp.nn;
   tmax=wfp.tmax;

   pp=amalloc(8*(nn+1)*sizeof(*pp),8);
   p=amalloc(8*(nn+1)*(tmax+1)*sizeof(*p),9);

   error_root((pp==NULL)||(p==NULL),1,"init_ms5dat [ms5_utils.c]",
              "Unable to allocate data arrays");

   data5.Usl=pp;
   data5.Msl=pp+nn+1;
   data5.Fsl[0]=pp+2*(nn+1);
   data5.Fsl[1]=pp+3*(nn+1);
   data5.Fsl[2]=pp+4*(nn+1);
   data5.Fsl[3]=pp+5*(nn+1);
   data5.Fsl[4]=pp+6*(nn+1);
   data5.Fsl[5]=pp+7*(nn+1);

   for (in=0;in<(8*(nn+1));in++)
   {
      *pp=p;
      pp+=1;
      p+=tmax;
   }

   Uact=p;
   p+=nn+1;
   Mact=p;
   p+=nn+1;
   Flux[0]=p;
   p+=nn+1;
   Flux[1]=p;
   p+=nn+1;
   Flux[2]=p;
   p+=nn+1;
   Flux[3]=p;
   p+=nn+1;
   Flux[4]=p;
   p+=nn+1;
   Flux[5]=p;
   
   init=1;
}


void write_ms5dat(FILE *fdat)
{
   int iw,nn,tmax;
   int in,k;

   init_ms5dat();

   iw=write_little_int(0,fdat,1,data5.nt);

   nn=wfp.nn;
   tmax=wfp.tmax;

   for (in=0;in<=nn;in++)
      iw+=write_little_dblearray(0,fdat,tmax,data5.Usl[in]);

   for (in=0;in<=nn;in++)
      iw+=write_little_dblearray(0,fdat,tmax,data5.Msl[in]);

   for (k=0;k<6;k++)
   {
      for (in=0;in<=nn;in++)
         iw+=write_little_dblearray(0,fdat,tmax,data5.Fsl[k][in]);
   }

   error_root(iw!=(1+8*(nn+1)*tmax),1,"write_ms5dat [ms5_utils.c]",
              "Incorrect write count");
}


int read_ms5dat(FILE *fdat,int *nt)
{
   int ir,nn,tmax;
   int in,k;

   init_ms5dat();

   ir=read_little_int(0,fdat,1,&(data5.nt));
   if (ir!=1) return 0;

   nn=wfp.nn;
   tmax=wfp.tmax;

   for (in=0;in<=nn;in++)
      ir+=read_little_dblearray(0,fdat,tmax,data5.Usl[in]);

   for (in=0;in<=nn;in++)
      ir+=read_little_dblearray(0,fdat,tmax,data5.Msl[in]);

   for (k=0;k<6;k++)
   {
      for (in=0;in<=nn;in++)
         ir+=read_little_dblearray(0,fdat,tmax,data5.Fsl[k][in]);
   }

   error_root(ir!=(1+8*(nn+1)*tmax),1,"read_ms5dat [ms5_utils.c]",
              "Read error or incomplete data record");

   (*nt)=data5.nt;

   return 1;
}


void write_ms5dat_head(FILE *fdat)
{
   int iw;

   init_ms5dat();

   iw=write_little_int(0,fdat,3,wfp.dn,wfp.nn,wfp.tmax);
   iw+=write_little_dble(0,fdat,1,wfp.eps);

   error_root(iw!=4,1,"write_ms5dat_head [ms5_utils.c]",
              "Incorrect write count");
}


void check_ms5dat_head(FILE *fdat)
{
   init_ms5dat();

   check_little_int("check_ms5dat_head [ms5_utils.c]",fdat,3,wfp.dn,wfp.nn,wfp.tmax);
   check_little_dble("check_ms5dat_head [ms5_utils.c]",fdat,1,wfp.eps);
}


void set_ms5dat(int nt)
{
   int flint,in,dn,nn,x0;
   double eps;
   double **asv;

   init_ms5dat();

   asv=reserve_wad(1);
   assign_dvec2dvec(4*VOLUME,adfld(),asv[0]);

   data5.nt=nt;
   dn=wfp.dn;
   nn=wfp.nn;
   eps=wfp.eps;
   flint=wfp.flint;

   for (in=0;in<nn;in++)
   {
      Uact[in]=u1_plaq_action_slices(data5.Usl[in]);
      Mact[in]=mxw_action_slices(data5.Msl[in]);
      Flux[0][in]=u1fluxes_slices(0,data5.Fsl[0][in]);
      Flux[1][in]=u1fluxes_slices(1,data5.Fsl[1][in]);
      Flux[2][in]=u1fluxes_slices(2,data5.Fsl[2][in]);
      Flux[3][in]=u1fluxes_slices(3,data5.Fsl[3][in]);
      Flux[4][in]=u1fluxes_slices(4,data5.Fsl[4][in]);
      Flux[5][in]=u1fluxes_slices(5,data5.Fsl[5][in]);
      if (bc_cstar()!=0)
      {
         for(x0=0;x0<N0;x0++)
         {
            data5.Usl[in][x0]*=0.5;
            data5.Msl[in][x0]*=0.5;
         }
         Uact[in]*=0.5;
         Mact[in]*=0.5;
      }

      if (flint==0)
         fwd_u1_euler(dn,eps);
      else if (flint==1)
         fwd_u1_rk2(dn,eps);
      else
         fwd_u1_rk3(dn,eps);
   }

   Uact[in]=u1_plaq_action_slices(data5.Usl[in]);
   Mact[in]=mxw_action_slices(data5.Msl[in]);
   Flux[0][in]=u1fluxes_slices(0,data5.Fsl[0][in]);
   Flux[1][in]=u1fluxes_slices(1,data5.Fsl[1][in]);
   Flux[2][in]=u1fluxes_slices(2,data5.Fsl[2][in]);
   Flux[3][in]=u1fluxes_slices(3,data5.Fsl[3][in]);
   Flux[4][in]=u1fluxes_slices(4,data5.Fsl[4][in]);
   Flux[5][in]=u1fluxes_slices(5,data5.Fsl[5][in]);
   if (bc_cstar()!=0)
   {
      for(x0=0;x0<N0;x0++)
      {
         data5.Usl[in][x0]*=0.5;
         data5.Msl[in][x0]*=0.5;
      }
      Uact[in]*=0.5;
      Mact[in]*=0.5;
   }

   assign_dvec2dvec(4*VOLUME,asv[0],adfld());
   set_flags(UPDATED_AD);
   release_wad();
}


void print_ms5dat(void)
{
   int in,dn,nn,din;
   double eps;

   init_ms5dat();

   dn=wfp.dn;
   nn=wfp.nn;
   eps=wfp.eps;

   din=nn/10;
   if (din<1)
      din=1;

   printf("U(1) WF observables:\n\n");

   for (in=0;in<=nn;in+=din)
      printf("n = %3d, t = %.2e, Uact = %.6e, Mact = %.6e, F(01,02,03,23,31,12) = ( %+.2e, %+.2e, %+.2e, %+.2e, %+.2e, %+.2e )\n",
             in*dn,eps*(double)(in*dn),Uact[in],Mact[in],Flux[0][in],Flux[1][in],Flux[2][in],Flux[3][in],Flux[4][in],Flux[5][in]);

   printf("\n");
}
