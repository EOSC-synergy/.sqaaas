
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2015 Martin Luescher
*               2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program dft4d() (random momenta).
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "utils.h"
#include "random.h"
#include "lattice.h"
#include "dft.h"
#include "global.h"

#define NTEST 8
#define NMOM 24

extern int dpgen(void);
extern void print_dft4d_parms(int id);

static int rcpr[4],nproc[4];
static int my_rank,csize,ism,ks[4];
static complex_dble *fs=NULL,*fts,*fs0,*fts0;


static void set_fs(int id)
{
   int m,nx[4];
   dft4d_parms_t *dp4d;

   dp4d=dft4d_parms(id);
   csize=(*dp4d).csize;
   nx[0]=(*dp4d).nx[0][rcpr[0]];
   nx[1]=(*dp4d).nx[1][rcpr[1]];
   nx[2]=(*dp4d).nx[2][rcpr[2]];
   nx[3]=(*dp4d).nx[3][rcpr[3]];

   if (fs!=NULL)
      afree(fs);

   m=csize*nx[0]*nx[1]*nx[2]*nx[3];
   fs=amalloc(2*(m+csize)*sizeof(*fs),4);
   error(fs==NULL,1,"set_f [check6.c]","Unable to allocate field array");

   fts=fs+m;
   fs0=fts+m;
   fts0=fs0+csize;
   gauss_dble((double*)(fs),4*(m+csize));
   ism=init_hsum(2*csize);
}


static void random_mom(int id)
{
   int mu,n[4];
   double r[4];
   dft4d_parms_t *dp4d;
   dft_parms_t *dp;

   dp4d=dft4d_parms(id);

   for (mu=0;mu<4;mu++)
   {
      dp=(*dp4d).dp[mu];
      n[mu]=(*dp).n+((*dp).type!=EXP);
   }

   ranlxd(r,4);

   for (mu=0;mu<4;mu++)
      ks[mu]=(int)((double)(n[mu])*r[mu]);

   MPI_Bcast(ks,4,MPI_INT,0,MPI_COMM_WORLD);
}


static void set_fs0(int id)
{
   int n[4],b[4],c[4],d[4];
   int i,mu,bo[4],nx[4],**nxa;
   int ix,iy,vol,x[4];
   double pi,r,r0[4];
   complex_dble pt,z,w;
   dft4d_parms_t *dp4d;
   dft_parms_t **dp;

   dp4d=dft4d_parms(id);
   dp=(*dp4d).dp;
   nxa=(*dp4d).nx;
   pi=4.0*atan(1.0);

   for (mu=0;mu<4;mu++)
   {
      n[mu]=dp[mu][0].n;
      b[mu]=dp[mu][0].b;
      c[mu]=dp[mu][0].c;
      d[mu]=(dp[mu][0].type==SIN);
      bo[mu]=0;

      for (i=0;i<rcpr[mu];i++)
         bo[mu]+=nxa[mu][i];

      nx[mu]=nxa[mu][rcpr[mu]];

      if (dp[mu][0].type==EXP)
         r0[mu]=0.5*pi/(double)(n[mu]);
      else
         r0[mu]=0.25*pi/(double)(n[mu]);
   }

   vol=nx[0]*nx[1]*nx[2]*nx[3];
   reset_hsum(ism);

   for (ix=0;ix<vol;ix++)
   {
      iy=ix;
      x[3]=bo[3]+iy%nx[3];
      iy/=nx[3];
      x[2]=bo[2]+iy%nx[2];
      iy/=nx[2];
      x[1]=bo[1]+iy%nx[1];
      x[0]=bo[0]+iy/nx[1];

      pt.re=1.0;
      pt.im=0.0;

      for (mu=0;mu<4;mu++)
      {
         z.re=0.0;
         z.im=0.0;
         r=r0[mu]*(double)((2*ks[mu]+b[mu])*(2*x[mu]+c[mu]));

         if (dp[mu][0].type==EXP)
         {
            z.re=cos(r);
            z.im=sin(r);
         }
         else if ((c[mu]==1)&&(x[mu]<n[mu]))
         {
            if (d[mu]==0)
               z.re=2.0*cos(r);
            else
               z.im=2.0*sin(r);
         }
         else if (c[mu]==0)
         {
            if (d[mu]==0)
            {
               if ((x[mu]>0)&&(x[mu]<n[mu]))
                  z.re=2.0*cos(r);
               else
                  z.re=cos(r);
            }
            else
            {
               if ((x[mu]>0)&&(x[mu]<n[mu]))
                  z.im=2.0*sin(r);
               else
                  z.im=sin(r);
            }
         }

         w.re=pt.re*z.re-pt.im*z.im;
         w.im=pt.re*z.im+pt.im*z.re;
         pt.re=w.re;
         pt.im=w.im;
      }

      for (i=0;i<csize;i++)
      {
         fs0[i].re=pt.re*fs[i+csize*ix].re-pt.im*fs[i+csize*ix].im;
         fs0[i].im=pt.re*fs[i+csize*ix].im+pt.im*fs[i+csize*ix].re;
      }

      add_to_hsum(ism,(double*)(fs0));
   }

   global_hsum(ism,(double*)(fs0));
}


static void set_fts0(int id)
{
   int **nxa;
   int i,mu,b,nx[4],n[4],k[4];
   int ip,mip,ik;
   dft4d_parms_t *dp4d;

   dp4d=dft4d_parms(id);
   nxa=(*dp4d).nx;

   for (mu=0;mu<4;mu++)
   {
      b=0;

      for (i=0;i<nproc[mu];i++)
      {
         k[mu]=ks[mu]-b;
         b+=nxa[mu][i];

         if (ks[mu]<b)
         {
            n[mu]=i;
            break;
         }
      }

      nx[mu]=nxa[mu][rcpr[mu]];
   }

   ip=ipr_global(n);

   if (my_rank==ip)
   {
      ik=k[3]+nx[3]*k[2]+nx[3]*nx[2]*k[1]+nx[3]*nx[2]*nx[1]*k[0];

      for (i=0;i<csize;i++)
      {
         fts0[i].re=fts[i+csize*ik].re;
         fts0[i].im=fts[i+csize*ik].im;
      }
   }
   MPI_Bcast(fts0,2*csize,MPI_DOUBLE,ip,MPI_COMM_WORLD);


   if(bc_cstar()!=0)
   {
      n[1]=(n[1]+NPROC1/2)%NPROC1;
      mip=ipr_global(n);
      
      if (my_rank==mip)
      {
         ik=k[3]+nx[3]*k[2]+nx[3]*nx[2]*k[1]+nx[3]*nx[2]*nx[1]*k[0];

         for (i=0;i<csize;i++)
         {
            fts0[i].re+=fts[i+csize*ik].re;
            fts0[i].im+=fts[i+csize*ik].im;
         }
      }
      MPI_Bcast(fts0,2*csize,MPI_DOUBLE,mip,MPI_COMM_WORLD);
   }
}


static double cmp1(int n,complex_dble *f0,complex_dble *f1)
{
   int i;
   double nrm,nmx,dev,dmx;

   nmx=0.0;
   dmx=0.0;

   for (i=0;i<n;i++)
   {
      nrm=fabs(f0[i].re)+fabs(f0[i].im);
      if (nrm>nmx)
         nmx=nrm;

      dev=fabs(f0[i].re-f1[i].re)+fabs(f0[i].im-f1[i].im);
      if (dev>dmx)
         dmx=dev;
   }

   if (nmx>0.0)
      return dmx/nmx;
   else
      return dmx;
}


int main(int argc,char *argv[])
{

   int mu,it,im,id,cs,bc;
   double phi[2],phi_prime[2],u1phi,u1phi_prime;
   double d,dmax,dall;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);

      printf("\n");
      printf("Check of the program dft4d() (random momenta)\n");
      printf("---------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
      
      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check6.c]",
                    "Syntax: check6 [-bc <type>] [-cs <cstar>]");

      fflush(flog);
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   u1phi=0.0;
   u1phi_prime=0.0;
   set_bc_parms(bc,cs,phi,phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();
   dall=0.0;

   nproc[0]=NPROC0;
   nproc[1]=NPROC1;
   nproc[2]=NPROC2;
   nproc[3]=NPROC3;
   if (bc_cstar()!=0)
      nproc[1]/=2;
   for (mu=0;mu<4;++mu)
      rcpr[mu]=cpr[mu]%nproc[mu];


   for (it=0;it<NTEST;it++)
   {
      id=dpgen();
      print_dft4d_parms(id);
      set_fs(id);
      inv_dft4d(id,fs,fs);
      dft4d(id,fs,fts);
      dmax=0.0;

      for (im=0;im<NMOM;im++)
      {
         random_mom(id);
         set_fs0(id);
         set_fts0(id);
         d=cmp1(csize,fts0,fs0);

         if (d>dmax)
            dmax=d;
      }

      if (my_rank==0)
         printf("Deviation = %.1e\n\n",dmax);

      if (dmax>dall)
         dall=dmax;
   }

   if (my_rank==0)
   {
      printf("Maximal deviation = %.1e\n\n",dall);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
