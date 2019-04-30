
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* Based on openQCD-1.4/develop/uflds/check1.c
* Copyright (C) 2009, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization of the link variables.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "su3fcts.h"
#include "uflds.h"
#include "u1flds.h"
#include "hflds.h"
#include "global.h"
#include "gflds_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

double phi[2],phi_prime[2];


static su3_dble pmat[4];
static complex_dble phase[4];
static complex_dble pdet[4];


static void set_phase(void)
{
   double r;
   bc_parms_t bcp;
   dirac_parms_t dp;

   pdet[0].re=1.0;
   pdet[0].im=0.0;

   bcp=bc_parms();
   dp=dirac_parms();

   if (bcp.type==3)
   {
      pdet[0].re=-1.0;
      pdet[0].im=0.0;
   }

   r=3.0*dp.theta[0]/(double)(N1);
   pdet[1].re=cos(r);
   pdet[1].im=sin(r);

   r=3.0*dp.theta[1]/(double)(N2);
   pdet[2].re=cos(r);
   pdet[2].im=sin(r);

   r=3.0*dp.theta[2]/(double)(N3);
   pdet[3].re=cos(r);
   pdet[3].im=sin(r);

   phase[0].re=1.0;
   phase[0].im=0.0;

   if (bcp.type==3)
   {
      phase[0].re=-1.0;
      phase[0].im=0.0;
   }

   r=dp.theta[0]/(double)(N1);
   phase[1].re=cos(r);
   phase[1].im=sin(r);

   r=dp.theta[1]/(double)(N2);
   phase[2].re=cos(r);
   phase[2].im=sin(r);

   r=dp.theta[2]/(double)(N3);
   phase[3].re=cos(r);
   phase[3].im=sin(r);
   
   cm3x3_zero(4,pmat);
   pmat[0].c11=pmat[0].c22=pmat[0].c33=phase[0];
   pmat[1].c11=pmat[1].c22=pmat[1].c33=phase[1];
   pmat[2].c11=pmat[2].c22=pmat[2].c33=phase[2];
   pmat[3].c11=pmat[3].c22=pmat[3].c33=phase[3];
}


static complex_dble udet(su3_dble *u)
{
   complex_dble det1,det2,det3,det;

   det1.re=
      ((*u).c22.re*(*u).c33.re-(*u).c22.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c32.re-(*u).c23.im*(*u).c32.im);
   det1.im=
      ((*u).c22.re*(*u).c33.im+(*u).c22.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c32.im+(*u).c23.im*(*u).c32.re);
   det2.re=
      ((*u).c21.re*(*u).c33.re-(*u).c21.im*(*u).c33.im)-
      ((*u).c23.re*(*u).c31.re-(*u).c23.im*(*u).c31.im);
   det2.im=
      ((*u).c21.re*(*u).c33.im+(*u).c21.im*(*u).c33.re)-
      ((*u).c23.re*(*u).c31.im+(*u).c23.im*(*u).c31.re);
   det3.re=
      ((*u).c21.re*(*u).c32.re-(*u).c21.im*(*u).c32.im)-
      ((*u).c22.re*(*u).c31.re-(*u).c22.im*(*u).c31.im);
   det3.im=
      ((*u).c21.re*(*u).c32.im+(*u).c21.im*(*u).c32.re)-
      ((*u).c22.re*(*u).c31.im+(*u).c22.im*(*u).c31.re);

   det.re=
      ((*u).c11.re*det1.re-(*u).c11.im*det1.im)-
      ((*u).c12.re*det2.re-(*u).c12.im*det2.im)+
      ((*u).c13.re*det3.re-(*u).c13.im*det3.im);
   det.im=
      ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
      ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
      ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   return det;
}


static double dev_zero(su3_dble *u)
{
   double d,dmax,*r,*rm;

   r=(double*)(u);
   rm=r+18;
   dmax=0.0;

   for (;r<rm;r++)
   {
      d=fabs(*r);
      if (d>dmax)
         dmax=d;
   }

   return dmax;
}


static double dev_unity(su3_dble *u)
{
   su3_dble v;

   v=(*u);
   v.c11.re-=1.0;
   v.c22.re-=1.0;
   v.c33.re-=1.0;

   return dev_zero(&v);
}


static double dev_bval(int k,double *phi,su3_dble *u)
{
   double s[3],phi3;

   su3_dble v;

   su3dagxsu3(pmat+k,u,&v);
   s[0]=(double)(N1);
   s[1]=(double)(N2);
   s[2]=(double)(N3);
   phi3=-phi[0]-phi[1];

   v.c11.re-=cos(phi[0]/s[k-1]);
   v.c11.im-=sin(phi[0]/s[k-1]);
   v.c22.re-=cos(phi[1]/s[k-1]);
   v.c22.im-=sin(phi[1]/s[k-1]);
   v.c33.re-=cos(phi3/s[k-1]);
   v.c33.im-=sin(phi3/s[k-1]);

   return dev_zero(&v);
}


static double dev_u3(su3_dble *u)
{
   su3_dble v;

   su3dagxsu3(u,u,&v);

   return dev_unity(&v);
}


static double dev_udu(su3_dble *ud,su3 *u)
{
   float *r,*rm;
   double d,dmax,*rd;

   rd=(double*)(ud);
   r=(float*)(u);
   rm=r+18;
   dmax=0.0;

   for (;r<rm;r++)
   {
      d=fabs((*rd)-(double)(*r));
      if (d>dmax)
         dmax=d;
      rd+=1;
   }

   return dmax;
}


static double dev_udud(su3_dble *u,su3_dble *v)
{
   double d,dmax,*s,*r,*rm;

   s=(double*)(u);
   r=(double*)(v);
   rm=r+18;
   dmax=0.0;

   for (;r<rm;r++)
   {
      d=fabs((*s)-(*r));
      if (d>dmax)
         dmax=d;
      s+=1;
   }

   return dmax;
}


static double dev_init(int t,int ifc,su3_dble *u)
{
   double dmax;
   su3_dble w;

   cm3x3_unity(1,&w);
   if ((ifc>1)||((ifc==0)&&(t==(N0-1)))||((ifc==1)&&(t==0)))
   {
      w=pmat[ifc/2];
   }

   dmax=dev_udud(u,&w);
   
   return dmax;
}


static double dev_hd(int t,int ifc,su3_dble *h,su3_dble *u,double *a)
{
   int q;
   double dmax;
   su3_dble w;
   complex_dble u1q;

   q=dirac_parms().qhat;
   
   u1q.re=cos(q*(*a));
   u1q.im=sin(q*(*a));
   if ((ifc>1)||((ifc==0)&&(t==(N0-1)))||((ifc==1)&&(t==0)))
      u1xu1(phase+ifc/2,&u1q,&u1q);
   
   u1xu1(&u1q,&((*u).c11),&(w.c11));
   u1xu1(&u1q,&((*u).c12),&(w.c12));
   u1xu1(&u1q,&((*u).c13),&(w.c13));
   u1xu1(&u1q,&((*u).c21),&(w.c21));
   u1xu1(&u1q,&((*u).c22),&(w.c22));
   u1xu1(&u1q,&((*u).c23),&(w.c23));
   u1xu1(&u1q,&((*u).c31),&(w.c31));
   u1xu1(&u1q,&((*u).c32),&(w.c32));
   u1xu1(&u1q,&((*u).c33),&(w.c33));
   
   dmax=dev_udud(h,&w);
   
   return dmax;
}


static double dev_hdet(int t,int ifc,su3_dble *h,double *a)
{
   int q;
   double d,dmax;
   complex_dble w,det;

   q=dirac_parms().qhat;

   w.re=cos(3*q*(*a));
   w.im=sin(3*q*(*a));
   if ((ifc>1)||((ifc==0)&&(t==(N0-1)))||((ifc==1)&&(t==0)))
      u1xu1(pdet+ifc/2,&w,&w);

   det=udet(h);
   dmax=0.0;

   d=fabs(w.re-det.re);
   if (d>dmax)
      dmax=d;
   d=fabs(w.im-det.im);
   if (d>dmax)
      dmax=d;

   return dmax;
}


static void check_hdinit(double *dev1,double *dev2)
{
   int iu,ix,ifc,x0,bc;
   double d1,d2,dmax1,dmax2;
   su3_dble *ud,*udb,*udm;
   double phi0[2]={0.0,0.0};

   bc=bc_type();
   
   udb=hdfld();

   udm=udb+4*VOLUME;
   dmax1=0.0;
   dmax2=0.0;

   for (ud=udb;ud<udm;ud++)
   {
      iu=(ud-udb);
      ix=iu/8+VOLUME/2;
      ifc=iu%8;
      x0=global_time(ix);

      if ((bc==0)&&(((x0==0)&&(ifc==1))||((x0==(N0-1))&&(ifc==0))))
      {
         d2=dev_zero(ud);
         if (d2>dmax2)
            dmax2=d2;
      }
      else if ((bc!=1)||(x0>0)||(ifc<2))
      {
         d1=dev_init(x0,ifc,ud);
         if (d1>dmax1)
            dmax1=d1;
      }
      else
      {
         if((gauge()&1)!=0)
            d2=dev_bval(ifc/2,phi,ud);
         else
            d2=dev_bval(ifc/2,phi0,ud);
         if (d2>dmax2)
            dmax2=d2;
      }
   }

   MPI_Reduce(&dmax1,dev1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,dev2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
}


static void check_hdh(double *dev1)
{
   double d1,dmax1;
   su3_dble *ud,*udb,*udm;
   su3 *u;

   udb=hdfld();
   udm=udb+4*VOLUME;
   u=hfld();
   dmax1=0.0;

   for (ud=udb;ud<udm;ud++)
   {
      d1=dev_udu(ud,u);

      if (d1>dmax1)
         dmax1=d1;

      u+=1;
   }

   MPI_Reduce(&dmax1,dev1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
}


static void check_hhdag(double *dev1,double *dev2)
{
   int iu,ix,ifc,x0,bc;
   double d1,d2,dmax1,dmax2;
   su3_dble *ud,*udb,*udm;

   bc=bc_type();

   udb=hdfld();
   udm=udb+4*VOLUME;
   dmax1=0.0;
   dmax2=0.0;

   for (ud=udb;ud<udm;ud++)
   {
      iu=(ud-udb);
      ix=iu/8+VOLUME/2;
      ifc=iu%8;
      x0=global_time(ix);

      if ((bc==0)&&(((x0==0)&&(ifc==1))||((x0==(N0-1))&&(ifc==0))))
      {
         d2=dev_zero(ud);
         if (d2>dmax2)
            dmax2=d2;
      }
      else
      {
         d1=dev_u3(ud);
         if (d1>dmax1)
            dmax1=d1;
      }

   }

   MPI_Reduce(&dmax1,dev1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,dev2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
}


static void check_deth(double *dev1,double *dev2)
{
   int iu,ix,ifc,x0,bc;
   double d1,d2,dmax1,dmax2;
   su3_dble *ud,*udb,*udm;
   double *ad;

   bc=bc_type();
   
   udb=hdfld();
   ad=adfld();

   udm=udb+4*VOLUME;
   dmax1=0.0;
   dmax2=0.0;

   for (ud=udb;ud<udm;ud++)
   {
      iu=(ud-udb);
      ix=iu/8+VOLUME/2;
      ifc=iu%8;
      x0=global_time(ix);

      if ((bc==0)&&(((x0==0)&&(ifc==1))||((x0==(N0-1))&&(ifc==0))))
      {
         d2=dev_zero(ud);
         if (d2>dmax2)
            dmax2=d2;
      }
      else
      {
         d1=dev_hdet(x0,ifc,ud,ad);
         if (d1>dmax1)
            dmax1=d1;
      }
      
      ad++;
   }

   MPI_Reduce(&dmax1,dev1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,dev2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
}


static double cmp_hd(su3_dble *usv)
{
   double *a,*b,*am;
   double d1,dmax1;
   
   a=(double*)hdfld();
   b=(double*)usv;
   am=a+18*4*VOLUME;
   dmax1=0.0;
   for(;a<am;a++)
   {
      d1=fabs((*a)-(*b));
      if (d1>dmax1)
         dmax1=d1;
      b++;
   }
   
   d1=dmax1;
   MPI_Reduce(&d1,&dmax1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(&dmax1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return dmax1;
}


static double cmp_hdud(void)
{
   int iu,ix,x0,ifc;
   double d1,dmax1,a;
   su3_dble *hdb,*hd,*hdm,*ud;
   
   hdb=hdfld();

   ud=udfld();
   hdm=hdb+4*VOLUME;
   dmax1=0.0;
   for(hd=hdb;hd<hdm;hd++)
   {
      iu=(hd-hdb);
      ix=iu/8+VOLUME/2;
      ifc=iu%8;
      x0=global_time(ix);
      
      a=0.0;
      d1=dev_hd(x0,ifc,hd,ud,&a);
      if (d1>dmax1)
         dmax1=d1;
      
      ud++;
   }
   
   d1=dmax1;
   MPI_Reduce(&d1,&dmax1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Bcast(&dmax1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return dmax1;
}


int main(int argc,char *argv[])
{
   int my_rank,cf,q,bc;
   int iu,ix,ifc,x0;
   double su3csw,u1csw,cF[2],theta[3];
   double dev1,dev2;
   FILE *flog=NULL;
   su3_dble *usv,*hsv,*u,*h,m,w;
   double *asv,*a;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Group properties of the U(3) gauge field\n");
      printf("----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check1.c]",
                  "Syntax: check1 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check1.c]",
                  "Syntax: check1 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      }
      else
         q=-3;
   }

   MPI_Bcast(&cf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&q,1,MPI_INT,0,MPI_COMM_WORLD);
   set_flds_parms(cf,0);
   print_flds_parms();
   if(gauge()==1) q=0;

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   if (gauge()==2)
   {
      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
   }   
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();
   
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   if ((gauge()&1)!=0) su3csw=0.95;
   if ((gauge()&2)!=0) u1csw=0.8;
   if (bc_type()!=3)
   {
      cF[0]=1.301;
      cF[1]=0.789;
   }
   if (bc_cstar()==0)
   {
      theta[0]=0.35;
      theta[1]=-1.25;
      theta[2]=0.78;
   }
   set_dirac_parms9(q,-0.0123,su3csw,u1csw,cF[0],cF[1],
                    theta[0],theta[1],theta[2]);
   print_dirac_parms();

   start_ranlux(0,123456);
   geometry();
   
   set_phase();

   print_flags();
   
   check_hdinit(&dev1,&dev2);
   if (my_rank==0)
   {
      printf("Initialization double-precision U(3) field\n");
      printf("|hd-e^(ith)| = %.2e\n",dev1);
      if (bc!=3)
         printf("|hd-bval| = %.2e\n",dev2);
      printf("\n");
   }

   print_flags();
   
   
   if((gauge()&1)!=0)
   {
      random_ud();

      dev1=cmp_hdud();

      if (my_rank==0)
      {
         printf("Random ud, and ad=0\n");
         printf("|hd-e^(ith)*ud| = %.2e\n\n",dev1);
      }

      print_flags();
   }


   random_ud();
   random_ad();

   check_hdh(&dev1);
   if (my_rank==0)
   {
      printf("Random fields\n");
      printf("|hd-h| = %.2e\n\n",dev1);
   }

   print_flags();


   random_ud();
   random_ad();

   check_hhdag(&dev1,&dev2);

   if (my_rank==0)
   {
      printf("Random fields\n");
      printf("|h^dag*h-1| = %.2e\n",dev1);
      if (bc==0)
         printf("|hd-bval| = %.2e\n",dev2);
      printf("\n");
   }
   
   print_flags();


   random_ud();
   random_ad();

   check_deth(&dev1,&dev2);

   if (my_rank==0)
   {
      printf("Random fields\n");
      printf("|det h-u1^(3q)*e^(3ith)| = %.2e\n",dev1);
      if (bc==0)
         printf("|det h-bval| = %.2e\n",dev2);
      printf("\n");
   }
   
   print_flags();


   random_ud();
   random_ad();
   u=udfld();
   a=adfld();
   h=hdfld();
   usv=malloc(sizeof(*usv)*4*VOLUME);
   asv=malloc(sizeof(*asv)*4*VOLUME);
   hsv=malloc(sizeof(*hsv)*4*VOLUME);
   for(iu=0;iu<4*VOLUME;iu++)
   {
      usv[iu]=u[iu];
      asv[iu]=a[iu];
      hsv[iu]=h[iu];
   }
   
   random_ud();
   random_ad();
   hdfld();
   for(iu=0;iu<4*VOLUME;iu++)
   {
      ix=iu/8+VOLUME/2;
      ifc=iu%8;
      x0=global_time(ix);
      
      cm3x3_unity(1,&w);
      if ((ifc>1)||((ifc==0)&&(x0==(N0-1)))||((ifc==1)&&(x0==0)))
      {
         w=pmat[ifc/2];
      }
      
      su3xsu3(usv+iu,u+iu,&m);
      u[iu]=m;
      
      a[iu]+=asv[iu];
      
      su3xsu3(hsv+iu,h+iu,&m);
      su3dagxsu3(&w,&m,hsv+iu);
   }
   set_flags(UPDATED_UD);
   set_flags(UPDATED_AD);
   
   dev1=cmp_hd(hsv);
   
   if (my_rank==0)
   {
      printf("Group property h(u,a)*h(u',a') = e^(ith) h(u*u',a+a')\n");
      printf("Maximum deviation = %.2e\n",dev1);
      printf("\n");
   }
   
   print_flags();


   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();
   exit(0);
}
