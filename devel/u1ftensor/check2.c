
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2016 Nazario Tantalo
*               2017 Agostino Patella
*
* Based on openQCD-1.6/devel/tcharge/check5.c
* Copyright (C) 2010, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Maxwell action of background fields with constant field tensor.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "u1ftensor.h"
#include "global.h"
#include "gflds_utils.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int bc,np[4],bo[4];
static double mt[4][4],inp[4],twopi;


static double afld(int *x,int mu)
{
   int nu;
   double xt[4],phi;

   xt[0]=(double)(safe_mod(x[0],N0));
   xt[1]=(double)(safe_mod(x[1],N1));
   xt[2]=(double)(safe_mod(x[2],N2));
   xt[3]=(double)(safe_mod(x[3],N3));
   
   phi=0.0;

   for (nu=0;nu<mu;nu++)
      phi-=inp[nu]*mt[mu][nu]*xt[nu];

   phi*=inp[mu];

   if (safe_mod(x[mu],np[mu])==(np[mu]-1))
   {
      for (nu=(mu+1);nu<4;nu++)
         phi-=inp[nu]*mt[mu][nu]*xt[nu];      
   }

   return twopi*phi;
}


static double ftplaq(int *x,int mu,int nu)
{
   double sm,om,*phi;
   bc_parms_t bcp;

   bcp=bc_parms();
   
   if ((x[0]==0)&&(mu==0)&&(bc==1))
   {
      sm=afld(x,mu);
      x[mu]+=1;
      sm+=afld(x,nu);
      x[mu]-=1;
      x[nu]+=1;
      sm-=afld(x,mu);
      x[nu]-=1;

      phi=bcp.phi3[0];
      om=sm-phi[0]*inp[nu];
   }
   else if ((x[0]==(N0-1))&&(mu==0)&&((bc==1)||(bc==2)))
   {
      sm=afld(x,mu)-afld(x,nu);
      x[nu]+=1;
      sm-=afld(x,mu);
      x[nu]-=1;

      phi=bcp.phi3[1];
      om=sm+phi[0]*inp[nu];
   }
   else
   {
      sm=afld(x,mu)-afld(x,nu);
      x[mu]+=1;
      sm+=afld(x,nu);
      x[mu]-=1;
      x[nu]+=1;
      sm-=afld(x,mu);
      x[nu]-=1;

      om=sm;
   }

   return sin(om);
}


static double Abnd(void)
{
   int ib,x1,x2,x3,x[4];
   int mu,nu;
   double ft;
   double r0,r1,r2,r3;
   double aloc,aall;

   aloc=0.0;

   for (ib=0;ib<2;ib++)
   {
      if (ib==0)
         x[0]=1;
      else
         x[0]=N0-1;
      
      if (((ib==0)&&(bc==1)&&(cpr[0]==0))||
          ((ib==1)&&((bc==1)||(bc==2))&&(cpr[0]==(NPROC0-1))))
      {
         for (x1=0;x1<L1;x1++)
         {
            for (x2=0;x2<L2;x2++)
            {
               for (x3=0;x3<L3;x3++)
               {
                  x[1]=bo[1]+x1;
                  x[2]=bo[2]+x2;
                  x[3]=bo[3]+x3;

                  for (mu=0;mu<3;mu++)
                  {
                     for (nu=(mu+1);nu<4;nu++)
                     {
                        r0=ftplaq(x,mu,nu);
                        x[mu]-=1;
                        r1=ftplaq(x,mu,nu);
                        x[nu]-=1;
                        r2=ftplaq(x,mu,nu);
                        x[mu]+=1;
                        r3=ftplaq(x,mu,nu);
                        x[nu]+=1;

                        ft=0.25*(r0+r1+r2+r3);

                        aloc+=ft*ft;
                     }
                  }
               }
            }
         }
      }
   }
   
   if (NPROC>1)
   {
      MPI_Reduce(&aloc,&aall,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&aall,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      return aall;
   }
   else
      return aloc;
}


static double Amt(void)
{
   int mu,nu;
   double sm,pi,iqel2;
   double xl[4],phi,ft;
   u1lat_parms_t u1lp;

   u1lp=u1lat_parms();
   iqel2 =u1lp.invqel;
   iqel2*=u1lp.invqel;

   xl[0]=(double)(NPROC0*L0);
   xl[1]=(double)(NPROC1*L1);
   xl[2]=(double)(NPROC2*L2);
   xl[3]=(double)(NPROC3*L3);

   pi=4.0*atan(1.0);
   sm=0.0;

   for (mu=1;mu<4;mu++)
   {
      for (nu=0;nu<mu;nu++)
      {
         phi=2.0*pi*mt[mu][nu]/(xl[mu]*xl[nu]);
         
         ft=sin(phi);

         sm+=ft*ft;
      }
   }

   if (bc==0)
      sm*=(double)((N0-2)*N1)*(double)(N2*N3);
   else if (bc==1)
   {
      sm*=(double)((N0-3)*N1)*(double)(N2*N3);
      sm+=Abnd();
   }
   else if (bc==2)
   {
      sm*=(double)((N0-2)*N1)*(double)(N2*N3);
      sm+=Abnd();
   }
   else
      sm*=(double)(N0*N1)*(double)(N2*N3);

   return 0.5*iqel2*sm; 
}


static void choose_mt(void)
{
   int mu,nu;
   double r[6];

   ranlxd(r,6);
   MPI_Bcast(r,6,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   mt[0][1]=(double)((int)(3.0*r[0])-1);
   mt[0][2]=(double)((int)(3.0*r[1])-1);
   mt[0][3]=(double)((int)(3.0*r[2])-1);
   mt[1][2]=(double)((int)(3.0*r[3])-1);
   mt[1][3]=(double)((int)(3.0*r[4])-1);
   mt[2][3]=(double)((int)(3.0*r[5])-1);   

   for (mu=0;mu<4;mu++)
   {
      mt[mu][mu]=0.0;
            
      for (nu=0;nu<mu;nu++)
         mt[mu][nu]=-mt[nu][mu];
   }
}


static void set_ad(void)
{
   int x[4];
   int x0,x1,x2,x3;
   int ix,ifc;
   double phi;
   double *adb,*a;

   adb=adfld();
   
   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               if (ix>=(VOLUME/2))
               {
                  x[0]=bo[0]+x0;
                  x[1]=bo[1]+x1;
                  x[2]=bo[2]+x2;
                  x[3]=bo[3]+x3;                  
                  
                  a=adb+8*(ix-(VOLUME/2));

                  for (ifc=0;ifc<8;ifc++)
                  {
                     if (ifc&0x1)
                        x[ifc/2]-=1;

                     phi=afld(x,ifc/2);

                     if (ifc&0x1)
                        x[ifc/2]+=1;
                     
		     (*a)=phi;
                     a+=1;
                  }
               }
            }
         }
      }
   }
   set_ad_bc();

   set_flags(UPDATED_AD);   
}


int main(int argc,char *argv[])
{
   int my_rank,i;
   double phi[2],phi_prime[2];   
   double A1,A2,d,dmax;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
      printf("\n");
      printf("Maxwell action of background fields with constant field tensor\n");
      printf("--------------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }

   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,phi,phi_prime,0.0,0.0);
   print_bc_parms();

   set_u1lat_parms(0,1.0/137.0,5.0,0.0,7.0,0.0,0.0,0);
   print_lat_parms();

   start_ranlux(0,123);
   geometry();

   twopi=8.0*atan(1.0);
   
   np[0]=N0;
   np[1]=N1;
   np[2]=N2;
   np[3]=N3;

   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;
   
   inp[0]=1.0/(double)(np[0]);
   inp[1]=1.0/(double)(np[1]);
   inp[2]=1.0/(double)(np[2]);
   inp[3]=1.0/(double)(np[3]);    
   
   dmax=0.0;
   
   for (i=0;i<10;i++)
   {
      choose_mt();
      set_ad();


      A1=Amt();
      A2=mxw_action();

      if (my_rank==0)
         printf("Field no = %2d, A1 = %12.6e, A2 = %12.6e\n",i+1,A1,A2);

      d=fabs(A1-A2)/A1;
      if (d>dmax)
         dmax=d;
   }


   if (my_rank==0)
   {
      printf("\n");
      printf("Maximal relative deviation = %.1e\n\n",dmax);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}

#undef N0
#undef N1
#undef N2
#undef N3
