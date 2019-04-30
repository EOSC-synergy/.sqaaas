
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* Based on NSPT-1.4/devel/aflds/check2.c
* Copyright (C) 2015 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs gather_u1mom() and scatter_u1mom().
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
#include "flags.h"
#include "lattice.h"
#include "flags.h"
#include "u1flds.h"
#include "mdflds.h"
#include "global.h"

static int bc;
static complex_dble *mom[4];

static void random_flds(void)
{
   int i;

   for (i=0;i<4;i++)
   {
      mom[i]=amalloc(VOLUME*sizeof(*(mom[i])),ALIGN);
      error(mom[i]==NULL,1,"random_flds [check3.c]",
            "Unable to allocate memory space for service arrays");

      gauss_dble((double*)(mom[i]),2*VOLUME);
   }
}


static void set_lks(void)
{
   int x0,x1,x2,x3,ix,mu;
   int n[4],bo[4],nb1,x[4],r;
   double *X;
   mdflds_t *mdf;

   mdf=mdflds();

   n[0]=NPROC0*L0;
   n[1]=NPROC1*L1;
   n[2]=NPROC2*L2;
   n[3]=NPROC3*L3;

   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;

   nb1=(cpr[1]+NPROC1/2)%NPROC1;
   nb1*=L1;

   for (x0=0;x0<L0;x0++)
   for (x1=0;x1<L1;x1++)
   for (x2=0;x2<L2;x2++)
   for (x3=0;x3<L3;x3++)
   {
      ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

      if (ix>=(VOLUME/2))
      {
         x[0]=x0;
         x[1]=x1;
         x[2]=x2;
         x[3]=x3;

         X=mdf->u1mom+8*(ix-(VOLUME/2));

         for (mu=0;mu<4;mu++)
         {
            X[0]=(double)(mu+bo[0]+x[0]+bo[1]+x[1]+bo[2]+x[2]+bo[3]+x[3]);
            X+=1;

            r=x[mu];
            x[mu]-=1;
            if ((bo[mu]==0)&&(x[mu]<0))
            {
               x[mu]+=n[mu];
               if (mu<=bc_cstar() && mu>1)
               x[mu]+=nb1-bo[1];
            }

            X[0]=(double)(mu+bo[0]+x[0]+bo[1]+x[1]+bo[2]+x[2]+bo[3]+x[3]);
            X+=1;

            x[mu]=r;
         }
      }
   }
}


static void check_gather(complex_dble **ad)
{
   int mu,ix,iy,x0,x1,x2,x3;
   int bo[4],ie;
   complex_dble *X;

   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;
   ie=0;

   for (mu=0;mu<4;mu++)
   {
      X=ad[mu];

      for (ix=0;ix<VOLUME;ix++)
      {
         iy=ix;
         x3=iy%L3;
         iy/=L3;
         x2=iy%L2;
         iy/=L2;
         x1=iy%L1;
         iy/=L1;
         x0=iy;

         ie|=(X[0].re!=(double)(mu+bo[0]+x0+bo[1]+x1+bo[2]+x2+bo[3]+x3));
         ie|=(X[0].im!=0.0);

         X+=1;
      }

      error(ie!=0,1,"check_gather [check3.c]",
            "Field variables are not correctly copied");
   }
}


static void check_scatter(void)
{
   int x0,x1,x2,x3,ix,mu;
   int n[4],bo[4],nb1,x[4],r,ie;
   double *X;
   mdflds_t *mdf;

   mdf=mdflds();

   n[0]=NPROC0*L0;
   n[1]=NPROC1*L1;
   n[2]=NPROC2*L2;
   n[3]=NPROC3*L3;

   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;

   nb1=(cpr[1]+NPROC1/2)%NPROC1;
   nb1*=L1;

   ie=0;

   for (x0=0;x0<L0;x0++)
   for (x1=0;x1<L1;x1++)
   for (x2=0;x2<L2;x2++)
   for (x3=0;x3<L3;x3++)
   {
      ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

      if (ix>=(VOLUME/2))
      {
         x[0]=x0;
         x[1]=x1;
         x[2]=x2;
         x[3]=x3;

         X=mdf->u1mom+8*(ix-(VOLUME/2));

         for (mu=0;mu<4;mu++)
         {
            ie|=((*X)!=(double)(mu+bo[0]+x[0]+bo[1]+x[1]+bo[2]+x[2]+bo[3]+x[3]));
            X+=1;

            r=x[mu];
            x[mu]-=1;

            if ((bo[mu]==0)&&(x[mu]<0))
            {
               x[mu]+=n[mu];
               if (mu<=bc_cstar() && mu>1)
               x[mu]+=nb1-bo[1];
            }

            ie|=((*X)!=(double)(mu+bo[0]+x[0]+bo[1]+x[1]+bo[2]+x[2]+bo[3]+x[3]));
            X+=1;
            x[mu]=r;
         }
      }
   }

   error(ie!=0,1,"check_scatter [check3.c]",
         "Field variables are not correctly copied");
}



int main(int argc,char *argv[])
{
   int my_rank,cs;
   double phi[2],phi_prime[2],u1phi,u1phi_prime;
   FILE *flog=NULL;
   mdflds_t *mdf;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Check of the programs gather_u1mom() and scatter_u1mom()\n");
      printf("------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
      error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check3.c]",
                 "Syntax: check3 [-bc <type>] [-cs <cstar>]");
   }

   set_flds_parms(2,0);
   print_flds_parms();

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

   start_ranlux(0,1236);
   geometry();

   random_flds();

   set_lks();

   mdf=mdflds();

   gather_u1mom(mdf->u1mom,mom);
   check_gather(mom);

   scatter_u1mom(mom,mdf->u1mom);
   check_scatter();

   if (my_rank==0)
   {
      printf("No errors discovered\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
