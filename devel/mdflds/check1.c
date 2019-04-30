
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/devel/uflds/check1.c
* Copyright (C) 2009, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check boundary conditions of random SU(3) momentum.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "su3fcts.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "mdflds.h"
#include "global.h"
#include "mdflds_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


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

   v=(*u);
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


static double check_cstar_ud(void)
{
   int size,tag;
   double d,dmax,dmax_all;
   su3_dble *rbuf,*udb;
   complex_dble *cp1,*cp2;
   MPI_Status stat;
   
   if(bc_cstar()==0) return 0.0;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
      size+=3;
   
   rbuf=malloc(size*sizeof(su3_dble));
   error(rbuf==NULL,1,"main [check1.c]",
         "Unable to allocate rbuf");
   
   udb=udfld();
   
   tag=mpi_tag();
   MPI_Sendrecv(udb,18*size,MPI_DOUBLE,get_mirror_rank(),tag,
                rbuf,18*size,MPI_DOUBLE,get_mirror_rank(),tag,
                MPI_COMM_WORLD,&stat);

   dmax=0.0;
   cp1=(complex_dble*)udb;
   cp2=(complex_dble*)rbuf;
   for (;cp1<(complex_dble*)udb+9*size;cp1++)
   {
      d=fabs((*cp1).re-(*cp2).re);
      if (d>dmax) dmax=d;
      d=fabs((*cp1).im+(*cp2).im);
      if (d>dmax) dmax=d;
      
      cp2++;
   }

   MPI_Allreduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   free(rbuf);
   
   return dmax_all;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cs;
   int iu,ix,ifc,x0,k;
   double d1,d2,dmax1,dmax2;
   double dmax1_all,dmax2_all;
   double su3phi[2],su3phi_prime[2];
   double u1phi,u1phi_prime;
   su3_dble *ud,*udb,*udm;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Boundary conditions of random SU(3) momentum\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <type>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <type>]");
   }
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   su3phi[0]=0.0;
   su3phi[1]=0.0;
   su3phi_prime[0]=0.0;
   su3phi_prime[1]=0.0;
   u1phi=0.0;
   u1phi_prime=0.0;
   if(cs==0)
   {
      su3phi[0]=0.123;
      su3phi[1]=-0.534;
      su3phi_prime[0]=0.912;
      su3phi_prime[1]=0.078;
      u1phi=0.573;
      u1phi_prime=-1.827;
   }
   set_bc_parms(bc,cs,su3phi,su3phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   start_ranlux(0,1236);
   geometry();
   
   
   udb=udfld();
   random_su3mom();
   rot_ud(1.0);
   
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
         d1=dev_unity(ud);
         if (d1>dmax1)
            dmax1=d1;
      }
      else
      {
         d2=dev_bval(ifc/2,su3phi,ud);
         if (d2>dmax2)
            dmax2=d2;
      }
   }

   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
   {
      ud=udb+4*VOLUME+7*(BNDRY/4);

      for (k=1;k<4;k++)
      {
         d2=dev_bval(k,su3phi_prime,ud);
         ud+=1;

         if (d2>dmax2)
            dmax2=d2;
      }
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Generate random SU(3) momenta and rotate unit SU(3) gauge field\n");
      printf("|ud-1| = %.2e (this should not be small)\n",dmax1_all);
      if (bc!=3)
         printf("|ud-bval| = %.2e\n",dmax2_all);
   }

   if(bc_cstar()!=0)
   {
      dmax1_all=check_cstar_ud();

      if (my_rank==0)
         printf("C* boundary conditions = %.2e\n",dmax1_all);
   }
   
   if((bc_cstar()==0)&&(bc_type()==3))
   {
      if (my_rank==0)
         printf("Nothing to be checked\n");
   }
   
   if (my_rank==0)
      printf("\n");

   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();
   exit(0);
}
