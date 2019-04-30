
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2009, 2011, 2013, 2016 Martin Luescher
*               2017 Agostino Patella
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
#include "su3fcts.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)



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


static double dev_u3(su3_dble *u)
{
   su3_dble v;

   su3dagxsu3(u,u,&v);

   return dev_unity(&v);
}


static double dev_udet(int t,int ifc,su3_dble *u)
{
   double d,dmax;
   complex_dble det;

   det=udet(u);
   dmax=0.0;

   d=fabs(det.re-1.0);
   if (d>dmax)
      dmax=d;
   d=fabs(det.im);
   if (d>dmax)
      dmax=d;

   return dmax;
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


static double check_cstar(void)
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
   int iu,ix,ifc,x0,k,ie;
   double d1,d2,dmax1,dmax2;
   double dmax1_all,dmax2_all;
   double phi[2],phi_prime[2];
   su3 *u,*ub,*um;
   su3_dble vd,*ud,*udb,*udm;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Initialization of the link variables\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <cstar>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   if(cs==0)
   {
      phi[0]=0.123;
      phi[1]=-0.534;
      phi_prime[0]=0.912;
      phi_prime[1]=0.078;
   }
   set_bc_parms(bc,cs,phi,phi_prime,0.0,0.0);
   print_bc_parms();

   start_ranlux(0,123456);
   geometry();

   cm3x3_unity(1,&vd);
   ub=ufld();
   um=ub+4*VOLUME;
   dmax1=0.0;

   for (u=ub;u<um;u++)
   {
      d1=dev_udu(&vd,u);

      if (d1>dmax1)
         dmax1=d1;
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Allocate single-precision gauge field\n");
      printf("|u-1| = %.2e\n\n",dmax1_all);
   }

   print_flags();

   udb=udfld();
   ie=check_bc(0.0);
   error_root(ie==0,1,"main [check1.c]","Boundary conditions not properly set");

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
         d2=dev_bval(ifc/2,phi,ud);
         if (d2>dmax2)
            dmax2=d2;
      }
   }

   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
   {
      ud=udb+4*VOLUME+7*(BNDRY/4);

      for (k=1;k<4;k++)
      {
         d2=dev_bval(k,phi_prime,ud);
         ud+=1;

         if (d2>dmax2)
            dmax2=d2;
      }
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Allocate double-precision gauge field\n");
      printf("|ud-1| = %.2e\n",dmax1_all);
      if (bc!=3)
         printf("|ud-bval| = %.2e\n",dmax2_all);
      printf("\n");
   }

   print_flags();

   random_ud();
   assign_ud2u();
   ie=check_bc(0.0);
   error_root(ie==0,1,"main [check1.c]","Boundary conditions changed");

   ud=udb;
   udm=udb+4*VOLUME;
   u=ub;
   dmax1=0.0;

   for (ud=udb;ud<udm;ud++)
   {
      d1=dev_udu(ud,u);

      if (d1>dmax1)
         dmax1=d1;

      u+=1;
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Random fields\n");
      printf("Assign double-precision to single-precision field\n");
      printf("Maximal deviation = %.2e\n\n",dmax1_all);
   }

   print_flags();

   random_ud();
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
         d1=dev_zero(ud);
         d2=d1;
      }
      else
      {
         d1=dev_u3(ud);
         d2=dev_udet(x0,ifc,ud);
      }

      if (d1>dmax1)
         dmax1=d1;
      if (d2>dmax2)
         dmax2=d2;
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Call random_ud\n");
      printf("|u^dag*u-1| = %.2e\n",dmax1_all);
      printf("|det{u}-1| = %.2e\n",dmax2_all);
   }

   if(bc_cstar()!=0)
   {
      dmax1_all=check_cstar();

      if (my_rank==0)
         printf("C* boundary conditions = %.2e\n",dmax1_all);
   }
   
   if (my_rank==0)
      printf("\n");

   print_flags();

   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();
   exit(0);
}
