
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2017 Nazario Tantalo, Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Checks of the program u1mom_Delta_no0 
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
#include "linalg.h"
#include "flags.h"
#include "dft.h"
#include "u1flds.h"
#include "mdflds.h"
#include "global.h"



int main(int argc,char *argv[])
{
   int bc,cs;
   size_t size;
   int ix,my_rank;
   double phi[2],phi_prime[2],u1phi,u1phi_prime;
   double *msv;
   double d,dev;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);

      printf("\n");
      printf("Check of the program u1mom_Delta_no0\n");
      printf("------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
      
      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check5.c]",
                    "Syntax: check5 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check5.c]",
                    "Syntax: check5 [-bc <type>] [-cs <cstar>]");
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

   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
      size+=3;
   msv=amalloc(size*sizeof(double),ALIGN);
   error(msv==NULL,1,"main [check5.c]",
         "Unable to allocate momentum fields");



   random_u1mom();
   u1mom_Delta_no0(0,mdflds()->u1mom,mdflds()->u1mom);
   assign_dvec2dvec(size,mdflds()->u1mom,msv);
   bnd_u1mom2zero();
   orbi_cpy_u1mom();
   muladd_assign_dvec(size,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
   {
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+0]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+1]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+2]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("\nboundary conditions after D(0)            dev= %e\n",dev);


   random_u1mom();
   u1mom_Delta_no0(1,mdflds()->u1mom,mdflds()->u1mom);
   assign_dvec2dvec(size,mdflds()->u1mom,msv);
   bnd_u1mom2zero();
   orbi_cpy_u1mom();
   muladd_assign_dvec(size,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
   {
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+0]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+1]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+2]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("boundary conditions after D(1)            dev= %e\n",dev);



   random_u1mom();
   u1mom_Delta_no0(2,mdflds()->u1mom,mdflds()->u1mom);
   assign_dvec2dvec(size,mdflds()->u1mom,msv);
   bnd_u1mom2zero();
   orbi_cpy_u1mom();
   muladd_assign_dvec(size,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
   {
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+0]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+1]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+2]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("boundary conditions after D(2)            dev= %e\n",dev);



   random_u1mom();
   u1mom_Delta_no0(3,mdflds()->u1mom,mdflds()->u1mom);
   assign_dvec2dvec(size,mdflds()->u1mom,msv);
   bnd_u1mom2zero();
   orbi_cpy_u1mom();
   muladd_assign_dvec(size,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
   {
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+0]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+1]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+2]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("boundary conditions after D(3)            dev= %e\n",dev);



   random_u1mom();
   u1mom_Delta_no0(4,mdflds()->u1mom,mdflds()->u1mom);
   assign_dvec2dvec(size,mdflds()->u1mom,msv);
   bnd_u1mom2zero();
   orbi_cpy_u1mom();
   muladd_assign_dvec(size,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
   {
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+0]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+1]);
      if(d>dev) dev=d;
      d=fabs(msv[4*VOLUME+7*(BNDRY/4)+2]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("boundary conditions after D(4)            dev= %e\n",dev);





   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,msv);
   u1mom_Delta_no0(0,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(0,mdflds()->u1frc,mdflds()->u1frc);
   u1mom_Delta_no0(0,msv,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,mdflds()->u1mom);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("\nlinearity of D(0)                         dev= %e\n",dev);



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,msv);
   u1mom_Delta_no0(1,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(1,mdflds()->u1frc,mdflds()->u1frc);
   u1mom_Delta_no0(1,msv,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,mdflds()->u1mom);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("linearity of D(1)                         dev= %e\n",dev);



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,msv);
   u1mom_Delta_no0(2,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(2,mdflds()->u1frc,mdflds()->u1frc);
   u1mom_Delta_no0(2,msv,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,mdflds()->u1mom);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("linearity of D(2)                         dev= %e\n",dev);



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,msv);
   u1mom_Delta_no0(3,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(3,mdflds()->u1frc,mdflds()->u1frc);
   u1mom_Delta_no0(3,msv,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,mdflds()->u1mom);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("linearity of D(3)                         dev= %e\n",dev);



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,msv);
   u1mom_Delta_no0(4,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(4,mdflds()->u1frc,mdflds()->u1frc);
   u1mom_Delta_no0(4,msv,msv);
   muladd_assign_dvec(4*VOLUME,1.0,mdflds()->u1frc,mdflds()->u1mom);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("linearity of D(4)                         dev= %e\n",dev);







   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   random_u1mom();
   u1mom_Delta_no0(0,mdflds()->u1mom,msv);
   dev=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1frc,msv);
   u1mom_Delta_no0(0,mdflds()->u1frc,msv);
   dev-=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1mom,msv);
   if (my_rank==0)
      printf("\nhermiticity of D(0)                       dev= %e\n",fabs(dev));



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   random_u1mom();
   u1mom_Delta_no0(1,mdflds()->u1mom,msv);
   dev=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1frc,msv);
   u1mom_Delta_no0(1,mdflds()->u1frc,msv);
   dev-=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1mom,msv);
   if (my_rank==0)
      printf("hermiticity of D(1)                       dev= %e\n",fabs(dev));



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   u1mom_Delta_no0(2,msv,msv);
   dev=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1frc,msv);
   assign_dvec2dvec(4*VOLUME,mdflds()->u1frc,msv);
   u1mom_Delta_no0(2,msv,msv);
   dev-=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1mom,msv);
   if (my_rank==0)
      printf("hermiticity of D(2)                       dev= %e\n",fabs(dev));



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   random_u1mom();
   u1mom_Delta_no0(3,mdflds()->u1mom,msv);
   dev=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1frc,msv);
   u1mom_Delta_no0(3,mdflds()->u1frc,msv);
   dev-=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1mom,msv);
   if (my_rank==0)
      printf("hermiticity of D(3)                       dev= %e\n",fabs(dev));



   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,mdflds()->u1frc);
   random_u1mom();
   u1mom_Delta_no0(4,mdflds()->u1mom,msv);
   dev=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1frc,msv);
   u1mom_Delta_no0(4,mdflds()->u1frc,msv);
   dev-=scalar_prod_dvec(4*VOLUME,1,mdflds()->u1mom,msv);
   if (my_rank==0)
      printf("hermiticity of D(4)                       dev= %e\n",fabs(dev));







   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(0,msv,msv);
   u1mom_Delta_no0(1,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("\nD(1) . D(0) == id                         dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(1,msv,msv);
   u1mom_Delta_no0(0,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(0) . D(1) == id                         dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(0,msv,msv);
   u1mom_Delta_no0(3,msv,msv);
   u1mom_Delta_no0(3,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(3)^2 . D(0) == id                       dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(3,msv,msv);
   u1mom_Delta_no0(3,msv,msv);
   u1mom_Delta_no0(0,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(0) . D(3)^2 == id                       dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(2,msv,msv);
   u1mom_Delta_no0(3,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(3) . D(2) == id                         dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(3,msv,msv);
   u1mom_Delta_no0(2,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(2) . D(3) == id                         dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(1,msv,msv);
   u1mom_Delta_no0(2,msv,msv);
   u1mom_Delta_no0(2,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(2)^2 . D(1) == id                       dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(2,msv,msv);
   u1mom_Delta_no0(2,msv,msv);
   u1mom_Delta_no0(1,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(1) . D(2)^2 == id                       dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(4,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(4) == id                                dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(0,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(2,msv,msv);
   u1mom_Delta_no0(2,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(0) == D(2)^2                            dev= %e\n",dev);
   
   
   
   random_u1mom();
   assign_dvec2dvec(4*VOLUME,mdflds()->u1mom,msv);
   
   u1mom_Delta_no0(1,mdflds()->u1mom,mdflds()->u1mom);
   u1mom_Delta_no0(3,msv,msv);
   u1mom_Delta_no0(3,msv,msv);
   muladd_assign_dvec(4*VOLUME,-1.0,mdflds()->u1mom,msv);
   
   dev=0.0;
   for(ix=0;ix<4*VOLUME;ix++)
   {
      d=fabs(msv[ix]);
      if(d>dev) dev=d;
   }
   
   d=dev;
   MPI_Allreduce(&d,&dev,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("D(1) == D(3)^2                            dev= %e\n",dev);



   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
