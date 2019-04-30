
/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/devel/forces/time1.c
* Copyright (C) 2005, 2008-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of elementary update steps().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "mdflds.h"
#include "linalg.h"
#include "global.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static void update_su3mom(void)
{
   int bc,cs,ix,t,ifc;
   int mirror, tag;
   su3_alg_dble *mom,*frc;
   mdflds_t *mdfs;
   MPI_Status stat;

   bc=bc_type();
   cs=bc_cstar();
   mdfs=mdflds();
   
   if (cs!=0)
   {
      mom=(*mdfs).su3mom;
      for (ix=0;ix<4*VOLUME;ix++)
      {
         _su3_alg_mul_assign((*mom),0.5);
         mom+=1;
      }
   }

   mom=(*mdfs).su3mom;
   frc=(*mdfs).su3frc;
   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         _su3_alg_sub_assign((*mom),(*frc));
         mom+=1;
         frc+=1;

         if (bc!=0)
         {
            _su3_alg_sub_assign((*mom),(*frc));
         }

         mom+=1;
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               _su3_alg_sub_assign((*mom),(*frc));
            }

            mom+=1;
            frc+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
         {
            _su3_alg_sub_assign((*mom),(*frc));
         }

         mom+=1;
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            _su3_alg_sub_assign((*mom),(*frc));
            mom+=1;
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            _su3_alg_sub_assign((*mom),(*frc));
            mom+=1;
            frc+=1;
         }
      }
   }

   if(cs>0) {
      mirror=get_mirror_rank();
      tag=mpi_tag();

      mom=mdflds()->su3mom;
      frc=mdflds()->su3frc;

      MPI_Sendrecv(mom,8*4*VOLUME,MPI_DOUBLE,mirror,tag,
                   frc,8*4*VOLUME,MPI_DOUBLE,mirror,tag,
                   MPI_COMM_WORLD,&stat);

      for (ix=0;ix<4*VOLUME;ix++)
      {
         (*mom).c1-=(*frc).c1;
         (*mom).c2-=(*frc).c2;
         (*mom).c3+=(*frc).c3;
         (*mom).c4-=(*frc).c4;
         (*mom).c5+=(*frc).c5;
         (*mom).c6-=(*frc).c6;
         (*mom).c7+=(*frc).c7;
         (*mom).c8-=(*frc).c8;
         mom+=1;
         frc+=1;
      }
   }
}


static void update_u1mom(void)
{
   int bc,cs,ix,t,ifc;
   int mirror, tag;
   double *mom,*frc;
   mdflds_t *mdfs;
   MPI_Status stat;

   bc=bc_type();
   cs=bc_cstar();
   mdfs=mdflds();
   
   if (cs!=0)
   {
      mom=(*mdfs).u1mom;
      for (ix=0;ix<4*VOLUME;ix++)
      {
         (*mom)*=0.5;
         mom+=1;
      }
   }

   mom=(*mdfs).u1mom;
   frc=(*mdfs).u1frc;
   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         (*mom)-=(*frc);
         mom+=1;
         frc+=1;

         if (bc!=0)
         {
            (*mom)-=(*frc);
         }

         mom+=1;
         frc+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               (*mom)-=(*frc);
            }

            mom+=1;
            frc+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
         {
            (*mom)-=(*frc);
         }

         mom+=1;
         frc+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            (*mom)-=(*frc);
            mom+=1;
            frc+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            (*mom)-=(*frc);
            mom+=1;
            frc+=1;
         }
      }
   }

   if(cs>0) {
      mirror=get_mirror_rank();
      tag=mpi_tag();

      mom=mdflds()->u1mom;
      frc=mdflds()->u1frc;

      MPI_Sendrecv(mom,4*VOLUME,MPI_DOUBLE,mirror,tag,
                   frc,4*VOLUME,MPI_DOUBLE,mirror,tag,
                   MPI_COMM_WORLD,&stat);

      for (ix=0;ix<4*VOLUME;ix++)
      {
         (*mom)-=(*frc);
         mom+=1;
         frc+=1;
      }
   }
}


static void update_ud(double eps)
{
   int bc,ix,t,ifc;
   su3_dble *u;
   su3_alg_dble *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).su3mom;
   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         if (bc!=0)
            expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            expXsu3(eps,mom,u);
         u+=1;
         mom+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            expXsu3(eps,mom,u);
            u+=1;
            mom+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


static void update_ad(double eps)
{
   int bc,ix,t,ifc;
   double *a;
   double *mom;
   mdflds_t *mdfs;

   bc=bc_type();
   mdfs=mdflds();
   mom=(*mdfs).u1mom;
   a=adfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         (*a)+=eps*(*mom);
         a+=1;
         mom+=1;

         if (bc!=0)
            (*a)+=eps*(*mom);
         a+=1;
         mom+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
               (*a)+=eps*(*mom);
            a+=1;
            mom+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc!=0)
            (*a)+=eps*(*mom);
         a+=1;
         mom+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            (*a)+=eps*(*mom);
            a+=1;
            mom+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            (*a)+=eps*(*mom);
            a+=1;
            mom+=1;
         }
      }
   }

   set_flags(UPDATED_AD);
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,n,count,k;
   double phi[2],phi_prime[2];
   double wt1,wt2,wdt,eps;
   mdflds_t *mdfs;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);

      printf("\n");
      printf("Timing of elementary update steps\n");
      printf("---------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   
   set_flds_parms(3,0);
   
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);

   print_flds_parms();
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();

   mdfs=mdflds();

   random_su3mom();
   assign_alg2alg(4*VOLUME,(*mdfs).su3mom,(*mdfs).su3frc);
   
   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
      {
         update_su3mom();
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
   {
      printf("Time per link\n");
      printf("update_su3mom():      %4.6f usec\n",wdt);
   }

   eps=1.0;
   for (k=0;k<4;k++)
   {
      random_su3mom();
      random_ud();
      
      n=(int)(3.0e6/(double)(4*VOLUME));
      if (n<2)
         n=2;
      wdt=0.0;

      while (wdt<5.0)
      {
         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();
         for (count=0;count<n;count++)
         {
            update_ud(eps);
         }
         MPI_Barrier(MPI_COMM_WORLD);
         wt2=MPI_Wtime();

         wdt=wt2-wt1;
         n*=2;
      }

      wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

      if (my_rank==0)
      {
         printf("update_ud(%.3f):     %4.6f usec\n",eps,wdt);
      }
      
      eps*=.1;
   }

   random_u1mom();
   assign_dvec2dvec(4*VOLUME,(*mdfs).u1mom,(*mdfs).u1frc);
   
   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
      {
         update_u1mom();
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
   {
      printf("update_u1mom():       %4.6f usec\n",wdt);
   }

   eps=1.0;
   random_u1mom();
   random_ad();
   
   n=(int)(3.0e6/(double)(4*VOLUME));
   if (n<2)
      n=2;
   wdt=0.0;

   while (wdt<5.0)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();
      for (count=0;count<n;count++)
      {
         update_ad(eps);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      wdt=wt2-wt1;
      n*=2;
   }

   wdt=2.0e6*wdt/((double)(n)*(double)(4*VOLUME));

   if (my_rank==0)
   {
      printf("update_ad(%.3f):     %4.6f usec\n",eps,wdt);
   }
   
   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
