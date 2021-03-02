
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2016, 2017 Agostino Patella
* 
* Based on openQCD-1.6/devel/uflds/check4.c
* Copyright (C) 2005, 2007-2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs for the plaquette sums of the double-precision
* U(1) gauge field.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1flds.h"
#include "global.h"
#include "gflds_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static void random_vec(int *svec)
{
   int mu,bs[4];
   double r[4];

   bs[0]=NPROC0*L0;
   bs[1]=NPROC1*L1;
   bs[2]=NPROC2*L2;
   bs[3]=NPROC3*L3;

   ranlxd(r,4);

   for (mu=0;mu<4;mu++)
   {
      svec[mu]=(int)((double)(bs[mu])*r[mu]);
      if (svec[mu]>(bs[mu]/2))
         svec[mu]-=bs[mu];
   }

   MPI_Bcast(svec,4,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int bc,cs,my_rank,n,t,s[4];
   static double asl1[N0],asl2[N0];
   double su3phi[2],su3phi_prime[2];
   double u1phi,u1phi_prime;
   double qel;
   double act1,nplaq1,nplaq2,p1,p2;
   double d1,d2,d3;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);

      printf("\n");
      printf("Plaquette sums of the double-precision U(1) gauge field\n");
      printf("-------------------------------------------------------\n\n");

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
   }

   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   su3phi[0]=0.0;
   su3phi[1]=0.0;
   su3phi_prime[0]=0.0;
   su3phi_prime[1]=0.0;
   u1phi=0.0;
   u1phi_prime=0.0;
   if (cs==0)
   {
      u1phi=0.573;
      u1phi_prime=-1.827;
   }
   set_bc_parms(bc,cs,su3phi,su3phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   qel=1.111;
   set_u1lat_parms(0,4.348,1.0/qel,0.0,7.0,0.0,0.0,0);
   print_lat_parms();

   start_ranlux(0,12345);
   geometry();

   p1=u1_plaq_sum_dble(1);
   p2=u1_plaq_wsum_dble(1);

   if (bc==0)
   {
      nplaq1=(double)((6*N0-3)*N1)*(double)(N2*N3);
      nplaq2=(double)((6*N0-6)*N1)*(double)(N2*N3);
   }
   else if (bc==3)
   {
      nplaq1=(double)(6*N0*N1)*(double)(N2*N3);
      nplaq2=nplaq1;
   }
   else
   {
      nplaq1=(double)((6*N0+3)*N1)*(double)(N2*N3);
      nplaq2=(double)(6*N0*N1)*(double)(N2*N3);
   }
   d1=0.0;
   d2=0.0;

   if (bc==1)
   {
      d1=cos(u1phi/(double)(N1))+
         cos(u1phi/(double)(N2))+
         cos(u1phi/(double)(N3));

      d1=(d1-3.0)*(double)(N1*N2*N3);
   }

   if ((bc==1)||(bc==2))
   {
      d2=cos(u1phi_prime/(double)(N1))+
         cos(u1phi_prime/(double)(N2))+
         cos(u1phi_prime/(double)(N3));

      d2=(d2-3.0)*(double)(N1*N2*N3);
   }

   if (my_rank==0)
   {
      printf("After field initialization:\n");
      printf("Deviation from expected value (plaq_sum)  = %.1e\n",
             fabs(1.0-p1/(nplaq1+d1+d2)));
      printf("Deviation from expected value (plaq_wsum) = %.1e\n\n",
             fabs(1.0-p2/(nplaq2+d1+d2)));
   }

   print_flags();
   random_gflds();

   p1=u1_plaq_sum_dble(1);
   p2=u1_plaq_wsum_dble(1);
   act1=u1_plaq_action_slices(asl1);
   d1=act1;

   if ((bc==0)||(bc==3))
   {
      for (t=0;t<N0;t++)
         d1-=asl1[t];
   }

   if (my_rank==0)
   {
      printf("Comparison of plaq_wsum_dble() with plaq_action_slices():\n");
      printf("Absolute difference of total action = %.1e\n",
             fabs(nplaq2-qel*qel*act1-p2));
      if ((bc==0)||(bc==3))
         printf("Deviation from sum of action slices = %.1e\n\n",
                fabs(d1));
      else
         printf("\n");
   }

   random_g();
   transform_gflds();
   d1=fabs(p1-u1_plaq_sum_dble(1));
   d2=fabs(p2-u1_plaq_wsum_dble(1));
   u1_plaq_action_slices(asl2);
   d3=0.0;

   for (t=0;t<N0;t++)
      d3+=fabs(asl1[t]-asl2[t]);

   if (my_rank==0)
   {
      printf("Gauge invariance:\n");
      printf("Relative difference (plaq_sum_dble)  = %.1e\n",d1/fabs(p1));
      printf("Relative difference (plaq_wsum_dble) = %.1e\n",d2/fabs(p2));
      printf("Relative difference (action slices)  = %.1e\n\n",
             d3/((double)(N0)*asl2[1]));
   }

   if (my_rank==0)
      printf("Translation invariance:\n");

   random_gflds();
   p1=u1_plaq_sum_dble(1);
   p2=u1_plaq_wsum_dble(1);
   u1_plaq_action_slices(asl1);

   for (n=0;n<8;n++)
   {
      random_vec(s);
      if (bc!=3)
         s[0]=0;
      shift_gflds(s);
      d1=fabs(p1-u1_plaq_sum_dble(1));
      d2=fabs(p2-u1_plaq_wsum_dble(1));
      u1_plaq_action_slices(asl2);
      d3=0.0;

      for (t=0;t<N0;t++)
         d3+=fabs(asl1[safe_mod(t-s[0],N0)]-asl2[t]);

      for (t=0;t<N0;t++)
         asl1[t]=asl2[t];

      if (my_rank==0)
      {
         printf("s=(% 3d,% 3d,% 3d,% 3d):\n",s[0],s[1],s[2],s[3]);
         printf("Absolute deviation (plaq_sum_dble)  = %.1e\n",d1);
         printf("Absolute deviation (plaq_wsum_dble) = %.1e\n",d2);
         printf("Absolute difference (action slices) = %.1e\n\n",
                d2/(double)(N0));
      }
   }

   if (bc==1)
   {
      random_gflds();
      p1=u1_plaq_sum_dble(1);
      p2=u1_plaq_wsum_dble(1);

      if (my_rank==0)
      {
         printf("\n");
         printf("Comparison of plaq_sum_dble() and plaq_wsum_dble():\n");
         printf("Absolute deviation = %.1e\n\n",
                fabs(p1-p2-3.0*(double)(N1*N2*N3)));
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
