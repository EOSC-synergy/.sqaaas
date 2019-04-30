
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*               2016 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge and translation invariance of the SU(3) gauge action.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "forces.h"
#include "global.h"
#include "gflds_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int bc;


static double bnd_action(void)
{
   int i,j;
   double c0,c1,*cG,*phi;
   double s[3],d0[2],d1[2],act;
   su3lat_parms_t lat;
   bc_parms_t bcp;

   if ((bc==1)||(bc==2))
   {
      lat=su3lat_parms();
      bcp=bc_parms();

      s[0]=(double)(N1);
      s[1]=(double)(N2);
      s[2]=(double)(N3);

      for (i=0;i<2;i++)
      {
         d0[i]=0.0;
         d1[i]=0.0;
         phi=bcp.phi3[i];

         for (j=0;j<3;j++)
         {
            d0[i]-=(cos(phi[0]/s[j])+cos(phi[1]/s[j])+
                    cos(phi[2]/s[j])-3.0);
            d1[i]-=(cos(2.0*phi[0]/s[j])+cos(2.0*phi[1]/s[j])+
                    cos(2.0*phi[2]/s[j])-3.0);
         }
      }

      c0=lat.c0;
      c1=lat.c1;
      cG=lat.cG;

      act=c0*cG[1]*d0[1]+c1*d0[1]+c1*1.5*d1[1];

      if (bc==1)
         act+=(c0*cG[0]*d0[0]+c1*d0[0]+c1*1.5*d1[0]);

      return (lat.beta/3.0)*(double)(N1*N2*N3)*act;
   }
   else
      return 0.0;
}


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
   int my_rank,n,s[4],sf,cs;
   double phi[2],phi_prime[2],p1,p2;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Gauge and translation invariance of the gauge action\n");
      printf("----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");

      sf=find_opt(argc,argv,"-sf");

      if (sf!=0)
         error_root(sscanf(argv[sf+1],"%d",&sf)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");
   }
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&sf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   set_su3lat_parms(3.50,0.95,0.82,1.32,sf);
   print_lat_parms();

   start_ranlux(0,12345);
   geometry();

   p1=action0(1);
   p2=bnd_action();

   if (my_rank==0)
   {
      printf("Action after initialization = %.15e\n",p1);
      printf("Expected value              = %.15e\n\n",p2);
   }

   random_gflds();
   p1=action0(1);
   random_g();
   transform_gflds();
   p2=action0(1);

   if (my_rank==0)
   {
      printf("Random gauge field:\n");
      printf("Action = %.12e\n",p1);
      printf("Gauge invariance: relative difference = %.1e\n\n",
             fabs(1.0-p2/p1));
   }

   if (my_rank==0)
      printf("Translation invariance:\n");

   p1=action0(1);

   for (n=0;n<8;n++)
   {
      random_vec(s);
      if (bc!=3)
         s[0]=0;
      shift_gflds(s);
      p2=action0(1);

      if (my_rank==0)
      {
         printf("s=(% d, % d,% d,% d), ",s[0],s[1],s[2],s[3]);
         printf("relative deviation = %.1e\n",fabs(1.0-p2/p1));
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
