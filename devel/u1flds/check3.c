
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency between U(1) and SU(3) shifts.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "su3fcts.h"
#include "u1flds.h"
#include "global.h"

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
   int my_rank,bc,cs;
   int s[4],n;
   double d,dmax,dmax_all;
   double phi[2],phi_prime[2];
   double *ad,*adb,*adm;
   su3_dble *ud;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Consistency between U(1) and SU(3) shifts\n");
      printf("-----------------------------------------\n\n");

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
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,123456);
   geometry();
   
   
   dmax=0.0;
   for (n=0;n<8;n++)
   {
      random_ad();
      random_vec(s);
      if (bc!=3)
         s[0]=0;

      adb=adfld();
      adm=adb+4*VOLUME;
      ud=udfld();
      cm3x3_zero(4*VOLUME,ud);
      for (ad=adb;ad<adm;ad++)
      {
         (*ud).c11.re=(*ad);
         (*ud).c11.im=1.0;
         
         (*ud).c22.re=(*ad)/((*ad)*(*ad)+1.0);
         (*ud).c22.im=-1.0/((*ad)*(*ad)+1.0);
         
         (*ud).c33.re=1.0;

         ud++;
      }
      set_flags(UPDATED_UD);
      
      shift_ad(s);
      shift_ud(s);

      ud=udfld();
      for (ad=adb;ad<adm;ad++)
      {
         d=fabs((*ud).c11.re-(*ad));
         if (d>dmax)
            dmax=d;
         ud++;
      }
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between shift_ad and shift_ud, dev = %.2e\n",dmax_all);
      printf("\n");
   }

   print_flags();

   if (my_rank==0)
      fclose(flog);
   MPI_Finalize();
   exit(0);
}
