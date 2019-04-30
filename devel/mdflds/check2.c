
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2017 Agostino Patella
* 
* Based on openQCD-1.6/devel/uflds/check1.c
* Copyright (C) 2009, 2011, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check boundary conditions of random U(1) momentum.
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
#include "u1flds.h"
#include "mdflds.h"
#include "global.h"
#include "mdflds_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static double check_cstar_ad(void)
{
   int size,tag;
   double d,dmax,dmax_all;
   double *rbuf,*adb;
   double *cp1,*cp2;
   MPI_Status stat;
   
   if(bc_cstar()==0) return 0.0;

   if (query_flags(ADBUF_UP2DATE)!=1)
      copy_bnd_ad();

   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
      size+=3;
   
   rbuf=malloc(size*sizeof(double));
   error(rbuf==NULL,1,"main [check2.c]",
         "Unable to allocate rbuf");
   
   adb=adfld();
   
   tag=mpi_tag();
   MPI_Sendrecv(adb,size,MPI_DOUBLE,get_mirror_rank(),tag,
                rbuf,size,MPI_DOUBLE,get_mirror_rank(),tag,
                MPI_COMM_WORLD,&stat);

   dmax=0.0;
   cp1=adb;
   cp2=rbuf;
   for (;cp1<adb+size;cp1++)
   {
      d=fabs((*cp1)+(*cp2));
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
   double phi[2],phi_prime[2];
   double *ad,*adb,*adm;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);

      printf("\n");
      printf("Initialization of the U(1) fields\n");
      printf("---------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>] [-cs <type>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>] [-cs <type>]");
   }
   
   set_flds_parms(3,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,0,cs,phi,phi_prime);
   print_bc_parms();

   start_ranlux(0,123456);
   geometry();



   adb=adfld();
   random_u1mom();
   rot_ad(1.0);

   adm=adb+4*VOLUME;
   dmax1=0.0;
   dmax2=0.0;

   for (ad=adb;ad<adm;ad++)
   {
      iu=(ad-adb);
      ix=iu/8+VOLUME/2;
      ifc=iu%8;
      x0=global_time(ix);

      if ((bc==0)&&(((x0==0)&&(ifc==1))||((x0==(N0-1))&&(ifc==0))))
      {
         d2=fabs(*ad);
         if (d2>dmax2)
            dmax2=d2;
      }
      else if ((bc!=1)||(x0>0)||(ifc<2))
      {
         d1=fabs(*ad);
         if (d1>dmax1)
            dmax1=d1;
      }
      else
      {
         d2=fabs(*ad);
         if (d2>dmax2)
            dmax2=d2;
      }
   }

   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
   {
      ad=adb+4*VOLUME+7*(BNDRY/4);

      for (k=1;k<4;k++)
      {
         d2=fabs(*ad);
         ad+=1;

         if (d2>dmax2)
            dmax2=d2;
      }
   }

   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("enerate random U(1) momenta and shift zero U(1) gauge field\n");
      printf("|ad| = %.2e (this should not be small)\n",dmax1_all);
      if (bc!=3)
         printf("|ad-bval| = %.2e\n",dmax2_all);
   }

   if(bc_cstar()!=0)
   {
      dmax1_all=check_cstar_ad();
      if (my_rank==0)
         printf("C* boundary conditions ad = %.2e\n",dmax1_all);
   }
   
   if((bc_cstar()==0)&&(bc_type()==3))
   {
      if (my_rank==0)
         printf("Nothing to be checked\n");
   }

   if (my_rank==0)
   {
      printf("\n");
   }

   if (my_rank==0)
      fclose(flog);
   
   MPI_Finalize();
   exit(0);
}
