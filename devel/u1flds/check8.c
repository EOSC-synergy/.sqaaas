
/*******************************************************************************
*
* File check8.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Renormalization of the noncompact U(1) field.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1flds.h"
#include "global.h"

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
   error(rbuf==NULL,1,"main [check1.c]",
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
   int i,n,n_all;
   double d,dmax,dmax_all;
   double pi;
   double phi[2],phi_prime[2];
   double *ad;
   complex_dble *u1d,*u1dsv;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check8.log","w",stdout);

      printf("\n");
      printf("Renormalization of the noncompact U(1) field\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check8.c]",
                    "Syntax: check8 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check8.c]",
                    "Syntax: check8 [-bc <type>] [-cs <cstar>]");
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

   pi=4.0*atan(1.0);

   ad=adfld();

   ranlxd(ad,4*VOLUME);
   for (i=0;i<4*VOLUME;i++)
      ad[i]=(ad[i]-500.)*1000.;
   set_flags(UPDATED_AD);
   set_ad_bc();
   orbi_cpy_ad();
   
   
   n=0;
   for (i=0;i<4*VOLUME;i++)
      if((ad[i]<-pi)||(ad[i]>=pi)) n++;
   MPI_Reduce(&n,&n_all,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   if (my_rank==0)
   {
      printf("Number of links outside of [-pi,pi) before renormalization = %d\n",n_all);
   }
   error_root(n_all==0,1,"main [check8.c]","The gauge field should not be in the [0,2*pi) range here");
   
   
   u1d=u1dfld(LOC);
   u1dsv=malloc(4*VOLUME*sizeof(complex_dble));
   memcpy(u1dsv,u1d,4*VOLUME*sizeof(complex_dble));
   
   
   renormalize_ad();
   u1d=u1dfld(LOC);
   
   n=0;
   for (i=0;i<4*VOLUME;i++)
      if((ad[i]<-pi)||(ad[i]>=pi)) n++;
   MPI_Reduce(&n,&n_all,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   
   dmax=0.0;
   for (i=0;i<4*VOLUME;i++)
   {
      d=fabs(u1d[i].re-u1dsv[i].re);
      if (d>dmax) dmax=d;
      d=fabs(u1d[i].im-u1dsv[i].im);
      if (d>dmax) dmax=d;
   }
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Number of links outside of [-pi,pi) after renormalization = %d\n\n",n_all);
      printf("|u1d(before)-u1d(after)| = %.2e\n\n",dmax_all);
   }
   
   if(bc_cstar()!=0)
   {
      dmax_all=check_cstar_ad();
      if (my_rank==0)
         printf("C* boundary conditions after renormalization = %.2e\n\n",dmax_all);
   }

   print_flags();


   if (my_rank==0)
      fclose(flog);
   
   MPI_Finalize();
   exit(0);
}
