
/*******************************************************************************
*
* File check5.c
*
* Copyright (C) 2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency between U(1) and SU(3) boundary staple fields.
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
#include "u1flds.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)



int main(int argc,char *argv[])
{
   int my_rank,bc,cs,nfc[8],ifc;
   size_t size;
   double phi[2],phi_prime[2];
   complex_dble *u1d,*u1db,*u1dm,*u1bstb,*u1bstm,*u1bst;
   su3_dble *ud,*bst;
   double d,dmax,dmax_all;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check5.log","w",stdout);

      printf("\n");
      printf("Consistency between U(1) and SU(3) boundary staple fields\n");
      printf("---------------------------------------------------------\n\n");

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
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();

   if (NPROC==1)
   {
      if (my_rank==0)
      {
         printf("Nothing to be checked!\n\n");
         fclose(flog);
      }

      MPI_Finalize();
      exit(0);
   }

   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
      size+=3;

   random_ad();
   u1db=u1dfld(EXT);
   u1dm=u1db+size;
   ud=udfld();
   cm3x3_zero(4*VOLUME,ud);
   for (u1d=u1db;u1d<u1dm;u1d++)
   {
      (*ud).c11=(*u1d);
      ud++;
   }
   set_flags(UPDATED_UD);
      
   set_bstap();
   set_u1_bstap();
   
   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   size=0;
   for (ifc=0;ifc<8;ifc+=2)
   {
      if (size<nfc[ifc])
         size=nfc[ifc];
   }
   size=3*(BNDRY+2*size);
   
   
   u1bstb=u1_bstap();
   u1bstm=u1bstb+size;
   
   dmax=0.0;
   bst=bstap();
   for (u1bst=u1bstb;u1bst<u1bstm;u1bst++)
   {
      d=fabs((*bst).c11.re-(*u1bst).re);
      if (d>dmax)
         dmax=d;
      d=fabs((*bst).c11.im-(*u1bst).im);
      if (d>dmax)
         dmax=d;
      bst++;
   }
   
   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Consistency between set_bstap() and u1_bstap(), dev = %.2e\n",dmax_all);
      printf("\n");
   }
   

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
