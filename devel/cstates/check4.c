
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Action of nabla_sq_dvec(), div_sym_dvec() and inv_nabla_sq()
* on plane waves.
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
#include "u1flds.h"
#include "u1ftensor.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "cstates.h"
#include "global.h"
#include "gflds_utils.h"
#include "sflds_utils.h"


int main(int argc,char *argv[])
{
   int my_rank,bc,cs;
   int n,i,ix,x0,x1,x2,x3;
   int np[4],bo[4];
   float ran[4];
   double pi,d[3],dmax[3],rs;
   double mp,pt,px,py,pz,p[4];
   double phi[2],phi_prime[2];
   double *phi1,*phi2,*phi3;
   double *phi4,*phi5,*phi6,*phi7;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      printf("\n");
      printf("Action of nabla_sq_dvec(), div_sym_dvec() and inv_nabla_sq() on plane waves\n");
      printf("--------------------------------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check5.c]",
                    "Syntax: check4 [-bc <type>] [-cs <cstar>]");
   }

   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,12345);
   geometry();

   phi1=amalloc(7*(VOLUME+BNDRY)*sizeof(*phi1),4);
   error(phi1==NULL,1,"main [check4.c]",
         "Unable to allocate service array");

   phi2=phi1+VOLUME+BNDRY;
   phi3=phi2+VOLUME+BNDRY;
   phi4=phi3+VOLUME+BNDRY;
   phi5=phi4+VOLUME+BNDRY;
   phi6=phi5+VOLUME+BNDRY;
   phi7=phi6+VOLUME+BNDRY;

   pi=4.0*atan(1.0);
   n=10;
   bo[0]=cpr[0]*L0;
   bo[1]=cpr[1]*L1;
   bo[2]=cpr[2]*L2;
   bo[3]=cpr[3]*L3;
   dmax[0]=0.0;
   dmax[1]=0.0;
   dmax[2]=0.0;

   for (i=0;i<n;i++)
   {
      ranlxs(ran,4);

      if (bc==0)
         np[0]=(int)(ran[0]*(float)(NPROC0*L0-1));
      else
         np[0]=(int)(ran[0]*(float)(NPROC0*L0));
      np[1]=(int)(ran[1]*(float)(NPROC1*L1));
      np[2]=(int)(ran[2]*(float)(NPROC2*L2));
      np[3]=(int)(ran[3]*(float)(NPROC3*L3));

      if (np[0]==0)
      np[0]=1;

      if (bc==0)
         p[0]=(double)(np[0])*pi/(double)(NPROC0*L0-1);
      else if (bc==3)
         p[0]=((double)(np[0])*2.0*pi+pi)/(double)(NPROC0*L0);
      else
         p[0]=(double)(np[0])*pi/(double)(NPROC0*L0);
      p[1]=(double)(np[1])*2.0*pi/(double)(NPROC1*L1);
      if((bc_cstar()<2)||(np[1]%2==0))
         p[2]=(double)(np[2])*2.0*pi/(double)(NPROC2*L2);
      else
         p[2]=(double)(2*np[2]+1)*pi/(double)(NPROC2*L2);
      if((bc_cstar()<3)||(np[1]%2==0))
         p[3]=(double)(np[3])*2.0*pi/(double)(NPROC3*L3);
      else
         p[3]=(double)(2*np[3]+1)*pi/(double)(NPROC3*L3);

      ranlxd(&rs,1);

      MPI_Bcast(p,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&rs,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      mp =(1.0-cos(p[1]));
      mp+=(1.0-cos(p[2]));
      mp+=(1.0-cos(p[3]));
      mp*=-2.0;

      for (x0=0;x0<L0;x0++)
      for (x1=0;x1<L1;x1++)
      for (x2=0;x2<L2;x2++)
      for (x3=0;x3<L3;x3++)
      {
         ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

         pt=p[0]*(double)(x0+bo[0]);
         px=p[1]*(double)(x1+bo[1]);
         py=p[2]*(double)(x2+bo[2]);
         pz=p[3]*(double)(x3+bo[3]);

         phi1[ix]=rs*cos(pt)*cos(px)*cos(py)*cos(pz);

         phi2[ix]=mp*phi1[ix];

         phi3[ix] =-sin(p[1])*rs*cos(pt)*sin(px)*cos(py)*cos(pz);
         phi3[ix]+=-sin(p[2])*rs*cos(pt)*cos(px)*sin(py)*cos(pz);
         phi3[ix]+=-sin(p[3])*rs*cos(pt)*cos(px)*cos(py)*sin(pz);

         phi4[ix]=phi1[ix]/mp;
      }

      nabla_sq_dvec(phi1,phi5);
      div_sym_dvec(phi1,phi1,phi1,phi6);

      inv_nabla_sq(phi7,phi1);

      muladd_assign_dvec(VOLUME,-1.0,phi5,phi2);
      d[0]=norm_square_dvec(VOLUME,1,phi2)/norm_square_dvec(VOLUME,1,phi1);
      d[0]=sqrt(d[0]);
      if (d[0]>dmax[0])
         dmax[0]=d[0];

      muladd_assign_dvec(VOLUME,-1.0,phi6,phi3);
      d[1]=norm_square_dvec(VOLUME,1,phi3)/norm_square_dvec(VOLUME,1,phi1);
      d[1]=sqrt(d[1]);
      if (d[1]>dmax[1])
         dmax[1]=d[1];

      muladd_assign_dvec(VOLUME,-1.0,phi7,phi4);
      d[2]=norm_square_dvec(VOLUME,1,phi4)/norm_square_dvec(VOLUME,1,phi1);
      d[2]=sqrt(d[2]);
      if (d[2]>dmax[2])
         dmax[2]=d[2];

      if (my_rank==0)
         printf("Normalized deviations at p=(%2d,%2d,%2d,%2d):\t d_lap= %.1e, d_div= %.1e, d_invlap= %.1e\n",
                  np[0],np[1],np[2],np[3],d[0],d[1],d[2]);
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Maximal normalized deviations: d_lap= %.1e, d_div= %.1e, d_invlap= %.1e\n\n",dmax[0],dmax[1],dmax[2]);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
