
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2008-2013, 2016 Martin Luescher, Filippo Palombi
*               2016 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the programs force0() and action0().
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
#include "mdflds.h"
#include "linalg.h"
#include "forces.h"
#include "global.h"
#include "gflds_utils.h"
#include "mdflds_utils.h"

#define N0 (NPROC0*L0)


static void unconstrained_random_su3mom(void)
{
   random_alg(4*VOLUME,mdflds()->su3mom);
   bnd_su3mom2zero();
}


static double dSdt(double c)
{
   int ie;
   mdflds_t *mdfs;

   mdfs=mdflds();
   ie=check_bnd_su3frc((*mdfs).su3mom);
   error(ie!=0,1,"dSdt [check3.c]",
         "Momentum field vanishes on an incorrect set of links");

   force0(c);
   ie=check_bnd_su3frc((*mdfs).su3frc);
   error(ie!=0,1,"dSdt [check3.c]",
         "Force field vanishes on an incorrect set of links");

   return scalar_prod_alg(4*VOLUME,0,(*mdfs).su3mom,(*mdfs).su3frc);
}


int main(int argc,char *argv[])
{
   int my_rank,k,ie,bc,sf,cs;
   double c,eps,act0,act1,dact,dsdt;
   double dev_frc,sig_loss,rdmy;
   double phi[2],phi_prime[2];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Check of the programs force0() and action0()\n");
      printf("--------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");

      sf=find_opt(argc,argv,"-sf");

      if (sf!=0)
         error_root(sscanf(argv[sf+1],"%d",&sf)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-sf <sf_type>] [-cs <cstar>]");
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
   set_bc_parms(bc,sf,cs,phi,phi_prime);
   print_bc_parms();

   set_su3lat_parms(3.50,0.95,0.82,1.32);

   start_ranlux(0,1234);
   geometry();
   alloc_wf3d(1);
   c=0.789;

   for (k=0;k<4;k++)
   {
      random_gflds();
      unconstrained_random_su3mom();
      dsdt=dSdt(c);

      eps=1.0e-4;
      rot_ud(eps);
      act0=2.0*action0(0)/3.0;
      rot_ud(-eps);

      rot_ud(-eps);
      act1=2.0*action0(0)/3.0;
      rot_ud(eps);

      rot_ud(2.0*eps);
      act0-=action0(0)/12.0;
      rot_ud(-2.0*eps);

      rot_ud(-2.0*eps);
      act1-=action0(0)/12.0;
      rot_ud(2.0*eps);

      act0*=c;
      act1*=c;

      dact=(act0-act1)/eps;
      dev_frc=dsdt-dact;
      sig_loss=-log10(fabs(1.0-act0/act1));

      rdmy=dsdt;
      MPI_Reduce(&rdmy,&dsdt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dsdt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=dev_frc;
      MPI_Reduce(&rdmy,&dev_frc,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&dev_frc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      rdmy=sig_loss;
      MPI_Reduce(&rdmy,&sig_loss,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&sig_loss,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      ie=check_bc(0.0);
      error_root(ie!=1,1,"main [check3.c]",
                 "Operations did not preserve boundary conditions");

      if (my_rank==0)
      {
         printf("Relative deviation of dS/dt = %.2e ",fabs(dev_frc/dsdt));
         printf("[significance loss = %d digits]\n",(int)(sig_loss));
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
