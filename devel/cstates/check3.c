
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check nabla_k A_k=0 for the Coulomb U(1) gauge field A^C_mu.
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
   double phi[2],phi_prime[2];
   double d,norm[3];
   double *camu,*div;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      printf("\n");
      printf("Check nabla_k A_k=0 for the Coulomb U(1) gauge field A^C_mu\n");
      printf("-----------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                     "Syntax: check3 -cs <cstar> [-bc <type>]");

      cs=find_opt(argc,argv,"-cs");
      error_root(cs==0 || sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check3.c]",
                  "Syntax: check3 -cs <cstar> [-bc <type>]");
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

   div=amalloc((VOLUME+BNDRY)*sizeof(*div),4);
   error(div==NULL,1,"main [check3.c]",
         "Unable to allocate service array");

   random_gflds();
   orbi_cpy_ad();

   camu=get_coulomb_amu();

   div_sym_dvec(camu+1*(VOLUME+BNDRY),camu+2*(VOLUME+BNDRY),camu+3*(VOLUME+BNDRY),div);
   d=norm_square_dvec(VOLUME,1,div);

   norm[0]=norm_square_dvec(VOLUME,1,camu+1*(VOLUME+BNDRY));
   norm[1]=norm_square_dvec(VOLUME,1,camu+2*(VOLUME+BNDRY));
   norm[2]=norm_square_dvec(VOLUME,1,camu+3*(VOLUME+BNDRY));

   if (my_rank==0)
   {   
      printf("|D_k A^{Coulomb}_k| / sqrt( sum_k |A^{Coulomb}_k|^2 ) = %.2e\n",sqrt(d/(norm[0]+norm[1]+norm[2])));
      printf("(should be around 1*10^(-15) or so)\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
