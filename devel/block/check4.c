
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of assign_swd2swbgr() and assign_swd2swdblk().
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "sw_term.h"
#include "block.h"
#include "global.h"
#include "gflds_utils.h"


static int cmp_sw(pauli *r,pauli *s)
{
   int i;

   for (i=0;i<36;i++)
   {
      if ((r[0].u[i]!=s[0].u[i])||(r[1].u[i]!=s[1].u[i]))
         return 1;
   }

   return 0;
}


static int cmp_swd(pauli_dble *r,pauli_dble *s)
{
   int i;

   for (i=0;i<36;i++)
   {
      if ((r[0].u[i]!=s[0].u[i])||(r[1].u[i]!=s[1].u[i]))
         return 1;
   }

   return 0;
}


static int check_sw(block_t *b)
{
   int x0,x1,x2,x3,x[4];
   int y0,y1,y2,y3,ix,iy;
   pauli *sw;

   sw=swfld();

   for (x0=0;x0<(*b).bs[0];x0++)
   {
      for (x1=0;x1<(*b).bs[1];x1++)
      {
         for (x2=0;x2<(*b).bs[2];x2++)
         {
            for (x3=0;x3<(*b).bs[3];x3++)
            {
               x[0]=x0;
               x[1]=x1;
               x[2]=x2;
               x[3]=x3;

               y0=(*b).bo[0]+x0;
               y1=(*b).bo[1]+x1;
               y2=(*b).bo[2]+x2;
               y3=(*b).bo[3]+x3;

               ix=ipt_blk(b,x);
               iy=ipt[y3+L3*y2+L2*L3*y1+L1*L2*L3*y0];

               if (cmp_sw((*b).sw+2*ix,sw+2*iy))
                     return 1;
            }
         }
      }
   }

   return 0;
}


static int check_swbgr(blk_grid_t grid)
{
   int nb,isw;
   block_t *b,*bm;

   b=blk_list(grid,&nb,&isw);
   bm=b+nb;

   for (;b<bm;b++)
   {
      if (check_sw(b))
         return 1;
   }

   return 0;
}


static int check_swd(block_t *b)
{
   int x0,x1,x2,x3,x[4];
   int y0,y1,y2,y3,ix,iy;
   pauli_dble *swd;

   swd=swdfld();

   for (x0=0;x0<(*b).bs[0];x0++)
   {
      for (x1=0;x1<(*b).bs[1];x1++)
      {
         for (x2=0;x2<(*b).bs[2];x2++)
         {
            for (x3=0;x3<(*b).bs[3];x3++)
            {
               x[0]=x0;
               x[1]=x1;
               x[2]=x2;
               x[3]=x3;

               y0=(*b).bo[0]+x0;
               y1=(*b).bo[1]+x1;
               y2=(*b).bo[2]+x2;
               y3=(*b).bo[3]+x3;

               ix=ipt_blk(b,x);
               iy=ipt[y3+L3*y2+L2*L3*y1+L1*L2*L3*y0];

               if (cmp_swd((*b).swd+2*ix,swd+2*iy))
                  return 1;
            }
         }
      }
   }

   return 0;
}


int main(int argc,char *argv[])
{
   int my_rank,bc,n,nb,isw,cf,qhat;
   int iset,ifail,bs[4];
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   ptset_t set;
   block_t *b;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of assign_swd2swbgr() and assign_swd2swdblk()\n");
      printf("---------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check4.c]",
                  "Syntax: check4 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;
   }

   MPI_Bcast(&cf,1,MPI_INT,0,MPI_COMM_WORLD);
   set_flds_parms(cf,0);
   print_flds_parms();

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,0,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,1234);
   geometry();

   set_sap_parms(bs,0,1,1);
   set_dfl_parms(bs,2);
   alloc_bgr(SAP_BLOCKS);
   alloc_bgr(DFL_BLOCKS);

   qhat=0;
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   if ((gauge()&2)!=0)
   {
      qhat=-3;
      u1csw=0.8;
   }
   if ((gauge()&1)!=0)
   {
      su3csw=0.95;
   }
   if (bc_type()!=3) cF[0]=cF[1]=1.0;
   if (bc_cstar()==0)
   {
      theta[0]=0.3;
      theta[1]=1.2;
      theta[2]=0.5;
   }
   set_dirac_parms9(qhat,-0.0123,su3csw,u1csw,cF[0],cF[1],
                    theta[0],theta[1],theta[2]);
   print_dirac_parms();

   random_gflds();
   ifail=0;

   for (iset=0;iset<(int)(PT_SETS);iset++)
   {
      if (iset==0)
         set=ALL_PTS;
      else if (iset==1)
         set=EVEN_PTS;
      else if (iset==2)
         set=ODD_PTS;
      else
         set=NO_PTS;

      sw_term(NO_PTS);
      ifail+=assign_swd2swbgr(SAP_BLOCKS,set);
      ifail+=sw_term(set);
      assign_swd2sw();
      error(check_swbgr(SAP_BLOCKS)!=0,1,"main [check4.c]",
            "assign_swd2swbgr() is incorrect");

      b=blk_list(DFL_BLOCKS,&nb,&isw);

      for (n=0;n<nb;n++)
      {
         sw_term(NO_PTS);
         ifail+=assign_swd2swdblk(DFL_BLOCKS,n,set);
         sw_term(set);
         error(check_swd(b+n)!=0,1,
               "main [check4.c]","assign_swd2swdblk() is incorrect");
      }
   }

   error(ifail!=0,1,"main [check4.c]",
         "Some of the inversions were not safe");

   if (my_rank==0)
   {
      printf("No errors detected\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
