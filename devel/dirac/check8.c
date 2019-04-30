
/*******************************************************************************
*
* File check8.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Comparison of Dw_blk_dble(),..,Dwhat_blk_dble with Dw_dble().
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
#include "u1flds.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "block.h"
#include "dirac.h"
#include "global.h"
#include "gflds_utils.h"


static void blk_sd2zero(int ic,spinor_dble *sd)
{
   int nb,isw;
   int nbh,n,nm,vol;
   block_t *b;

   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   if (ic^isw)
      n=nbh;
   else
      n=0;

   nm=n+nbh;

   for (;n<nm;n++)
   {
      set_sd2zero(vol,b[n].sd[0]);
      assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,0,sd);
   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cf,q;
   int nb,isw,nbh,ic,itm;
   int bs[4],n,nm,vol,volh,ie;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   double mu,d,dmax;
   spinor_dble **psd;
   block_t *b;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check8.log","w",stdout);
      fin=freopen("check7.in","r",stdin);

      printf("\n");
      printf("Comparison of Dw_blk_dble(),..,Dwhat_blk_dble() "
             "with Dw_dble()\n");
      printf("------------------------------------------------"
             "--------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check8.c]",
                    "Syntax: check8 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check8.c]",
                  "Syntax: check8 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check8.c]",
                  "Syntax: check8 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      }
      else
         q=-3;
   }

   MPI_Bcast(&cf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&q,1,MPI_INT,0,MPI_COMM_WORLD);
   set_flds_parms(cf,0);
   print_flds_parms();
   if(gauge()==1) q=0;

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   if (gauge()==2)
   {
      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
   }   
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();

   start_ranlux(0,1234);
   geometry();
   set_dfl_parms(bs,2);
   alloc_bgr(DFL_BLOCKS);
   alloc_wsd(4);

   mu=0.123;
   su3csw=u1csw=0.0;
   cF[0]=cF[1]=0.0;
   theta[0]=theta[1]=theta[2]=0.0;
   if ((gauge()&1)!=0) su3csw=0.95;
   if ((gauge()&2)!=0) u1csw=0.8;
   if (bc_type()!=3)
   {
      cF[0]=1.301;
      cF[1]=0.789;
   }
   if (bc_cstar()==0)
   {
      theta[0]=0.35;
      theta[1]=-1.25;
      theta[2]=0.78;
   }
   set_dirac_parms9(q,-0.0123,su3csw,u1csw,cF[0],cF[1],
                       theta[0],theta[1],theta[2]);
   print_dirac_parms();

   random_gflds();
   sw_term(NO_PTS);

   psd=reserve_wsd(4);
   b=blk_list(DFL_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;
   volh=vol/2;

   for (itm=0;itm<2;itm++)
   {
      ie=0;
      dmax=0.0;
      set_tm_parms(itm);

      if (my_rank==0)
         printf("Twisted-mass flag = %d\n",itm);

      for (ic=0;ic<2;ic++)
      {
         random_sd(VOLUME,psd[0],1.0);
         random_sd(VOLUME,psd[2],1.0);
         blk_sd2zero(ic^0x1,psd[0]);
         blk_sd2zero(ic^0x1,psd[2]);

         if (ic^isw)
            n=nbh;
         else
            n=0;

         nm=n+nbh;

         for (;n<nm;n++)
         {
            assign_ud2udblk(DFL_BLOCKS,n);
            assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);

            random_sd(vol,b[n].sd[1],1.0);
            assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
            Dw_blk_dble(DFL_BLOCKS,n,mu,0,1);
            assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,0,psd[2]);
            assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,1,psd[3]);
         }

         Dw_dble(mu,psd[0],psd[1]);
         blk_sd2zero(ic^0x1,psd[1]);
         blk_sd2zero(ic^0x1,psd[3]);

         mulr_spinor_add_dble(VOLUME,psd[0],psd[2],-1.0);
         mulr_spinor_add_dble(VOLUME,psd[1],psd[3],-1.0);

         if (norm_square_dble(VOLUME,0,psd[0])!=0.0)
            ie=1;

         d=norm_square_dble(VOLUME,1,psd[1])/
            norm_square_dble(VOLUME,1,psd[3]);

         if (d>dmax)
            dmax=d;
      }

      error(ie,1,"main [check8.c]",
            "Dw_blk_dble() changes the fields where it should not");

      dmax=sqrt(dmax);

      if (my_rank==0)
      {
         printf("The maximal relative deviations are:\n\n");
         printf("Dw_blk_dble():                    %.1e\n",dmax);
      }

      dmax=0.0;
      random_sd(VOLUME,psd[0],1.0);
      random_sd(VOLUME,psd[1],1.0);

      for (n=0;n<nb;n++)
      {
         assign_ud2udblk(DFL_BLOCKS,n);
         assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);

         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[1],1);
         Dwee_blk_dble(DFL_BLOCKS,n,mu,0,1);
         assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,0,psd[2]);
         assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,1,psd[3]);

         assign_sd2sdblk(DFL_BLOCKS,n,EVEN_PTS,psd[0],0);
         assign_sd2sdblk(DFL_BLOCKS,n,ODD_PTS,psd[1],0);
         Dwee_blk_dble(DFL_BLOCKS,n,mu,0,0);
         mulr_spinor_add_dble(vol,b[n].sd[0],b[n].sd[1],-1.0);
         if (norm_square_dble(vol,0,b[n].sd[0])!=0.0)
            ie=1;
      }

      Dwee_dble(mu,psd[0],psd[1]);
      mulr_spinor_add_dble(VOLUME,psd[0],psd[2],-1.0);
      mulr_spinor_add_dble(VOLUME,psd[1],psd[3],-1.0);

      if (norm_square_dble(VOLUME,0,psd[0])!=0.0)
         ie=1;

      d=norm_square_dble(VOLUME,1,psd[1])/
         norm_square_dble(VOLUME,1,psd[3]);

      if (d>dmax)
         dmax=d;

      random_sd(VOLUME,psd[0],1.0);
      random_sd(VOLUME,psd[1],1.0);

      for (n=0;n<nb;n++)
      {
         assign_ud2udblk(DFL_BLOCKS,n);
         assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);

         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[1],1);
         Dwoo_blk_dble(DFL_BLOCKS,n,mu,0,1);
         assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,0,psd[2]);
         assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,1,psd[3]);

         assign_sd2sdblk(DFL_BLOCKS,n,ODD_PTS,psd[0],0);
         assign_sd2sdblk(DFL_BLOCKS,n,EVEN_PTS,psd[1],0);
         Dwoo_blk_dble(DFL_BLOCKS,n,mu,0,0);
         mulr_spinor_add_dble(vol,b[n].sd[0],b[n].sd[1],-1.0);
         if (norm_square_dble(vol,0,b[n].sd[0])!=0.0)
            ie=1;
      }

      Dwoo_dble(mu,psd[0],psd[1]);
      mulr_spinor_add_dble(VOLUME,psd[0],psd[2],-1.0);
      mulr_spinor_add_dble(VOLUME,psd[1],psd[3],-1.0);

      if (norm_square_dble(VOLUME,0,psd[0])!=0.0)
         ie=1;

      d=norm_square_dble(VOLUME,1,psd[1])/
         norm_square_dble(VOLUME,1,psd[3]);

      if (d>dmax)
         dmax=d;

      error(ie,1,"main [check8.c]","Dwee_blk_dble() or Dwoo_blk_dble() "
            "changes the fields where it should not");

      dmax=sqrt(dmax);

      if (my_rank==0)
         printf("Dwee_blk_dble(), Dwoo_blk_dble(): %.1e\n",dmax);

      dmax=0.0;

      for (ic=0;ic<2;ic++)
      {
         random_sd(VOLUME,psd[0],1.0);
         random_sd(VOLUME,psd[1],1.0);
         random_sd(VOLUME,psd[2],1.0);
         blk_sd2zero(ic^0x1,psd[0]);
         blk_sd2zero(ic^0x1,psd[2]);

         if (ic^isw)
            n=nbh;
         else
            n=0;

         nm=n+nbh;

         for (;n<nm;n++)
         {
            assign_ud2udblk(DFL_BLOCKS,n);
            assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);

            assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
            assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[1],1);
            Dweo_blk_dble(DFL_BLOCKS,n,0,1);
            assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,0,psd[2]);
            assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,1,psd[3]);
         }

         Dweo_dble(psd[0],psd[1]);
         blk_sd2zero(ic^0x1,psd[1]);
         blk_sd2zero(ic^0x1,psd[3]);

         mulr_spinor_add_dble(VOLUME,psd[0],psd[2],-1.0);
         mulr_spinor_add_dble(VOLUME,psd[1],psd[3],-1.0);

         if (norm_square_dble(VOLUME,0,psd[0])!=0.0)
            ie=1;

         d=norm_square_dble(VOLUME,1,psd[1])/
            norm_square_dble(VOLUME,1,psd[3]);

         if (d>dmax)
            dmax=d;
      }

      error(ie,1,"main [check8.c]",
            "Dweo_blk_dble() changes the fields where it should not");

      dmax=sqrt(dmax);

      if (my_rank==0)
         printf("Dweo_blk_dble():                  %.1e\n",dmax);

      dmax=0.0;

      for (ic=0;ic<2;ic++)
      {
         random_sd(VOLUME,psd[0],1.0);
         random_sd(VOLUME,psd[1],1.0);
         random_sd(VOLUME,psd[2],1.0);
         blk_sd2zero(ic^0x1,psd[0]);
         blk_sd2zero(ic^0x1,psd[2]);

         if (ic^isw)
            n=nbh;
         else
            n=0;

         nm=n+nbh;

         for (;n<nm;n++)
         {
            assign_ud2udblk(DFL_BLOCKS,n);
            assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);

            assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
            assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[1],1);
            Dwoe_blk_dble(DFL_BLOCKS,n,0,1);
            assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,0,psd[2]);
            assign_sdblk2sd(DFL_BLOCKS,n,ALL_PTS,1,psd[3]);
         }

         Dwoe_dble(psd[0],psd[1]);
         blk_sd2zero(ic^0x1,psd[1]);
         blk_sd2zero(ic^0x1,psd[3]);

         mulr_spinor_add_dble(VOLUME,psd[0],psd[2],-1.0);
         mulr_spinor_add_dble(VOLUME,psd[1],psd[3],-1.0);

         if (norm_square_dble(VOLUME,0,psd[0])!=0.0)
            ie=1;

         d=norm_square_dble(VOLUME,1,psd[1])/
            norm_square_dble(VOLUME,1,psd[3]);

         if (d>dmax)
            dmax=d;
      }

      error(ie,1,"main [check8.c]",
            "Dwoe_blk_dble() changes the fields where it should not");

      dmax=sqrt(dmax);

      if (my_rank==0)
         printf("Dwoe_blk_dble():                  %.1e\n",dmax);

      dmax=0.0;
      random_sd(VOLUME,psd[0],1.0);
      random_sd(VOLUME,psd[1],1.0);

      for (n=0;n<nb;n++)
      {
         assign_ud2udblk(DFL_BLOCKS,n);
         assign_swd2swdblk(DFL_BLOCKS,n,NO_PTS);

         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
         Dwoe_blk_dble(DFL_BLOCKS,n,0,1);
         Dwee_blk_dble(DFL_BLOCKS,n,mu,0,0);
         Dwoo_blk_dble(DFL_BLOCKS,n,0.0,1,1);
         Dweo_blk_dble(DFL_BLOCKS,n,1,0);

         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],1);
         Dwhat_blk_dble(DFL_BLOCKS,n,mu,1,2);
         mulr_spinor_add_dble(volh,b[n].sd[0],b[n].sd[2],-1.0);
         d=norm_square_dble(volh,0,b[n].sd[0]);
         if (d>dmax)
            dmax=d;

         assign_sd2sdblk(DFL_BLOCKS,n,ALL_PTS,psd[0],0);
         mulr_spinor_add_dble(volh,b[n].sd[0]+volh,b[n].sd[1]+volh,-1.0);
         if (norm_square_dble(volh,0,b[n].sd[0]+volh)!=0.0)
            ie=1;
      }

      error(ie,1,"main [check8.c]",
            "Dwhat_blk_dble() changes the fields where it should not");

      dmax=sqrt(dmax);

      if (NPROC>1)
      {
         d=dmax;
         MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
         MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }

      if (my_rank==0)
         printf("Dwhat_blk_dble():                 %.1e\n\n",dmax);
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
