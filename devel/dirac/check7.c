
/*******************************************************************************
*
* File check7.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Comparison of Dw_blk(),..,Dwhat_blk() with Dw().
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
#include "u1flds.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "block.h"
#include "dirac.h"
#include "global.h"
#include "gflds_utils.h"


static void blk_s2zero(int ic,spinor *s)
{
   int nb,isw;
   int nbh,n,nm,vol;
   block_t *b;

   b=blk_list(SAP_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;

   if (ic^isw)
      n=nbh;
   else
      n=0;

   nm=n+nbh;

   for (;n<nm;n++)
   {
      set_s2zero(vol,b[n].s[0]);
      assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,s);
   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cf,q;
   int nb,isw,nbh,ic,itm;
   int bs[4],n,nm,vol,volh,ie;
   float mu,d,dmax;
   double phi[2],phi_prime[2];
   double su3csw,u1csw,cF[2],theta[3];
   spinor **ps;
   block_t *b;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check7.log","w",stdout);
      fin=freopen("check7.in","r",stdin);

      printf("\n");
      printf("Comparison of Dw_blk(),..,Dwhat_blk() with Dw()\n");
      printf("-----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check7.c]",
                    "Syntax: check7 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check7.c]",
                  "Syntax: check7 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
      else
         cf=1;

      q=find_opt(argc,argv,"-q");

      if (q!=0)
      {
         error_root(sscanf(argv[q+1],"%d",&q)!=1,1,"main [check7.c]",
                  "Syntax: check7 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
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
   set_sap_parms(bs,0,1,1);
   alloc_bgr(SAP_BLOCKS);
   alloc_ws(4);

   mu=0.123f;
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

   assign_swd2sw();
   assign_ud2ubgr(SAP_BLOCKS);
   assign_swd2swbgr(SAP_BLOCKS,NO_PTS);

   ps=reserve_ws(4);
   b=blk_list(SAP_BLOCKS,&nb,&isw);
   nbh=nb/2;
   vol=(*b).vol;
   volh=vol/2;

   for (itm=0;itm<2;itm++)
   {
      ie=0;
      dmax=0.0f;
      set_tm_parms(itm);

      if (my_rank==0)
         printf("Twisted-mass flag = %d\n",itm);

      for (ic=0;ic<2;ic++)
      {
         random_s(VOLUME,ps[0],1.0f);
         random_s(VOLUME,ps[2],1.0f);
         blk_s2zero(ic^0x1,ps[0]);
         blk_s2zero(ic^0x1,ps[2]);

         if (ic^isw)
            n=nbh;
         else
            n=0;

         nm=n+nbh;

         for (;n<nm;n++)
         {
            random_s(vol,b[n].s[1],1.0f);
            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
            Dw_blk(SAP_BLOCKS,n,mu,0,1);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,ps[2]);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,1,ps[3]);
         }

         Dw(mu,ps[0],ps[1]);
         blk_s2zero(ic^0x1,ps[1]);
         blk_s2zero(ic^0x1,ps[3]);

         mulr_spinor_add(VOLUME,ps[0],ps[2],-1.0f);
         mulr_spinor_add(VOLUME,ps[1],ps[3],-1.0f);

         if (norm_square(VOLUME,0,ps[0])!=0.0f)
            ie=1;

         d=norm_square(VOLUME,1,ps[1])/
            norm_square(VOLUME,1,ps[3]);

         if (d>dmax)
            dmax=d;
      }

      error(ie,1,"main [check7.c]",
            "Dw_blk() changes the fields where it should not");

      dmax=(float)(sqrt((double)(dmax)));

      if (my_rank==0)
      {
         printf("The maximal relative deviations are:\n\n");
         printf("Dw_blk():               %.1e\n",dmax);
      }

      dmax=0.0f;
      random_s(VOLUME,ps[0],1.0f);
      random_s(VOLUME,ps[1],1.0f);

      for (n=0;n<nb;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[1],1);
         Dwee_blk(SAP_BLOCKS,n,mu,0,1);
         assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,ps[2]);
         assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,1,ps[3]);

         assign_s2sblk(SAP_BLOCKS,n,EVEN_PTS,ps[0],0);
         assign_s2sblk(SAP_BLOCKS,n,ODD_PTS,ps[1],0);
         Dwee_blk(SAP_BLOCKS,n,mu,0,0);
         mulr_spinor_add(vol,b[n].s[0],b[n].s[1],-1.0f);
         if (norm_square(vol,0,b[n].s[0])!=0.0f)
            ie=1;
      }

      Dwee(mu,ps[0],ps[1]);
      mulr_spinor_add(VOLUME,ps[0],ps[2],-1.0f);
      mulr_spinor_add(VOLUME,ps[1],ps[3],-1.0f);

      if (norm_square(VOLUME,0,ps[0])!=0.0f)
         ie=1;

      d=norm_square(VOLUME,1,ps[1])/
         norm_square(VOLUME,1,ps[3]);

      if (d>dmax)
         dmax=d;

      random_s(VOLUME,ps[0],1.0f);
      random_s(VOLUME,ps[1],1.0f);

      for (n=0;n<nb;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[1],1);
         Dwoo_blk(SAP_BLOCKS,n,mu,0,1);
         assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,ps[2]);
         assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,1,ps[3]);

         assign_s2sblk(SAP_BLOCKS,n,ODD_PTS,ps[0],0);
         assign_s2sblk(SAP_BLOCKS,n,EVEN_PTS,ps[1],0);
         Dwoo_blk(SAP_BLOCKS,n,mu,0,0);
         mulr_spinor_add(vol,b[n].s[0],b[n].s[1],-1.0f);
         if (norm_square(vol,0,b[n].s[0])!=0.0f)
            ie=1;
      }

      Dwoo(mu,ps[0],ps[1]);
      mulr_spinor_add(VOLUME,ps[0],ps[2],-1.0f);
      mulr_spinor_add(VOLUME,ps[1],ps[3],-1.0f);

      if (norm_square(VOLUME,0,ps[0])!=0.0f)
         ie=1;

      d=norm_square(VOLUME,1,ps[1])/
         norm_square(VOLUME,1,ps[3]);

      if (d>dmax)
         dmax=d;

      error(ie,1,"main [check7.c]",
            "Dwee_blk() or Dwoo_blk() changes the fields where it should not");

      dmax=(float)(sqrt((double)(dmax)));

      if (my_rank==0)
         printf("Dwee_blk(), Dwoo_blk(): %.1e\n",dmax);

      dmax=0.0f;

      for (ic=0;ic<2;ic++)
      {
         random_s(VOLUME,ps[0],1.0f);
         random_s(VOLUME,ps[1],1.0f);
         random_s(VOLUME,ps[2],1.0f);
         blk_s2zero(ic^0x1,ps[0]);
         blk_s2zero(ic^0x1,ps[2]);

         if (ic^isw)
            n=nbh;
         else
            n=0;

         nm=n+nbh;

         for (;n<nm;n++)
         {
            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[1],1);
            Dweo_blk(SAP_BLOCKS,n,0,1);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,ps[2]);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,1,ps[3]);
         }

         Dweo(ps[0],ps[1]);
         blk_s2zero(ic^0x1,ps[1]);
         blk_s2zero(ic^0x1,ps[3]);

         mulr_spinor_add(VOLUME,ps[0],ps[2],-1.0f);
         mulr_spinor_add(VOLUME,ps[1],ps[3],-1.0f);

         if (norm_square(VOLUME,0,ps[0])!=0.0f)
            ie=1;

         d=norm_square(VOLUME,1,ps[1])/
            norm_square(VOLUME,1,ps[3]);

         if (d>dmax)
            dmax=d;
      }

      error(ie,1,"main [check7.c]",
            "Dweo_blk() changes the fields where it should not");

      dmax=(float)(sqrt((double)(dmax)));

      if (my_rank==0)
         printf("Dweo_blk():             %.1e\n",dmax);

      dmax=0.0f;

      for (ic=0;ic<2;ic++)
      {
         random_s(VOLUME,ps[0],1.0f);
         random_s(VOLUME,ps[1],1.0f);
         random_s(VOLUME,ps[2],1.0f);
         blk_s2zero(ic^0x1,ps[0]);
         blk_s2zero(ic^0x1,ps[2]);

         if (ic^isw)
            n=nbh;
         else
            n=0;

         nm=n+nbh;

         for (;n<nm;n++)
         {
            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
            assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[1],1);
            Dwoe_blk(SAP_BLOCKS,n,0,1);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,0,ps[2]);
            assign_sblk2s(SAP_BLOCKS,n,ALL_PTS,1,ps[3]);
         }

         Dwoe(ps[0],ps[1]);
         blk_s2zero(ic^0x1,ps[1]);
         blk_s2zero(ic^0x1,ps[3]);

         mulr_spinor_add(VOLUME,ps[0],ps[2],-1.0f);
         mulr_spinor_add(VOLUME,ps[1],ps[3],-1.0f);

         if (norm_square(VOLUME,0,ps[0])!=0.0f)
            ie=1;

         d=norm_square(VOLUME,1,ps[1])/
            norm_square(VOLUME,1,ps[3]);

         if (d>dmax)
            dmax=d;
      }

      error(ie,1,"main [check7.c]",
            "Dwoe_blk() changes the fields where it should not");

      dmax=(float)(sqrt((double)(dmax)));

      if (my_rank==0)
         printf("Dwoe_blk():             %.1e\n",dmax);

      dmax=0.0f;
      random_s(VOLUME,ps[0],1.0f);
      random_s(VOLUME,ps[1],1.0f);

      for (n=0;n<nb;n++)
      {
         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
         Dwoe_blk(SAP_BLOCKS,n,0,1);
         Dwee_blk(SAP_BLOCKS,n,mu,0,0);
         Dwoo_blk(SAP_BLOCKS,n,0.0f,1,1);
         Dweo_blk(SAP_BLOCKS,n,1,0);

         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],1);
         Dwhat_blk(SAP_BLOCKS,n,mu,1,2);
         mulr_spinor_add(volh,b[n].s[0],b[n].s[2],-1.0f);
         d=norm_square(volh,0,b[n].s[0]);
         if (d>dmax)
            dmax=d;

         assign_s2sblk(SAP_BLOCKS,n,ALL_PTS,ps[0],0);
         mulr_spinor_add(volh,b[n].s[0]+volh,b[n].s[1]+volh,-1.0f);
         if (norm_square(volh,0,b[n].s[0]+volh)!=0.0f)
            ie=1;
      }

      error(ie,1,"main [check7.c]",
            "Dwhat_blk() changes the fields where it should not");

      dmax=(float)(sqrt((double)(dmax)));

      if (NPROC>1)
      {
         d=dmax;
         MPI_Reduce(&d,&dmax,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
         MPI_Bcast(&dmax,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      }

      if (my_rank==0)
         printf("Dwhat_blk():            %.1e\n\n",dmax);
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
