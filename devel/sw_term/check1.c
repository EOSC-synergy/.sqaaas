
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2011-2013, 2016 Martin Luescher
*               2017 Nazario Tantalo, Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation, assignment and inversion of the global SW arrays.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sw_term.h"
#include "global.h"
#include "gflds_utils.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

typedef union
{
   weyl_dble w;
   complex_dble c[6];
} spin_dble_t;

static pauli_dble *sswd=NULL;
static spin_dble_t vd ALIGNED32;
static const weyl_dble vd0={{{0.0}}};


static void save_swd(void)
{
   pauli_dble *pa,*pb,*pm;

   if (sswd==NULL)
   {
      sswd=amalloc(2*VOLUME*sizeof(*sswd),ALIGN);
      error(sswd==NULL,1,"save_swd [check1.c]",
            "Unable to allocate auxiliary array");
   }

   pa=swdfld();
   pb=sswd;
   pm=pa+2*VOLUME;

   for (;pa<pm;pa++)
   {
      (*pb)=(*pa);
      pb+=1;
   }
}


static int is_unity(pauli *p)
{
   int i,ie;
   float *u;

   ie=1;
   u=(*p).u;

   for (i=0;i<6;i++)
      ie|=(u[i]==1.0f);

   for (i=7;i<36;i++)
      ie|=(u[i]==0.0f);

   return ie;
}


static int is_unity_dble(pauli_dble *p)
{
   int i,ie;
   double *u;

   ie=1;
   u=(*p).u;

   for (i=0;i<6;i++)
      ie|=(u[i]==1.0);

   for (i=7;i<36;i++)
      ie|=(u[i]==0.0);

   return ie;
}


static int check_swbnd(void)
{
   int bc,ix,t,ie;
   pauli_dble *swd;

   bc=bc_type();
   swd=swdfld();
   ie=1;

   for (ix=0;ix<(2*VOLUME);ix++)
   {
      t=global_time(ix/2);

      if (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0)))
         ie|=is_unity_dble(swd);
      else
         ie|=(is_unity_dble(swd)^0x1);

      swd+=1;
   }

   return ie;
}


static double cmp_swd(ptset_t set)
{
   int k;
   double d,dmax;
   pauli_dble *pa,*pb,*pm;

   pa=swdfld();
   pb=sswd;
   pm=pa;

   if (set==EVEN_PTS)
      pm=pa+VOLUME;
   else if (set==ODD_PTS)
   {
      pa+=VOLUME;
      pb+=VOLUME;
      pm=pa+VOLUME;
   }
   else if (set==ALL_PTS)
      pm=pa+2*VOLUME;

   dmax=0.0;

   for (;pa<pm;pa++)
   {
      for (k=0;k<36;k++)
      {
         d=fabs((*pa).u[k]-(*pb).u[k]);

         if (d>dmax)
            dmax=d;
      }

      pb+=1;
   }

   return dmax;
}


static double cmp_iswd(ptset_t set)
{
   int k,l;
   double d,dmax;
   pauli_dble *pa,*pb,*pm;

   pa=swdfld();
   pb=sswd;
   pm=pa;

   if (set==EVEN_PTS)
      pm=pa+VOLUME;
   else if (set==ODD_PTS)
   {
      pa+=VOLUME;
      pb+=VOLUME;
      pm=pa+VOLUME;
   }
   else if (set==ALL_PTS)
      pm=pa+2*VOLUME;

   dmax=0.0;

   for (;pa<pm;pa++)
   {
      for (k=0;k<6;k++)
      {
         vd.w=vd0;
         vd.c[k].re=1.0;

         mul_pauli_dble(0.0,pa,&(vd.w),&(vd.w));
         mul_pauli_dble(0.0,pb,&(vd.w),&(vd.w));
         vd.c[k].re-=1.0;

         for (l=0;l<6;l++)
         {
            d=vd.c[l].re*vd.c[l].re+vd.c[l].im*vd.c[l].im;
            if (d>dmax)
               dmax=d;
         }
      }

      pb+=1;
   }

   return sqrt(dmax);
}


static double cmp_sw2swd(ptset_t set)
{
   int k;
   double d,dmax;
   pauli *pa,*pm;
   pauli_dble *pb;

   pa=swfld();
   pb=swdfld();
   pm=pa;

   if (set==EVEN_PTS)
      pm=pa+VOLUME;
   else if (set==ODD_PTS)
   {
      pa+=VOLUME;
      pb+=VOLUME;
      pm=pa+VOLUME;
   }
   else if (set==ALL_PTS)
      pm=pa+2*VOLUME;

   dmax=0.0;

   for (;pa<pm;pa++)
   {
      for (k=0;k<36;k++)
      {
         d=fabs((double)((*pa).u[k])-(*pb).u[k]);

         if (d>dmax)
            dmax=d;
      }

      pb+=1;
   }

   return dmax;
}


static void chs_ud0(void)
{
   int nlks,*lks,*lkm;
   su3_dble *ud,*vd;

   lks=bnd_lks(&nlks);
   lkm=lks+nlks;
   ud=udfld();

   for (;lks<lkm;lks++)
   {
      vd=ud+(*lks);

      (*vd).c11.re=-(*vd).c11.re;
      (*vd).c11.im=-(*vd).c11.im;
      (*vd).c12.re=-(*vd).c12.re;
      (*vd).c12.im=-(*vd).c12.im;
      (*vd).c13.re=-(*vd).c13.re;
      (*vd).c13.im=-(*vd).c13.im;

      (*vd).c21.re=-(*vd).c21.re;
      (*vd).c21.im=-(*vd).c21.im;
      (*vd).c22.re=-(*vd).c22.re;
      (*vd).c22.im=-(*vd).c22.im;
      (*vd).c23.re=-(*vd).c23.re;
      (*vd).c23.im=-(*vd).c23.im;

      (*vd).c31.re=-(*vd).c31.re;
      (*vd).c31.im=-(*vd).c31.im;
      (*vd).c32.re=-(*vd).c32.re;
      (*vd).c32.im=-(*vd).c32.im;
      (*vd).c33.re=-(*vd).c33.re;
      (*vd).c33.im=-(*vd).c33.im;
   }
}


static void change_ud_phase(void)
{
   int k,bc;
   double p;
   su3_dble *ud,*um;
   dirac_parms_t dp;
   complex_dble phase[3] ALIGNED16;   
   
   dp=dirac_parms();
   bc=bc_type();
      
   p=dp.theta[0]/(double)(N1);
   phase[0].re=cos(p);
   phase[0].im=sin(p);

   p=dp.theta[1]/(double)(N2);
   phase[1].re=cos(p);
   phase[1].im=sin(p);

   p=dp.theta[2]/(double)(N3);
   phase[2].re=cos(p);
   phase[2].im=sin(p);
   
   ud=udfld();
   um=ud+4*VOLUME;

   for (;ud<um;)
   {
      ud+=2;

      for (k=0;k<3;k++)
      {
         cm3x3_mulc(phase+k,ud,ud);
         ud+=1;
         cm3x3_mulc(phase+k,ud,ud);
         ud+=1;
      }
   }

   if ((cpr[0]==(NPROC0-1))&&((bc==1)||(bc==2)))
   {
      ud=udfld()+4*VOLUME+7*(BNDRY/4);
      cm3x3_mulc(phase,ud,ud);
      ud+=1;
      cm3x3_mulc(phase+1,ud,ud);
      ud+=1;
      cm3x3_mulc(phase+2,ud,ud);
   }

   if (bc==3)
      chs_ud0();

   set_flags(UPDATED_UD);   
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cf,ix,ie,q;
   double phi[2],phi_prime[2],csw[2],cF[2],theta[3];
   double d,dmax;
   pauli *sw;
   pauli_dble *swd;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      printf("\n");
      printf("Initialization and inversion of the global SW arrays\n");
      printf("----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-gg <gauge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check1.c]",
                  "Syntax: check1 [-bc <type>] [-gg <gauge>]");
      else
         cf=1;
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cf,1,MPI_INT,0,MPI_COMM_WORLD);
   set_flds_parms(cf,0);
   print_flds_parms();

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;
   set_bc_parms(bc,0,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,123456);
   geometry();
   
   q=-3;
   csw[0]=0.95;
   csw[1]=0.8;
   cF[0]=1.301;
   cF[1]=0.789;
   theta[0]=0.35;
   theta[1]=-1.25;
   theta[2]=0.78;
   if(gauge()==1)
   {
      q=0;
      csw[1]=0.0;
   }
   if(gauge()==2)
   {
      csw[0]=0.0;
   }
   if (bc_type()==3) cF[0]=cF[1]=0.0;
   if (bc_cstar()!=0)
   {
      theta[0]=0.0;
      theta[1]=0.0;
      theta[2]=0.0;
   }
   set_dirac_parms9(q,-0.0123,csw[0],csw[1],cF[0],cF[1],theta[0],theta[1],theta[2]);
   print_dirac_parms();

   sw=swfld();
   swd=swdfld();
   ie=1;

   for (ix=0;ix<(2*VOLUME);ix++)
   {
      ie|=is_unity(sw);
      ie|=is_unity_dble(swd);
      sw+=1;
      swd+=1;
   }

   error(ie!=1,1,"main [check1.c]","SW fields are not correctly initialized");

   random_gflds();
   print_flags();
   sw_term(NO_PTS);
   ie=check_swbnd();
   error(ie!=1,1,"main [check1.c]","SW field has incorrect boundary values");
   save_swd();

   change_ud_phase();
   sw_term(NO_PTS);
   d=cmp_swd(ALL_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("After changing the spatial boundary conditions\n");
      printf("Maximal deviation of swd = %.1e\n\n",dmax);
   }

   print_flags();
   ie=sw_term(EVEN_PTS);
   error(ie!=0,1,"main [check1.c]","Unsafe inversion of swd_e");
   d=cmp_iswd(EVEN_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Inverted swd_e\n");
      printf("Maximal deviation of swd_e = %.1e\n",dmax);
   }

   d=cmp_swd(ODD_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("Maximal deviation of swd_o = %.1e\n\n",dmax);

   print_flags();
   random_gflds();
   sw_term(NO_PTS);
   save_swd();

   ie=sw_term(ODD_PTS);
   error(ie!=0,1,"main [check1.c]","Unsafe inversion of swd_o");
   d=cmp_swd(EVEN_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Inverted swd_o\n");
      printf("Maximal deviation of swd_e = %.1e\n",dmax);
   }

   d=cmp_iswd(ODD_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
      printf("Maximal deviation of swd_o = %.1e\n\n",dmax);

   print_flags();
   assign_swd2sw();
   d=cmp_sw2swd(ALL_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Assigned swd to sw\n");
      printf("Maximal deviation = %.1e\n\n",dmax);
   }

   print_flags();
   random_gflds();
   sw_term(NO_PTS);
   save_swd();

   ie=sw_term(ALL_PTS);
   error(ie!=0,1,"main [check1.c]","Unsafe inversion of swd");
   d=cmp_iswd(ALL_PTS);
   MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Inverted swd\n");
      printf("Maximal deviation = %.1e\n\n",dmax);
   }

   print_flags();

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}

#undef N0
#undef N1
#undef N2
#undef N3
