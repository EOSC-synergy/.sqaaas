
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2007, 2010-2013, 2016 Martin Luescher
*               2017 Marina Marinkovic, Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Writing and reading gauge configurations.
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
#include "su3fcts.h"
#include "linalg.h"
#include "archive.h"
#include "global.h"

static int *rlxs_state[2],*rlxd_state[2];
static int ns,nd,*state=NULL;


static void alloc_state(void)
{
   int n;

   ns=rlxs_size();
   nd=rlxd_size();

   if (ns<nd)
      n=nd;
   else
      n=ns;

   state=malloc(n*sizeof(int));
   error(state==NULL,1,"alloc_state [u1archive.c]",
         "Unable to allocate auxiliary array");
}


static void save_ranlux(void)
{
   int nlxs,nlxd;
   int *p;

   nlxs=rlxs_size();
   nlxd=rlxd_size();

   p=malloc(2*(nlxs+nlxd)*sizeof(*p));
   error(p==NULL,1,"save_ranlux [check1.c]",
         "Unable to allocate state arrays");
   rlxs_state[0]=p;
   p+=nlxs;
   rlxs_state[1]=p;
   p+=nlxs;
   rlxd_state[0]=p;
   p+=nlxd;
   rlxd_state[1]=p;

   rlxs_get(rlxs_state[0]);
   rlxd_get(rlxd_state[0]);
}


static int check_ranlux(void)
{
   int nlxs,nlxd,k,ie;

   nlxs=rlxs_size();
   nlxd=rlxd_size();

   rlxs_get(rlxs_state[1]);
   rlxd_get(rlxd_state[1]);
   ie=0;

   for (k=0;k<nlxs;k++)
      ie|=(rlxs_state[0][k]!=rlxs_state[1][k]);

   for (k=0;k<nlxd;k++)
      ie|=(rlxd_state[0][k]!=rlxd_state[1][k]);

   return ie;
}


static int cmp_ud(su3_dble *u,su3_dble *v)
{
   int it;

   it =((*u).c11.re!=(*v).c11.re);
   it|=((*u).c11.im!=(*v).c11.im);
   it|=((*u).c12.re!=(*v).c12.re);
   it|=((*u).c12.im!=(*v).c12.im);
   it|=((*u).c13.re!=(*v).c13.re);
   it|=((*u).c13.im!=(*v).c13.im);

   it|=((*u).c21.re!=(*v).c21.re);
   it|=((*u).c21.im!=(*v).c21.im);
   it|=((*u).c22.re!=(*v).c22.re);
   it|=((*u).c22.im!=(*v).c22.im);
   it|=((*u).c23.re!=(*v).c23.re);
   it|=((*u).c23.im!=(*v).c23.im);

   it|=((*u).c31.re!=(*v).c31.re);
   it|=((*u).c31.im!=(*v).c31.im);
   it|=((*u).c32.re!=(*v).c32.re);
   it|=((*u).c32.im!=(*v).c32.im);
   it|=((*u).c33.re!=(*v).c33.re);
   it|=((*u).c33.im!=(*v).c33.im);

   return it;
}


static int check_ud(su3_dble *usv)
{
   int it;
   su3_dble *u,*um;

   u=udfld();
   um=u+4*VOLUME;
   it=0;

   for (;u<um;u++)
   {
      it|=cmp_ud(u,usv);
      usv+=1;
   }

   return it;
}


static int cmp_ad(double *u,double *v)
{
   return ((*u)!=(*v));
}


static int check_ad(double *asv)
{
   int it;
   double *a,*am;

   a=adfld();
   am=a+4*VOLUME;
   it=0;

   for (;a<am;a++)
   {
      it|=cmp_ad(a,asv);
      asv+=1;
   }

   return it;
}


static double check_cstar_ud(void)
{
   int size,tag;
   double d,dmax,dmax_all;
   su3_dble *rbuf,*udb;
   complex_dble *cp1,*cp2;
   MPI_Status stat;

   if(bc_cstar()==0) return 0.0;

   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();

   size=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc_type()==1)||(bc_type()==2)))
      size+=3;

   rbuf=malloc(size*sizeof(su3_dble));
   error(rbuf==NULL,1,"main [check1.c]",
         "Unable to allocate rbuf");

   udb=udfld();

   tag=mpi_tag();
   MPI_Sendrecv(udb,18*size,MPI_DOUBLE,get_mirror_rank(),tag,
                rbuf,18*size,MPI_DOUBLE,get_mirror_rank(),tag,
                MPI_COMM_WORLD,&stat);

   dmax=0.0;
   cp1=(complex_dble*)udb;
   cp2=(complex_dble*)rbuf;
   for (;cp1<(complex_dble*)udb+9*size;cp1++)
   {
      d=fabs((*cp1).re-(*cp2).re);
      if (d>dmax) dmax=d;
      d=fabs((*cp1).im+(*cp2).im);
      if (d>dmax) dmax=d;

      cp2++;
   }

   MPI_Allreduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   free(rbuf);

   return dmax_all;
}


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


static double loc_avg_su3_plaq(void)
{
   double plaq;
   plaq=plaq_sum_dble(0);
   return plaq/((double)(6)*(double)(VOLUME));
}


static double loc_avg_u1_plaq(void)
{
   double plaq;
   plaq=u1_plaq_sum_dble(0);
   return plaq/((double)(6)*(double)(VOLUME));
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,nsize,ie,gg,ir;
   double phi[2],phi_prime[2];
   su3_dble *udb,**usv;
   double *adb,**asv;
   char loc_dir[NAME_SIZE],cnfg[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;
   int cnfg_type;
   double dmax1_all,dmax2_all;
   double u1plaq,su3plaq;
   int ldat[16];

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Writing and reading gauge configurations of type 1/2/3 (QCD/QED/QED+QCD)\n");
      printf("----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("loc_dir","%s\n",loc_dir);
      fclose(fin);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      gg=find_opt(argc,argv,"-gg");

      if (gg!=0)
         error_root(sscanf(argv[gg+1],"%d",&gg)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");
      else
         gg=1;
   }

   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;

   set_flds_parms(gg,0);
   print_flds_parms();

   set_bc_parms(bc,0,cs,phi,phi_prime);
   print_bc_parms();

   start_ranlux(0,123456);
   geometry();

   check_dir(loc_dir);
   nsize=name_size("%s/testcnfg_%d",loc_dir,NPROC);
   error_root(nsize>=NAME_SIZE,1,"main [check1.c]","loc_dir name is too long");
   sprintf(cnfg,"%s/testcnfg_%d",loc_dir,my_rank);

   if (my_rank==0)
   {
      printf("Write random field configuration to the files\n"
             "%s/testcnfg_*\n"
             "on the local disks.\n\n",loc_dir);
      printf("Then read the field from there, compare with the saved field\n"
             "and remove all files.\n\n");
   }

   if ((gauge()&1)!=0)
   {
      alloc_wud(1);
      usv=reserve_wud(1);
      udb=udfld();
      random_ud();
      cm3x3_assign(4*VOLUME,udb,usv[0]);
   }
   if ((gauge()&2)!=0)
   {
      alloc_wad(1);
      asv=reserve_wad(1);
      adb=adfld();
      random_ad();
      assign_dvec2dvec(4*VOLUME,adb,asv[0]);
   }

   save_ranlux();
   write_cnfg(cnfg);

   if ((gauge()&1)!=0) random_ud();
   if ((gauge()&2)!=0) random_ad();

   cnfg_type=read_cnfg(cnfg);
   error(cnfg_type!=gauge(),1,"main [check1.c]","Incorrect configuration type");

   if (state==NULL)
      alloc_state();

   if (my_rank==0)
   {
      fin=fopen(cnfg,"rb");
      error_root(fin==NULL,1,"main [check1.c]","Unable to open input file");

      ir=fread(ldat,sizeof(int),16,fin);
      ir+=fread(state,sizeof(int),ns,fin);
      rlxs_reset(state);
      ir+=fread(state,sizeof(int),nd,fin);
      rlxd_reset(state);
      error_root(ir!=16+ns+nd,1,"main [check1.c]","Incorrect read count");

      if ((gauge()&1)!=0)
      {
         ir=fread(&su3plaq,sizeof(double),1,fin);
         error_root(ir!=1,1,"main [check1.c]","Incorrect read count");
         printf("SU(3) plaquette on node 0 from cnfg file  = %.15e\n",su3plaq/(double)(6*VOLUME));
         printf("Computed SU(3) plaquette on node 0        = %.15e\n\n",loc_avg_su3_plaq());
      }

      if ((gauge()&2)!=0)
      {
         if (gauge()==3)
         {
            ir=fread(udb,sizeof(su3_dble),4*VOLUME,fin);
            error_root(ir!=4*VOLUME,1,"main [check3.c]","Incorrect read count.");
         }
         ir=fread(&u1plaq,sizeof(double),1,fin);
         error_root(ir!=1,1,"main [check3.c]","Incorrect read count.");
         printf("U(1) plaquette on node 0 from cnfg file   = %.15e\n",u1plaq/(double)(6*VOLUME));
         printf("Computed U(1) plaquette on node 0         = %.15e\n\n",loc_avg_u1_plaq());
      }
      fclose(fin);
   }

   remove(cnfg);

   if ((gauge()&1)!=0)
   {
      ie=(check_bc(0.0)^0x1);
      ie|=check_ud(usv[0]);
      error(ie!=0,1,"main [check1.c]","The SU(3) gauge field is not properly restored");
   }

   if ((gauge()&2)!=0)
   {
      ie=(check_ad_bc(0.0)^0x1);
      ie|=check_ad(asv[0]);
      error(ie!=0,1,"main [check1.c]","The U(1) gauge field is not properly restored");
   }

   if(bc_cstar()!=0)
   {
      dmax1_all=check_cstar_ud();
      dmax2_all=check_cstar_ad();
      if (my_rank==0)
      {
         printf("C* boundary conditions ud = %.2e\n",dmax1_all);
         printf("C* boundary conditions ad = %.2e\n\n",dmax2_all);
      }
   }

   ie=check_ranlux();
   error(ie!=0,1,"main [check1.c]",
                 "The random number generator is not properly restored");

   print_flags();

   if (my_rank==0)
   {
      printf("No errors detected --- the fields are correctly written\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
