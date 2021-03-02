
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2007, 2008, 2010-2013, 2016 Martin Luescher
*               2017 Marina Marinkovic, Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Importing a previously exported configuration.
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
#include "su3fcts.h"
#include "linalg.h"
#include "archive.h"
#include "u1flds.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static su3_dble *vbuf,*ubuf;


static void alloc_ubuf(int my_rank)
{
   if (my_rank==0)
   {
      ubuf=amalloc(4*(L3+N3)*sizeof(su3_dble),ALIGN);
      vbuf=ubuf+4*L3;
   }
   else
   {
      ubuf=amalloc(4*L3*sizeof(su3_dble),ALIGN);
      vbuf=NULL;
   }

   error(ubuf==NULL,1,"alloc_ubuf [check3.c]",
         "Unable to allocate auxiliary array for SU(3) gauge fields");
}


static void cstar_su3_dble(int vol,su3_dble* u)
{
  su3_dble *um; 

  um=u+vol;
  for (;u<um;u++)
  {    
     (*u).c11.im *= -1.0;
     (*u).c12.im *= -1.0;
     (*u).c13.im *= -1.0;
     (*u).c21.im *= -1.0;
     (*u).c22.im *= -1.0;
     (*u).c23.im *= -1.0;
     (*u).c31.im *= -1.0;
     (*u).c32.im *= -1.0;
     (*u).c33.im *= -1.0;
  }
}


static void cstar_double(int vol,double* a)
{
  double *am; 

  am=a+vol;
  for (;a<am;a++)  (*a) *= -1.0;
}


static double avg_su3_plaq(void)
{
   double plaq;
   
   if ((bc_cstar()!=0))
   {
      if (cpr[1]>=NPROC1/2)
         cstar_su3_dble(4*VOLUME,udfld());
      set_flags(UPDATED_UD);
   }

   plaq=plaq_sum_dble(1);

   if ((bc_cstar()!=0))
   {
      if (cpr[1]>=NPROC1/2)
         cstar_su3_dble(4*VOLUME,udfld());
      set_flags(UPDATED_UD);
   }
   
   return plaq/((double)(6*NPROC)*(double)(VOLUME));
}


static double avg_u1_plaq(void)
{
   double plaq;

   if ((bc_cstar()!=0))
   {
      if (cpr[1]>=NPROC1/2)
         cstar_double(4*VOLUME,adfld());
      set_flags(UPDATED_AD);
   }

   plaq=u1_plaq_sum_dble(1);

   if ((bc_cstar()!=0))
   {
      if (cpr[1]>=NPROC1/2)
         cstar_double(4*VOLUME,adfld());
      set_flags(UPDATED_AD);
   }

   return plaq/((double)(6*NPROC)*(double)(VOLUME));
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
   error(rbuf==NULL,1,"main [check3.c]",
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
   error(rbuf==NULL,1,"main [check3.c]",
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
   int my_rank,bc,nsize,ie,cnfg_type,cs,gg,itest;
   long long int ir=0,n0,n1,n2,n3,ix,globvol;
   stdint_t l[4];
   double phi[2],phi_prime[2],eps;
   double su3plaq0=0.0,su3plaq1=0.0,su3plaq2=0.0,u1plaq0=0.0,u1plaq1=0.0,u1plaq2=0.0;
   char cnfg_dir[NAME_SIZE],cnfg[NAME_SIZE];
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check3.in","r",stdin);

      printf("\n");
      printf("Importing gauge fields exported by check3\n");
      printf("-----------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("cnfg_dir","%s\n",cnfg_dir);
      fclose(fin);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");

      gg=find_opt(argc,argv,"-gg");

      if (gg!=0)
         error_root(sscanf(argv[gg+1],"%d",&gg)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>] [-gg <gauge>]");
      else
         gg=1;
   }

   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);

   phi[0]=0.123;
   phi[1]=-0.534;
   phi_prime[0]=0.912;
   phi_prime[1]=0.078;

   set_flds_parms(gg,0);
   print_flds_parms();

   set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
   print_bc_parms();

   start_ranlux(0,9876);
   geometry();
   
   eps=sqrt((double)(6*N0*N1)*(double)(N2*N3))*DBL_EPSILON;


   if (vbuf==NULL)
      alloc_ubuf(my_rank);

   su3plaq0=u1plaq0=0.0;
   if ((gauge()&1)!=0)
   {
      udfld();
      random_ud();
      su3plaq0=avg_su3_plaq();
   }
   if ((gauge()&2)!=0)
   {
      adfld();
      random_ad();
      u1plaq0=avg_u1_plaq();
   }

   check_dir_root(cnfg_dir);
   nsize=name_size("%s/testcnfg",cnfg_dir);
   error_root(nsize>=NAME_SIZE,1,"main [check3.c]","cnfg_dir name is too long");
   sprintf(cnfg,"%s/testcnfg",cnfg_dir);

   cnfg_type=import_cnfg(cnfg);

   itest=0;
   if (my_rank==0)
   {
      fin=fopen(cnfg,"rb");
      error_root(fin==NULL,1,"main [check3.c]","Unable to open input file");

      ir=fread(l,sizeof(stdint_t),4,fin);
      error_root(ir!=4,1,"main [check3.c]","Incorrect read count");
      if (endianness()==BIG_ENDIAN)
         bswap_int(4,l);

      if ((gauge()&1)!=0)
         printf("Random gauge field, average SU(3) plaquette = %.15e\n",su3plaq0);
      if ((gauge()&2)!=0)
         printf("Random gauge field, average U(1) plaquette  = %.15e\n\n",u1plaq0);
      printf("Now read ");
      if (cnfg_type==1)
         printf("SU(3)");
      else if(cnfg_type==2)
         printf("U(1)");
      else
         printf("SU(3)xU(1)");
      printf(" gauge field from file\n"
             "%s:\n",cnfg);
      printf("%dx%dx%dx%d lattice\n\n",
             (int)(l[0]),(int)(l[1]),(int)(l[2]),(int)(l[3]));
   }

   if ((gauge()&1)!=0)
   {
      if (my_rank==0)
      {
         if ((cnfg_type&1)!=0)
         {
            ir+=fread(&su3plaq1,sizeof(double),1,fin);
            error_root(ir!=5,1,"main [check3.c]","Incorrect read count");
            if (endianness()==BIG_ENDIAN)
               bswap_double(1,&su3plaq1);
         }
         else
            su3plaq1=su3plaq0;
      }

      ie=check_bc(0.0);
      error(ie!=1,1,"main [check3.c]","Boundary conditions of the SU(3) field are not preserved");
      su3plaq2=avg_su3_plaq();
   
      if (bc_cstar()!=0)
      {
         error(check_cstar_ud()!=0.0,1,"main [check2.c]",
               "The SU(3) gauge field does not satisfy C* boundary conditions");
      }
   
      if (my_rank==0)
      {
         printf("Calculated SU(3) plaquette = %.15e\n",su3plaq2);
         printf("SU(3) plaquette should be  = %.15e\n\n",su3plaq1);
         
         itest|=(su3plaq2-su3plaq1)>eps;
      }
   }

   if ((gauge()&2)!=0)
   {
      if (my_rank==0)
      {
         if ((cnfg_type&2)!=0)
         {
            if (cnfg_type==3)
            {
               n0=l[0];
               n1=l[1];
               n2=l[2];
               n3=l[3];

               ir=0;
               for (ix=0;ix<(n0*n1*n2);ix++)
                  ir+=fread(vbuf,sizeof(su3_dble),4*n3,fin);
               globvol=n0*n1*n2*n3;

               error_root(ir!=4*globvol,1,"main [check3.c]","Incorrect read count.");
            }
            ir=fread(&u1plaq1,sizeof(double),1,fin);
            error_root(ir!=1,1,"main [check3.c]","Incorrect read count.");

            if (endianness()==BIG_ENDIAN)
               bswap_double(1,&u1plaq1);
         }
         else
            u1plaq1=u1plaq0;
      }

      ie=check_ad_bc(0.0);
      error(ie!=1,1,"main [check3.c]","Boundary conditions of the U(1) field are not preserved");
      u1plaq2=avg_u1_plaq();

      if (bc_cstar()!=0)
      {
         error(check_cstar_ad()!=0.0,1,"main [check2.c]",
               "The U(1) gauge field does not satisfy C* boundary conditions");
      }
      
      if (my_rank==0)
      {
         printf("Calculated U(1) plaquette = %.15e\n",u1plaq2);
         printf("U(1) plaquette should be  = %.15e\n\n",u1plaq1);
         
         itest|=(u1plaq2-u1plaq1)>eps;
      }
   }
   
   if (my_rank==0)
      fclose(fin);

   print_flags();
   
   error_root(itest==1,1,"main [check3.c]","Something went wrong.");

   if (my_rank==0)
   {
      printf("No errors detected --- the fields are correctly exported\n\n");
      remove(cnfg);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
