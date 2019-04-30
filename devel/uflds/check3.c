
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program that translates the double-precision gauge field.
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
#include "global.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)
#define GLB_VOL (N0*N1*N2*N3)


/* given a point ix and a direction mu,
      udfld()+ret[ix][mu]
   returns the gauge link in (ix,mu)
*/
int** lidx(void) {
   int ix;
   int ip[4];
   int **ptr;

   ptr=malloc(VOLUME*sizeof(int*));
   ptr[0]=malloc(4*VOLUME*sizeof(int));
   for(ix=1;ix<VOLUME;ix++) ptr[ix]=ptr[ix-1]+4;
   
   /* Memo plaquette
    *  n      ={0    ,1    ,2    ,3    ,4    ,5    } 
    *  (mu,nu)={(0,1),(0,2),(0,3),(2,3),(3,1),(1,2)}
    *
    *   ip[0] -> U(x,mu)
    *   ip[1] -> U(x+mu,nu)         
    *   ip[2] -> U(x,nu)
    *   ip[3] -> U(x+nu,mu)
    */
   for(ix=0;ix<VOLUME;ix++) {
      plaq_uidx(0,ix,ip);
      ptr[ix][0]=ip[0];
      ptr[ix][1]=ip[2];
      plaq_uidx(3,ix,ip);
      ptr[ix][2]=ip[0];
      ptr[ix][3]=ip[2];
   }
   
   return ptr;
}


static double gf[2][GLB_VOL];
void save_ud_cmp(int mu,int c)
{
   static int **ln=NULL;
   su3_dble *ud;
   int ix,gx,x0,x1,x2,x3,gx0,gx1,gx2,gx3;
   
   if (ln==NULL) ln=lidx();
   
   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();
   
   ud=udfld();
   
   for (ix=0;ix<GLB_VOL;ix++)
      gf[0][ix]=0.0;

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
               
               gx0=cpr[0]*L0+x0;
               gx1=cpr[1]*L1+x1;
               gx2=cpr[2]*L2+x2;
               gx3=cpr[3]*L3+x3;
               
               gx=gx3+N3*gx2+N2*N3*gx1+N1*N2*N3*gx0;
               
               gf[0][gx]=((double*)(ud+ln[ix][mu]))[c];
            }
         }
      }
   }
   
   MPI_Allreduce(gf[0],gf[1],GLB_VOL,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

double dev_shifted_ud_cmp(int mu,int c,int s[4])
{
   static int **ln=NULL;
   su3_dble *ud;
   int ix,gx,x0,x1,x2,x3,gx0,gx1,gx2,gx3;
   double d,dmax;
   
   if (ln==NULL) ln=lidx();
   
   if (query_flags(UDBUF_UP2DATE)!=1)
      copy_bnd_ud();
   
   ud=udfld();
   
   dmax=0.0;
   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
               
               gx0=cpr[0]*L0+x0-s[0];
               gx1=cpr[1]*L1+x1-s[1];
               gx2=cpr[2]*L2+x2-s[2];
               gx3=cpr[3]*L3+x3-s[3];
               
               gx0=safe_mod(gx0,N0);
               gx2=safe_mod(gx2,2*N2);
               if (gx2>=N2)
               {
                  gx2-=N2;
                  if(bc_cstar()>=2) gx1+=N1/2;
               }
               gx3=safe_mod(gx3,2*N3);
               if (gx3>=N3)
               {
                  gx3-=N3;
                  if(bc_cstar()>=3) gx1+=N1/2;
               }
               gx1=safe_mod(gx1,N1);
               
               gx=gx3+N3*gx2+N2*N3*gx1+N1*N2*N3*gx0;
               
               d=fabs(gf[1][gx]-((double*)(ud+ln[ix][mu]))[c]);
               if (d>dmax) dmax=d;
            }
         }
      }
   }
   
   d=dmax;
   MPI_Allreduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   
   return dmax;
}


static double check_cstar(void)
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


static void random_vec(int *svec)
{
   int bs[4];
   double r[4];

   bs[0]=NPROC0*L0;
   bs[1]=NPROC1*L1;
   bs[2]=NPROC2*L2;
   bs[3]=NPROC3*L3;

   ranlxd(r,4);

   svec[0]=(int)((r[0]-0.5)*(3*bs[0]));
   svec[1]=(int)((r[1]-0.5)*(3*bs[1]));
   svec[2]=(int)((r[2]-0.5)*(3*bs[2]));
   svec[3]=(int)((r[3]-0.5)*(3*bs[3]));
   
   MPI_Bcast(svec,4,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int my_rank,bc,cs,ie,itest;
   int ifc,mu,c,n,s[4];
   double phi[2],phi_prime[2],dev;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);

      printf("\n");
      printf("Translation of the double-precision gauge field\n");
      printf("-----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-cs <cstar>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   if(cs==0)
   {
      phi[0]=0.123;
      phi[1]=-0.534;
      phi_prime[0]=0.912;
      phi_prime[1]=0.078;
   }
   set_bc_parms(bc,cs,phi,phi_prime,0.0,0.0);
   print_bc_parms();

   start_ranlux(0,4325);
   geometry();

   if (my_rank==0)
      printf("Elementary shift vectors:\n\n");

   for (ifc=0;ifc<8;ifc++)
   {
      if ((ifc>1)||(bc==3))
      {
         s[0]=0;
         s[1]=0;
         s[2]=0;
         s[3]=0;

         if ((ifc&0x1)==0)
            s[ifc/2]=1;
         else
            s[ifc/2]=-1;

         itest=0;
         for (mu=0;mu<4;mu++)
         {
            for (c=0;c<18;c++)
            {
               random_ud();
               save_ud_cmp(mu,c);
               shift_ud(s);
               dev=dev_shifted_ud_cmp(mu,c,s);
               if (dev!=0.0) itest=1;

               ie=check_bc(0.0);
               error_root(ie==0,1,"main [check3.c]","Boundary conditions changed");
            }
         }

         if (my_rank==0)
         {
            printf("Shift vector (% 4d,% 4d,% 4d,% 4d): ",
                   s[0],s[1],s[2],s[3]);

            if (itest==0)
               printf("ok\n");
            else
               printf("failed\n");
         }
      }
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Random shift vectors:\n\n");
   }

   for (n=0;n<8;n++)
   {
      random_vec(s);
      if (bc!=3)
         s[0]=0;

      itest=0;
      for (mu=0;mu<4;mu++)
      {
         for (c=0;c<18;c++)
         {
            random_ud();
            save_ud_cmp(mu,c);
            shift_ud(s);
            dev=dev_shifted_ud_cmp(mu,c,s);
            if (dev!=0.0) itest=1;

            ie=check_bc(0.0);
            error_root(ie==0,1,"main [check3.c]","Boundary conditions changed");
         }
      }

      if (my_rank==0)
      {
         printf("Shift vector (% 4d,% 4d,% 4d,% 4d): ",
                s[0],s[1],s[2],s[3]);

         if (itest==0)
            printf("ok\n");
         else
            printf("failed\n");
      }
   }
   
   if(bc_cstar()!=0)
   {
      dev=check_cstar();
      if (my_rank==0)
         printf("\nCheck C* boundary conditions after shift = %.2e\n",dev);
   }

   if (my_rank==0)
   {
      printf("\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
