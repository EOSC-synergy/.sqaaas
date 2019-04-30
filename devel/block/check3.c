
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2005, 2011, 2013, 2016 Martin Luescher
*               2016, 2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of assign_ud2ubgr() and assign_ud2udblk().
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
#include "block.h"
#include "global.h"
#include "gflds_utils.h"

typedef union
{
   su3 u;
   float r[18];
} umat_t;

typedef union
{
   su3_dble u;
   double r[18];
} umat_dble_t;


static void set_ud(void)
{
   int x0,x1,x2,x3,ix;
   int y0,y1,y2,y3,ifc;
   su3_dble *udb,*ud;
   double *adb,*ad;

   random_gflds();
   udb=udfld();

   if ((gauge()&1)!=0)
   {
      for (x0=0;x0<L0;x0++)
      {
         for (x1=0;x1<L1;x1++)
         {
            for (x2=0;x2<L2;x2++)
            {
               for (x3=0;x3<L3;x3++)
               {
                  if ((x0+x1+x2+x3)&0x1)
                  {
                     y0=cpr[0]*L0+x0;
                     y1=cpr[1]*L1+x1;
                     y2=cpr[2]*L2+x2;
                     y3=cpr[3]*L3+x3;

                     ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
                     ud=udb+8*(ix-(VOLUME/2));

                     for (ifc=0;ifc<8;ifc++)
                     {
                        (*ud).c11.re=(double)(y0);
                        (*ud).c11.im=(double)(y1);

                        (*ud).c22.re=(double)(y2);
                        (*ud).c22.im=(double)(y3);

                        (*ud).c33.re=(double)(ifc);
                        (*ud).c33.im=0.0;

                        ud+=1;
                     }
                  }
               }
            }
         }
      }
   }

   if ((gauge()&2)!=0)
   {
      adb=adfld();
      
      for (x0=0;x0<L0;x0++)
      {
         for (x1=0;x1<L1;x1++)
         {
            for (x2=0;x2<L2;x2++)
            {
               for (x3=0;x3<L3;x3++)
               {
                  if ((x0+x1+x2+x3)&0x1)
                  {
                     y0=cpr[0]*L0+x0;
                     y1=cpr[1]*L1+x1;
                     y2=cpr[2]*L2+x2;
                     y3=cpr[3]*L3+x3;

                     ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
                     ad=adb+8*(ix-(VOLUME/2));

                     for (ifc=0;ifc<8;ifc++)
                     {
                        (*ad)=(double)(y0+2*y1+3*y2+4*y3+5*ifc);
                        ad+=1;
                     }
                  }
               }
            }
         }
      }
   }

   set_bc();
   set_ad_bc();

   set_flags(UPDATED_UD);
   set_flags(UPDATED_AD);
}


static double check_u(int y0,int y1,int y2,int y3,int ifc,su3 *u)
{
   double d,dmax;
   int sign=1;
   complex_dble v,w,z;
   int q;
   
   if(bc_type()==3&&((y0==0&&ifc==1)||(y0==NPROC0*L0-1&&ifc==0)))
      sign=-1;
   
   if(bc_type()==0&&((y0==0&&ifc==1)||(y0==NPROC0*L0-1&&ifc==0)))
      sign=0;
   
   if(bc_type()==1&&(y0==0&&ifc>1))
      return 0.0;
   
   q=dirac_parms().qhat;

   dmax=0.0;
   
   if(gauge()==1)
   {
      d=fabs((*u).c11.re-(float)(sign*y0));
      if(d>dmax) dmax=d;
      d=fabs((*u).c11.im-(float)(sign*y1));
      if(d>dmax) dmax=d;
      d=fabs((*u).c22.re-(float)(sign*y2));
      if(d>dmax) dmax=d;
      d=fabs((*u).c22.im-(float)(sign*y3));
      if(d>dmax) dmax=d;
      d=fabs((*u).c33.re-(float)(sign*ifc));
      if(d>dmax) dmax=d;
      d=fabs((*u).c33.im-0.0f);
      if(d>dmax) dmax=d;
   }
   else if(gauge()==3)
   {
      w.re=cos(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));
      w.im=sin(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));

      #define cmpl_mul(a,b,c) \
         (a).re=(b).re*(c).re-(b).im*(c).im; \
         (a).im=(b).re*(c).im+(b).im*(c).re

      z.re=(double)(sign*y0);
      z.im=(double)(sign*y1);
      cmpl_mul(v,w,z);
      d=fabs((*u).c11.re-(float)(v.re));
      if(d>dmax) dmax=d;
      d=fabs((*u).c11.im-(float)(v.im));
      if(d>dmax) dmax=d;

      z.re=(double)(sign*y2);
      z.im=(double)(sign*y3);
      cmpl_mul(v,w,z);
      d=fabs((*u).c22.re-(float)(v.re));
      if(d>dmax) dmax=d;
      d=fabs((*u).c22.im-(float)(v.im));
      if(d>dmax) dmax=d;
      
      z.re=(double)(sign*ifc);
      z.im=0.;
      cmpl_mul(v,w,z);
      d=fabs((*u).c33.re-(float)(v.re));
      if(d>dmax) dmax=d;
      d=fabs((*u).c33.im-(float)(v.im));
      if(d>dmax) dmax=d;
   }
   else
   {
      v.re=sign*cos(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));
      v.im=sign*sin(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));

      d=fabs((*u).c11.re-(float)(v.re));
      if(d>dmax) dmax=d;
      d=fabs((*u).c11.im-(float)(v.im));
      if(d>dmax) dmax=d;
      d=fabs((*u).c22.re-(float)(v.re));
      if(d>dmax) dmax=d;
      d=fabs((*u).c22.im-(float)(v.im));
      if(d>dmax) dmax=d;
      d=fabs((*u).c33.re-(float)(v.re));
      if(d>dmax) dmax=d;
      d=fabs((*u).c33.im-(float)(v.im));
      if(d>dmax) dmax=d;
   }

   return dmax;
}


static double check_ublk(block_t *b)
{
   double d,dmax;
   int *bo,*bs;
   int x0,x1,x2,x3,ix;
   int y0,y1,y2,y3,ifc;
   su3 *ub,*u;

   bo=(*b).bo;
   bs=(*b).bs;
   ub=(*b).u;
   dmax=0.0;

   for (x0=0;x0<bs[0];x0++)
   {
      for (x1=0;x1<bs[1];x1++)
      {
         for (x2=0;x2<bs[2];x2++)
         {
            for (x3=0;x3<bs[3];x3++)
            {
               if ((x0+x1+x2+x3)&0x1)
               {
                  y0=cpr[0]*L0+bo[0]+x0;
                  y1=cpr[1]*L1+bo[1]+x1;
                  y2=cpr[2]*L2+bo[2]+x2;
                  y3=cpr[3]*L3+bo[3]+x3;

                  ix=(*b).ipt[x3+bs[3]*x2+bs[2]*bs[3]*x1+bs[1]*bs[2]*bs[3]*x0];
                  u=ub+8*(ix-((*b).vol/2));

                  for (ifc=0;ifc<8;ifc++)
                  {
                     d=check_u(y0,y1,y2,y3,ifc,u);
                     if(d>dmax) dmax=d;
                     u+=1;
                  }
               }
            }
         }
      }
   }
   
   if(NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static double check_ud(int y0,int y1,int y2,int y3,int ifc,su3_dble *ud)
{
   double d,dmax;
   int sign=1;
   complex_dble v,w,z;
   int q;
   
   if(bc_type()==3&&((y0==0&&ifc==1)||(y0==NPROC0*L0-1&&ifc==0)))
      sign=-1;
   
   if(bc_type()==0&&((y0==0&&ifc==1)||(y0==NPROC0*L0-1&&ifc==0)))
      sign=0;
   
   if(bc_type()==1&&(y0==0&&ifc>1))
      return 0.0;
   
   q=dirac_parms().qhat;

   dmax=0.0;
   
   if(gauge()==1)
   {
      d=fabs((*ud).c11.re-(double)(sign*y0));
      if(d>dmax) dmax=d;
      d=fabs((*ud).c11.im-(double)(sign*y1));
      if(d>dmax) dmax=d;
      d=fabs((*ud).c22.re-(double)(sign*y2));
      if(d>dmax) dmax=d;
      d=fabs((*ud).c22.im-(double)(sign*y3));
      if(d>dmax) dmax=d;
      d=fabs((*ud).c33.re-(double)(sign*ifc));
      if(d>dmax) dmax=d;
      d=fabs((*ud).c33.im-0.0f);
      if(d>dmax) dmax=d;
   }
   else if(gauge()==3)
   {
      w.re=cos(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));
      w.im=sin(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));

      #define cmpl_mul(a,b,c) \
         (a).re=(b).re*(c).re-(b).im*(c).im; \
         (a).im=(b).re*(c).im+(b).im*(c).re

      z.re=(double)(sign*y0);
      z.im=(double)(sign*y1);
      cmpl_mul(v,w,z);
      d=fabs((*ud).c11.re-v.re);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c11.im-v.im);
      if(d>dmax) dmax=d;

      z.re=(double)(sign*y2);
      z.im=(double)(sign*y3);
      cmpl_mul(v,w,z);
      d=fabs((*ud).c22.re-v.re);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c22.im-v.im);
      if(d>dmax) dmax=d;
      
      z.re=(double)(sign*ifc);
      z.im=0.;
      cmpl_mul(v,w,z);
      d=fabs((*ud).c33.re-v.re);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c33.im-v.im);
      if(d>dmax) dmax=d;
   }
   else
   {
      v.re=sign*cos(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));
      v.im=sign*sin(q*(double)(y0+2*y1+3*y2+4*y3+5*ifc));

      d=fabs((*ud).c11.re-v.re);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c11.im-v.im);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c22.re-v.re);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c22.im-v.im);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c33.re-v.re);
      if(d>dmax) dmax=d;
      d=fabs((*ud).c33.im-v.im);
      if(d>dmax) dmax=d;
   }

   return dmax;
}


static double check_udblk(block_t *b)
{
   double d,dmax;
   int *bo,*bs;
   int x0,x1,x2,x3,ix;
   int y0,y1,y2,y3,ifc;
   su3_dble *udb,*ud;

   bo=(*b).bo;
   bs=(*b).bs;
   udb=(*b).ud;
   dmax=0.0;

   for (x0=0;x0<bs[0];x0++)
   {
      for (x1=0;x1<bs[1];x1++)
      {
         for (x2=0;x2<bs[2];x2++)
         {
            for (x3=0;x3<bs[3];x3++)
            {
               if ((x0+x1+x2+x3)&0x1)
               {
                  y0=cpr[0]*L0+bo[0]+x0;
                  y1=cpr[1]*L1+bo[1]+x1;
                  y2=cpr[2]*L2+bo[2]+x2;
                  y3=cpr[3]*L3+bo[3]+x3;

                  ix=(*b).ipt[x3+bs[3]*x2+bs[2]*bs[3]*x1+bs[1]*bs[2]*bs[3]*x0];
                  ud=udb+8*(ix-((*b).vol/2));

                  for (ifc=0;ifc<8;ifc++)
                  {
                     d=check_ud(y0,y1,y2,y3,ifc,ud);
                     if(d>dmax) dmax=d;
                     ud+=1;
                  }
               }
            }
         }
      }
   }

   if(NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static void fnd_coord(block_t *b,int ix,int *x0,int *x1,int *x2,int *x3)
{
   int iz;

   (*x0)=0;
   iz=(*b).idn[ix][0];

   while (iz<(*b).vol)
   {
      (*x0)+=1;
      iz=(*b).idn[iz][0];
   }

   (*x1)=0;
   iz=(*b).idn[ix][1];

   while (iz<(*b).vol)
   {
      (*x1)+=1;
      iz=(*b).idn[iz][1];
   }

   (*x2)=0;
   iz=(*b).idn[ix][2];

   while (iz<(*b).vol)
   {
      (*x2)+=1;
      iz=(*b).idn[iz][2];
   }

   (*x3)=0;
   iz=(*b).idn[ix][3];

   while (iz<(*b).vol)
   {
      (*x3)+=1;
      iz=(*b).idn[iz][3];
   }
}


static double cmp_u(su3 *u,su3 *v)
{
   double d,dmax;
   int i;
   umat_t *uu,*uv;

   uu=(umat_t*)(u);
   uv=(umat_t*)(v);
   
   dmax=0.0;

   for (i=0;i<18;i++)
   {
      d=fabs((*uu).r[i]-(*uv).r[i]);
      if(d>dmax) dmax=d;
   }

   return dmax;
}


static double is_zero(su3 *u)
{
   double d,dmax;
   int i;
   umat_t *um;

   um=(umat_t*)(u);
   dmax=0.0;

   for (i=0;i<18;i++)
   {
      d=fabs((*um).r[i]);
      if(d>dmax) dmax=d;
   }

   return dmax;
}


static int is_on_bnd(int ix)
{
   int bc,t;

   bc=bc_type();
   t=global_time(ix);

   if (((t==0)&&(bc!=3))||((t==(NPROC0*L0-1))&&(bc==0)))
      return 1;
   else
      return 0;
}


static double check_ubnd(block_t *b)
{
   double d,dmax;
   int ix,iz,ifc,*bo;
   int x0,x1,x2,x3;
   int y0,y1,y2,y3;
   su3 *u;
   bndry_t *bb;

   bo=(*b).bo;
   bb=(*b).bb;
   dmax=0.0;

   if ((*bb).u==NULL)
      return 0.0;

   for (ifc=0;ifc<8;ifc++)
   {
      u=(*bb).u;

      for (iz=0;iz<(*bb).vol;iz++)
      {
         ix=(*bb).ipp[iz];

         if ((ifc<=1)&&(is_on_bnd((*b).imb[ix])))
         {
            d=is_zero(u);
            if(d>dmax) dmax=d;
         }
         else
         {
            if (iz<((*bb).vol/2))
            {
               d=cmp_u(u,(*b).u+8*(ix-((*b).vol/2))+(ifc^0x1));
               if(d>dmax) dmax=d;
            }

            fnd_coord(b,ix,&x0,&x1,&x2,&x3);

            y0=cpr[0]*L0+bo[0]+x0;
            y1=cpr[1]*L1+bo[1]+x1;
            y2=cpr[2]*L2+bo[2]+x2;
            y3=cpr[3]*L3+bo[3]+x3;

            if (iz<((*bb).vol/2))
            {
               d=check_u(y0,y1,y2,y3,ifc^0x1,u);
               if(d>dmax) dmax=d;
            }
            else
            {
               if (ifc==0)
                  y0=safe_mod(y0-1,NPROC0*L0);
               if (ifc==1)
                  y0=safe_mod(y0+1,NPROC0*L0);

               if (ifc==2)
                  y1=safe_mod(y1-1,NPROC1*L1);
               if (ifc==3)
                  y1=safe_mod(y1+1,NPROC1*L1);

               if (ifc==4)
                  y2=safe_mod(y2-1,NPROC2*L2);
               if (ifc==5)
                  y2=safe_mod(y2+1,NPROC2*L2);

               if (ifc==6)
                  y3=safe_mod(y3-1,NPROC3*L3);
               if (ifc==7)
                  y3=safe_mod(y3+1,NPROC3*L3);

               d=check_u(y0,y1,y2,y3,ifc,u);
               if(d>dmax) dmax=d;
            }
         }

         u+=1;
      }

      bb+=1;
   }
   
   if(NPROC>1)
   {
      d=dmax;
      MPI_Reduce(&d,&dmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      MPI_Bcast(&dmax,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   return dmax;
}


static double check_ubgr(blk_grid_t grid)
{
   double d,dmax;
   int nb,isw,n;
   block_t *b;

   b=blk_list(grid,&nb,&isw);
   dmax=0;

   for (n=0;n<nb;n++)
   {
      d=check_ublk(b+n);
      if(d>dmax) dmax=d;
      d=check_ubnd(b+n);
      if(d>dmax) dmax=d;
   }

   return dmax;
}


int main(int argc,char *argv[])
{
   double d,dmax;
   int my_rank,bc,bs[4],cf,qhat;
   int nb,isw,n;
   double phi[2],phi_prime[2];
   double cF[2];
   block_t *b;
   FILE *flog=NULL,*fin=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
      fin=freopen("check1.in","r",stdin);

      printf("\n");
      printf("Check of assign_ud2ubgr() and assign_ud2udblk()\n");
      printf("-----------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      read_line("bs","%d %d %d %d",&bs[0],&bs[1],&bs[2],&bs[3]);
      fclose(fin);

      printf("bs = %d %d %d %d\n\n",bs[0],bs[1],bs[2],bs[3]);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>] [-gg <gauge>] [-q <echarge>]");

      cf=find_opt(argc,argv,"-gg");

      if (cf!=0)
         error_root(sscanf(argv[cf+1],"%d",&cf)!=1,1,"main [check3.c]",
                  "Syntax: check3 [-bc <type>] [-gg <gauge>] [-q <echarge>]");
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
   if (gauge()==2)
   {
      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
   }   
   set_bc_parms(bc,0,0,phi,phi_prime);
   print_bc_parms();
   
   qhat=0;
   cF[0]=cF[1]=0.0;
   if ((gauge()&2)!=0) qhat=-3;
   if (bc_type()!=3) cF[0]=cF[1]=1.0;
   set_dirac_parms9(qhat,DBL_MAX,0.0,0.0,cF[0],cF[1],0.0,0.0,0.0);
   print_dirac_parms();

   start_ranlux(0,1234);
   geometry();
   set_sap_parms(bs,0,1,1);
   set_dfl_parms(bs,2);
   alloc_bgr(SAP_BLOCKS);
   alloc_bgr(DFL_BLOCKS);

   set_ud();
   assign_ud2ubgr(SAP_BLOCKS);
   print_flags();
   print_grid_flags(SAP_BLOCKS);

   dmax=check_ubgr(SAP_BLOCKS);
   if (my_rank==0)
   {
      printf("Deviation in blocks after assign_ud2ubgr(SAP_BLOCKS) = %.2e\n\n",dmax);
   }

   b=blk_list(DFL_BLOCKS,&nb,&isw);
   random_gflds();
   assign_ud2udblk(DFL_BLOCKS,0);
   set_ud();
   dmax=0.0;

   for (n=0;n<nb;n++)
   {
      assign_ud2udblk(DFL_BLOCKS,n);
      d=check_udblk(b+n);
      if(d>dmax) dmax=d;
   }

   print_flags();
   print_grid_flags(DFL_BLOCKS);

   if (my_rank==0)
   {
      printf("Deviation in blocks after assign_ud2udblk(DFL_BLOCKS,*) = %.2e\n\n",dmax);
   }

   if (my_rank==0)
   {
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
