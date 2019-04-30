
/*******************************************************************************
*
* File check4.c
*
* Copyright (C) 2017 Nazario Tantalo, Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the program u1mom_Delta_no0 with plane waves 
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "utils.h"
#include "random.h"
#include "flags.h"
#include "lattice.h"
#include "linalg.h"
#include "flags.h"
#include "dft.h"
#include "u1flds.h"
#include "mdflds.h"
#include "global.h"

static int bc,cs;
static int b[4],n[4];
static int idp[4],rcpr[4];
static double phat[4],p[4][4];
static double phat0[4],p0[4][4];
static double *inmom,*outmom,*expmom;
static complex_dble *mom[4],*pmom[4],*spmom[4],*ispmom[4],*ipmom[4];


static void alloc_flds(void)
{
   int i;

   mom[0]=amalloc(20*VOLUME*sizeof(*(mom[0])),ALIGN);
   error(mom[0]==NULL,1,"alloc_flds [check4.c]",
         "Unable to allocate memory space for service arrays");
   mom[1]=mom[0]+VOLUME;
   mom[2]=mom[1]+VOLUME;
   mom[3]=mom[2]+VOLUME;

   pmom[0]=mom[3]+VOLUME;
   pmom[1]=pmom[0]+VOLUME;
   pmom[2]=pmom[1]+VOLUME;
   pmom[3]=pmom[2]+VOLUME;

   ipmom[0]=pmom[3]+VOLUME;
   ipmom[1]=ipmom[0]+VOLUME;
   ipmom[2]=ipmom[1]+VOLUME;
   ipmom[3]=ipmom[2]+VOLUME;

   spmom[0]=ipmom[3]+VOLUME;
   spmom[1]=spmom[0]+VOLUME;
   spmom[2]=spmom[1]+VOLUME;
   spmom[3]=spmom[2]+VOLUME;

   ispmom[0]=spmom[3]+VOLUME;
   ispmom[1]=ispmom[0]+VOLUME;
   ispmom[2]=ispmom[1]+VOLUME;
   ispmom[3]=ispmom[2]+VOLUME;

   i=4*VOLUME+7*(BNDRY/4);
   inmom=amalloc(3*i*sizeof(*inmom),ALIGN);
   error(inmom==NULL,1,"alloc_flds [check4.c]",
         "Unable to allocate memory space for service arrays");   
   outmom=inmom+i;
   expmom=outmom+i;
}


static void set_idp(void)
{
   int mu,id[4],nproc[4],nx[4];
   
   b[0]=0;
   b[1]=(cs>0);
   b[2]=(cs>1);
   b[3]=(cs>2);
   if(bc==2)
      b[0]=1;
   
   nproc[0]=NPROC0;
   nproc[1]=NPROC1;
   nproc[2]=NPROC2;
   nproc[3]=NPROC3;
   if(b[1]==1)
      nproc[1]/=2;
   
   nx[0]=L0;
   nx[1]=L1;
   nx[2]=L2;
   nx[3]=L3;
   for (mu=0;mu<4;++mu)
   {
      rcpr[mu]=cpr[mu]%nproc[mu];
      n[mu]=nproc[mu]*nx[mu];
   }
   if(bc==0)
      n[0]-=1;   
   if ( (rcpr[0]==NPROC0-1)&&(bc>0)&&(bc<3) )
      nx[0]+=1;
   
   if (bc==0)
   {
      id[0]=set_dft_parms(SIN,n[0],b[0],1);
      id[1]=set_dft_parms(EXP,n[1],b[1],0);
      id[2]=set_dft_parms(EXP,n[2],b[2],0);
      id[3]=set_dft_parms(EXP,n[3],b[3],0);
      
      idp[0]=set_dft4d_parms(id,nx,1);
      
      for (mu=1;mu<4;mu++)
      {
         id[0]=set_dft_parms(COS,n[0],b[0],0);
         id[1]=set_dft_parms(EXP,n[1],b[1],(mu==1));
         id[2]=set_dft_parms(EXP,n[2],b[2],(mu==2));
         id[3]=set_dft_parms(EXP,n[3],b[3],(mu==3));
    
         idp[mu]=set_dft4d_parms(id,nx,1);
      }
   }
   else if (bc==1)
   {
      id[0]=set_dft_parms(COS,n[0],b[0],1);
      id[1]=set_dft_parms(EXP,n[1],b[1],0);
      id[2]=set_dft_parms(EXP,n[2],b[2],0);
      id[3]=set_dft_parms(EXP,n[3],b[3],0);
      
      idp[0]=set_dft4d_parms(id,nx,1);
      
      for (mu=1;mu<4;mu++)
      {
         id[0]=set_dft_parms(SIN,n[0],b[0],0);
         id[1]=set_dft_parms(EXP,n[1],b[1],(mu==1));
         id[2]=set_dft_parms(EXP,n[2],b[2],(mu==2));
         id[3]=set_dft_parms(EXP,n[3],b[3],(mu==3));
    
         idp[mu]=set_dft4d_parms(id,nx,1);
      }
   }
   else if (bc==2)
   {
      id[0]=set_dft_parms(SIN,n[0],b[0],1);
      id[1]=set_dft_parms(EXP,n[1],b[1],0);
      id[2]=set_dft_parms(EXP,n[2],b[2],0);
      id[3]=set_dft_parms(EXP,n[3],b[3],0);
      
      idp[0]=set_dft4d_parms(id,nx,1);
      
      for (mu=1;mu<4;mu++)
      {
         id[0]=set_dft_parms(COS,n[0],b[0],0);
         id[1]=set_dft_parms(EXP,n[1],b[1],(mu==1));
         id[2]=set_dft_parms(EXP,n[2],b[2],(mu==2));
         id[3]=set_dft_parms(EXP,n[3],b[3],(mu==3));
    
         idp[mu]=set_dft4d_parms(id,nx,1);
      }
   }
   else
      for (mu=0;mu<4;mu++)
      {
         id[0]=set_dft_parms(EXP,n[0],b[0],(mu==0));
         id[1]=set_dft_parms(EXP,n[1],b[1],(mu==1));
         id[2]=set_dft_parms(EXP,n[2],b[2],(mu==2));
         id[3]=set_dft_parms(EXP,n[3],b[3],(mu==3));
    
         idp[mu]=set_dft4d_parms(id,nx,1);
      }
}


static void set_phat(void)
{
   int mu,a,k;
   float ran[4];
   double pi,r0,r1;
   dft4d_parms_t *dp;

   pi=4.0*atan(1.0);
      
   for (a=0;a<4;a++)
   {
      dp=dft4d_parms(idp[a]);
      
      ranlxs(ran,4);

      phat[a]=0.0;
      for (mu=0;mu<4;mu++)
      {
         k=(int)(ran[mu]*(float)(n[mu]));
         if ((*dp).dp[mu]->type==EXP)
            r0=2.0*pi/(double)(n[mu]);
         else
            r0=pi/(double)(n[mu]);

         if (2*k+b[mu]!=0)
         {
            p[a][mu]=0.5*r0*(double)(2*k+b[mu]);
            r1=2.0*sin(0.25*r0*(double)(2*k+b[mu]));
            phat[a]+=r1*r1;
         }
         else
            p[a][mu]=0.0;
      }
      if (phat[a]==0.0)
         phat[a]=1.0;

      MPI_Bcast(p[a],4,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   MPI_Bcast(phat,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


static void set_moms(void)
{
   int x0,x1,x2,x3,ix,a,mu;   
   double w[4],px,x[4];
   dft4d_parms_t *dp;
   
   gauss_dble(w,4);
   MPI_Bcast(w,4,MPI_DOUBLE,0,MPI_COMM_WORLD);

   for (x0=0;x0<L0;x0++)
   for (x1=0;x1<L1;x1++)
   for (x2=0;x2<L2;x2++)
   for (x3=0;x3<L3;x3++)
   {
      ix=x3+L3*x2+L2*L3*x1+L1*L2*L3*x0;

      for (a=0;a<4;a++)
      {
         mom[a][ix].im=0.0;
         pmom[a][ix].im=0.0;
         ipmom[a][ix].im=0.0;

         dp=dft4d_parms(idp[a]);

         x[0]=(double)(2*x0+2*rcpr[0]*L0+(*dp).dp[0]->c);
         x[1]=(double)(2*x1+2*rcpr[1]*L1+(*dp).dp[1]->c);
         x[2]=(double)(2*x2+2*rcpr[2]*L2+(*dp).dp[2]->c);
         x[3]=(double)(2*x3+2*rcpr[3]*L3+(*dp).dp[3]->c);

         mom[a][ix].re=w[a];
         for (mu=0;mu<4;mu++)
         {
            px= 0.5*p[a][mu]*x[mu];        

            if ((*dp).dp[mu]->type==SIN)
               mom[a][ix].re*=sin(px);
            else
               mom[a][ix].re*=cos(px);
         }
         pmom[a][ix].re=mom[a][ix].re*phat[a];
         ipmom[a][ix].re=mom[a][ix].re/phat[a];
         spmom[a][ix].re=mom[a][ix].re/sqrt(phat[a]);
         ispmom[a][ix].re=mom[a][ix].re*sqrt(phat[a]);
      }
   }
}


static void set_phat0(void)
{
   int mu,a,k;
   double pi,r0,r1;
   dft4d_parms_t *dp;

   pi=4.0*atan(1.0);
      
   for (a=0;a<4;a++)
   {
      dp=dft4d_parms(idp[a]);

      k=0;
      phat0[a]=0.0;
      for (mu=0;mu<4;mu++)
      {
         if ((*dp).dp[mu]->type==EXP)
            r0=2.0*pi/(double)(n[mu]);
         else
            r0=pi/(double)(n[mu]);

         if (2*k+b[mu]!=0)
         {
            p0[a][mu]=0.5*r0*(double)(2*k+b[mu]);
            r1=2.0*sin(0.25*r0*(double)(2*k+b[mu]));
            phat0[a]+=r1*r1;
         }
         else
            p0[a][mu]=0.0;
      }
      if (phat0[a]==0.0)
         phat0[a]=1.0;

      MPI_Bcast(p0[a],4,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   MPI_Bcast(phat0,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


static void add_moms0(void)
{
   int x0,x1,x2,x3,ix,a,mu;   
   double w[4],mm,px,x[4];
   dft4d_parms_t *dp;
   
   gauss_dble(w,4);
   MPI_Bcast(w,4,MPI_DOUBLE,0,MPI_COMM_WORLD);

   for (x0=0;x0<L0;x0++)
   for (x1=0;x1<L1;x1++)
   for (x2=0;x2<L2;x2++)
   for (x3=0;x3<L3;x3++)
   {
      ix=x3+L3*x2+L2*L3*x1+L1*L2*L3*x0;

      for (a=0;a<4;a++)
      {
         dp=dft4d_parms(idp[a]);

         x[0]=(double)(2*x0+2*rcpr[0]*L0+(*dp).dp[0]->c);
         x[1]=(double)(2*x1+2*rcpr[1]*L1+(*dp).dp[1]->c);
         x[2]=(double)(2*x2+2*rcpr[2]*L2+(*dp).dp[2]->c);
         x[3]=(double)(2*x3+2*rcpr[3]*L3+(*dp).dp[3]->c);

         mm=w[a];
         for (mu=0;mu<4;mu++)
         {
            px= 0.5*p0[a][mu]*x[mu];        

            if ((*dp).dp[mu]->type==SIN)
               mm*=sin(px);
            else
               mm*=cos(px);
         }
         mom[a][ix].re+=mm;
         pmom[a][ix].re+=mm*phat0[a];
         ipmom[a][ix].re+=mm/phat0[a];
         spmom[a][ix].re+=mm/sqrt(phat0[a]);
         ispmom[a][ix].re+=mm*sqrt(phat0[a]);
      }
   }
}


static void Gop(double *mdf,int flag)
{
   int ix,ifc;
   double w;
   
   if(flag==0) w=sqrt(0.5);
   else w=sqrt(2.0);
   
   if(((bc==0)||(bc==2))&&(cpr[0]==0))
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==0)
         {
            for(ifc=2;ifc<8;ifc++)
               mdf[8*(ix-VOLUME/2)+ifc]*=w;
         }
      }
   }
   else if((bc==0)&&(cpr[0]==NPROC0-1))
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==NPROC0*L0-1)
         {
            for(ifc=2;ifc<8;ifc++)
               mdf[8*(ix-VOLUME/2)+ifc]*=w;
         }
      }
   }
}


static void bnd_u1mom(double *mom)
{
   int bc,ifc;
   int nlks,*lks,*lkm,npts,*pts,*ptm;
   double *m;

   bc=bc_type();

   if ((bc==0)||(bc==1))
   {
      if (bc==0)
      {
         lks=bnd_lks(&nlks);
         lkm=lks+nlks;

         for (;lks<lkm;lks++)
            mom[*lks]=0.0;
      }
      else if (bc==1)
      {
         pts=bnd_pts(&npts);
         ptm=pts+npts;
         pts+=(npts/2);

         for (;pts<ptm;pts++)
         {
            m=mom+8*(pts[0]-(VOLUME/2));

            for (ifc=2;ifc<8;ifc++)
               m[ifc]=0.0;
         }
      }
   }
}


int main(int argc,char *argv[])
{
   int it,mu,my_rank;
   double d0,d1,d2,d3;
   double phi[2],phi_prime[2],u1phi,u1phi_prime;
   double d0max,d1max,d2max,d3max;
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check4.log","w",stdout);

      printf("\n");
      printf("Check of the program u1mom_Delta_no0 with plane waves\n");
      printf("-----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
      
      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-cs <cstar>]");

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check4.c]",
                    "Syntax: check4 [-bc <type>] [-cs <cstar>]");
   }
   
   set_flds_parms(2,0);
   print_flds_parms();

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   u1phi=0.0;
   u1phi_prime=0.0;
   set_bc_parms(bc,cs,phi,phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   start_ranlux(0,1236);
   geometry();

   alloc_flds();
   set_idp();

   d0max=0.0;
   d1max=0.0;
   d2max=0.0;
   d3max=0.0;

   for(it=0;it<20;++it)
   {
      set_phat();
      set_moms();
      
      set_phat0();
      add_moms0();
      
      scatter_u1mom(mom,expmom);
      bnd_u1mom(expmom);
      
      scatter_u1mom(ipmom,inmom);
      bnd_u1mom(inmom);
      Gop(inmom,0);
      u1mom_Delta_no0(0,inmom,outmom);
      Gop(outmom,1);
      
      muladd_assign_dvec(4*VOLUME,-1.0,expmom,outmom);
      d0=norm_square_dvec(4*VOLUME,1,outmom)/norm_square_dvec(4*VOLUME,1,inmom);
      d0=sqrt(d0);
      if (d0>d0max)
         d0max=d0;
      
      scatter_u1mom(pmom,inmom);
      bnd_u1mom(inmom);
      Gop(inmom,0);
      u1mom_Delta_no0(1,inmom,outmom);
      Gop(outmom,1);
      
      muladd_assign_dvec(4*VOLUME,-1.0,expmom,outmom);
      d1=norm_square_dvec(4*VOLUME,1,outmom)/norm_square_dvec(4*VOLUME,1,inmom);
      d1=sqrt(d1);
      if (d1>d1max)
         d1max=d1;
      
      scatter_u1mom(spmom,inmom);
      bnd_u1mom(inmom);
      Gop(inmom,0);
      u1mom_Delta_no0(2,inmom,outmom);
      Gop(outmom,1);
      
      muladd_assign_dvec(4*VOLUME,-1.0,expmom,outmom);
      d2=norm_square_dvec(4*VOLUME,1,outmom)/norm_square_dvec(4*VOLUME,1,inmom);
      d2=sqrt(d2);
      if (d2>d2max)
         d2max=d2;
      
      scatter_u1mom(ispmom,inmom);
      bnd_u1mom(inmom);
      Gop(inmom,0);
      u1mom_Delta_no0(3,inmom,outmom);
      Gop(outmom,1);
      
      muladd_assign_dvec(4*VOLUME,-1.0,expmom,outmom);
      d3=norm_square_dvec(4*VOLUME,1,outmom)/norm_square_dvec(4*VOLUME,1,inmom);
      d3=sqrt(d3);
      if (d3>d3max)
         d3max=d3;
      
      if (my_rank==0)
      {
         printf("\n\nPlane wave momenta\n");
         for (mu=0;mu<4;++mu)
            printf("Pi_%d, p={ %.2e, %.2e, %.2e, %.2e }\n",mu,p[mu][0],p[mu][1],p[mu][2],p[mu][3]);

         printf("Normalized deviations: d0= %.1e, d1= %.1e, d2= %.1e, d3= %.1e\n",d0,d1,d2,d3);
      }
   }

   if (my_rank==0)
   {
      printf("\n\nMaximanl deviations: d0= %.1e, d1= %.1e, d2= %.1e, d3= %.1e\n",d0max,d1max,d2max,d3max);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
