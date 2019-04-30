
/*******************************************************************************
*
* File u1mom_facc.c
*
* Copyright (C) 2017 Nazario Tantalo, Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge-field propagator and generation of random gauge potentials.
*
* The externally accessible functions are
*
*   void u1mom_Delta_no0(int flag,double* inmom,double* outmom)
*     The double array inmom is assumed to contain a vector field IN_mu(x)
*     stored as the U(1) gauge potentials and to include the corresponding 
*     communication buffers. Similarly, the double array outmom will contain
*     the vector field OUT_mu(x) where 
*       OUT_mu(x)= G^{1/2} Delta_no0 G^{-1/2} IN_mu(x)           if flag=0,
*       OUT_mu(x)= G^{1/2} Delta_no0^{-1} G^{-1/2} IN_mu(x)      if flag=1,
*       OUT_mu(x)= G^{1/2} Delta_no0^{1/2} G^{-1/2} IN_mu(x)     if flag=2,
*       OUT_mu(x)= G^{1/2} Delta_no0^{-1/2} G^{-1/2} IN_mu(x)    if flag=3.
*       OUT_mu(x)= IN_mu(x)                           otherwise.
*     The operators Delta_no0 and G are described in the notes below.
*     The pointers inmom and outmom may be equal.
*
* Notes:
*
* The operator Delta_no0 is diagonal in Fourier space. The allowed momenta
* appearing in the Fourier transform of IN_mu(x),
*
* FIN_mu(p_0,p_1,p_2,p_3)= FT{ IN_mu(x) } ,
*
* depend upon the boundary conditions (see doc/fourier.pdf). For example,
*
* bc=0, mu=0, bc_cstar()=1, 
*
* p_0=   pi*n0/(NPROC0*L0-1),   n0=1,...,NPROC0*L0-1,
* p_1= 4*pi*(n1+1/2)/NPROC1*L1, n1=0,...,NPROC1*L1/2-1,
* p_2= 2*pi*n2/NPROC2*L2,       n2=0,...,NPROC2*L2-1,
* p_3= 2*pi*n3/NPROC3*L3,       n3=0,...,NPROC3*L3-1,
*
* or,
*
* bc=1, mu=0, bc_cstar()=0, 
*
* p_0=   pi*n0/NPROC0*L0,       n0=0,...,NPROC0*L0-1,
* p_1= 2*pi*n1/NPROC1*L1,       n1=0,...,NPROC1*L1-1,
* p_2= 2*pi*n2/NPROC2*L2,       n2=0,...,NPROC2*L2-1,
* p_3= 2*pi*n3/NPROC3*L3,       n3=0,...,NPROC3*L3-1,
*
* The operator Delta_no0 is such that
*
* Delta_no0(p_0,p_1,p_2,p_3) * FIN_mu(p_0,p_1,p_2,p_3) =
*                                     = FT{ Delta_no0 IN_mu(x) } ,
*
* phat^2= sum_{mu}( 2*sin(p_mu/2) )^2 ,
*
* Delta_no0(p_0,p_1,p_2,p_3)= 1,      phat=0,
* Delta_no0(p_0,p_1,p_2,p_3)= phat^2, phat!=0.
*
* The operator G is diagonal in coordinate space. It leaves all link untouched,
* except the spatial links on open boundaries, which are multiplied times 1/2.
* The operator G^{1/2} Delta_no0 G^{-1/2} is symmetric and positive-definite
* with respect to the canonical scalar product.
*
* In case of open boundary conditions, the non-active links of outmom are set
* to zero.
* 
* The size of the arrays inmom and outmom must be at least 4*VOLUME+7*(BNDRY/4).
*
* The programs in this module act globally and must be called on all MPI
* processes simultaneously.
*
*******************************************************************************/

#define U1MOM_FACC_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "flags.h"
#include "lattice.h"
#include "mdflds.h"
#include "dft.h"
#include "u1flds.h"
#include "global.h"


static int b[4],n[4],nx[4];
static int vol,idp[4],rcpr[4];
static double *phatsq;
static complex_dble *mom[4];


static void alloc_mom(void)
{
   int bc,vol,i;
   
   bc=bc_type();
   
   if ( (cpr[0]==NPROC0-1)&&(bc>0)&&(bc<3) )
      vol=VOLUME+L1*L2*L3;
   else
      vol=VOLUME;
   
   mom[0]=amalloc(4*vol*sizeof(*mom[0]),ALIGN);
   mom[1]=mom[0]+vol;
   mom[2]=mom[1]+vol;
   mom[3]=mom[2]+vol;
   
   error(mom[0]==NULL,1,"alloc_mom [u1mom_facc.c]",
         "Unable to allocate field array");
   
   for (i=0;i<(4*vol);i++)
   {
      mom[0][i].re=0.0;
      mom[0][i].im=0.0;
   }
}


static void set_idp(void)
{
   int bc,mu,cs;
   int id[4],nproc[4];
   
   error(iup[0][0]==0,1,"set_idp [u1mom_facc.c]",
         "Geometry arrays are not set");
   
   bc=bc_type();
   cs=bc_cstar();
   
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

   vol=nx[0]*nx[1]*nx[2]*nx[3];
   
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
   {
      for (mu=0;mu<4;mu++)
      {
         id[0]=set_dft_parms(EXP,n[0],b[0],(mu==0));
         id[1]=set_dft_parms(EXP,n[1],b[1],(mu==1));
         id[2]=set_dft_parms(EXP,n[2],b[2],(mu==2));
         id[3]=set_dft_parms(EXP,n[3],b[3],(mu==3));

         idp[mu]=set_dft4d_parms(id,nx,1);
      }
   }
}


static void set_phat(void)
{
   int mu,k;
   int k0,k1,k2,k3;
   double pi,r0,r1,*phat[4];
   double p0,p1,p2,p3;
   dft4d_parms_t *dp;
   
   dp=dft4d_parms(idp[0]);

   k=n[0]+n[1]+n[2]+n[3];
   phat[0]=malloc(k*sizeof(*(phat[0])));
   error(phat[0]==NULL,1,"set_phat [u1mom_facc.c]",
         "Unable to allocate auxiliary array");
   phat[1]=phat[0]+n[0];
   phat[2]=phat[1]+n[1];
   phat[3]=phat[2]+n[2];
   
   phatsq=malloc(vol*sizeof(*phatsq));
   error(phatsq==NULL,1,"set_phat [u1mom_facc.c]",
         "Unable to allocate auxiliary array");

   pi=4.0*atan(1.0);
   
   for (mu=0;mu<4;mu++)
   {
      if ((*dp).dp[mu][0].type==EXP)
         r0=2.0*pi/(double)(n[mu]);
      else
         r0=pi/(double)(n[mu]);
      
      for (k=0;k<n[mu];k++)
      {
         if (2*k+b[mu]!=0)
         {
            r1=2.0*sin(0.25*r0*(double)(2*k+b[mu]));
            phat[mu][k]=r1*r1;
         }
         else
            phat[mu][k]=0.0;
      }
   }

   for (k0=0;k0<nx[0];k0++)
   for (k1=0;k1<nx[1];k1++)
   for (k2=0;k2<nx[2];k2++)
   for (k3=0;k3<nx[3];k3++)
   {
      p0=phat[0][k0+rcpr[0]*L0];
      p1=phat[1][k1+rcpr[1]*L1];
      p2=phat[2][k2+rcpr[2]*L2];
      p3=phat[3][k3+rcpr[3]*L3];

      k=k3+nx[3]*k2+nx[2]*nx[3]*k1+nx[1]*nx[2]*nx[3]*k0;
      phatsq[k]=p0+p1+p2+p3;
      if (phatsq[k]==0.0)
         phatsq[k]=1.0;
   }
   
   free(phat[0]);
}


void u1mom_Delta_no0(int flag,double* inmom,double* outmom)
{
   static int init=0;
   int mu,ik,ix;
   double p,wgt;

   if (init==0)
   {
      alloc_mom();
      set_idp();
      set_phat();
      init=1;
   }

   check_global_int("u1mom_Delta_no0",1,flag);   

   wgt=sqrt(2.0);
      
   if ((bc_type()==0)||(bc_type()==2))
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==0)
         {
            inmom[8*(ix-VOLUME/2)+2]*=wgt;
            inmom[8*(ix-VOLUME/2)+3]*=wgt;
            inmom[8*(ix-VOLUME/2)+4]*=wgt;
            inmom[8*(ix-VOLUME/2)+5]*=wgt;
            inmom[8*(ix-VOLUME/2)+6]*=wgt;
            inmom[8*(ix-VOLUME/2)+7]*=wgt;
         }
      }
   }
      
   if (bc_type()==0)
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==L0*NPROC0-1)
         {
            inmom[8*(ix-VOLUME/2)+2]*=wgt;
            inmom[8*(ix-VOLUME/2)+3]*=wgt;
            inmom[8*(ix-VOLUME/2)+4]*=wgt;
            inmom[8*(ix-VOLUME/2)+5]*=wgt;
            inmom[8*(ix-VOLUME/2)+6]*=wgt;
            inmom[8*(ix-VOLUME/2)+7]*=wgt;
         }
      }
   }
   
   gather_u1mom(inmom,mom);
   for (mu=0;mu<4;mu++)
      dft4d(idp[mu],mom[mu],mom[mu]);

   if ((bc_type()==0)||(bc_type()==2))
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==0)
         {
            inmom[8*(ix-VOLUME/2)+2]/=wgt;
            inmom[8*(ix-VOLUME/2)+3]/=wgt;
            inmom[8*(ix-VOLUME/2)+4]/=wgt;
            inmom[8*(ix-VOLUME/2)+5]/=wgt;
            inmom[8*(ix-VOLUME/2)+6]/=wgt;
            inmom[8*(ix-VOLUME/2)+7]/=wgt;
         }
      }
   }
      
   if (bc_type()==0)
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==L0*NPROC0-1)
         {
            inmom[8*(ix-VOLUME/2)+2]/=wgt;
            inmom[8*(ix-VOLUME/2)+3]/=wgt;
            inmom[8*(ix-VOLUME/2)+4]/=wgt;
            inmom[8*(ix-VOLUME/2)+5]/=wgt;
            inmom[8*(ix-VOLUME/2)+6]/=wgt;
            inmom[8*(ix-VOLUME/2)+7]/=wgt;
         }
      }
   }
   
   if (flag==0)
   {
      for (ik=0;ik<vol;ik++)
      {
         p=phatsq[ik];
         for (mu=0;mu<4;mu++)
         {
            mom[mu][ik].re*=p;
            mom[mu][ik].im*=p;
         }
      }
   }
   else if (flag==1)
   {
      for (ik=0;ik<vol;ik++)
      {
         p=1.0/phatsq[ik];
         for (mu=0;mu<4;mu++)
         {
            mom[mu][ik].re*=p;
            mom[mu][ik].im*=p;
         }
      }
   }
   else if (flag==2)
   {
      for (ik=0;ik<vol;ik++)
      {
         p=sqrt(phatsq[ik]);
         for (mu=0;mu<4;mu++)
         {
            mom[mu][ik].re*=p;
            mom[mu][ik].im*=p;
         }
      }
   }
   else if (flag==3)
   {
      for (ik=0;ik<vol;ik++)
      {
         p=1.0/sqrt(phatsq[ik]);
         for (mu=0;mu<4;mu++)
         {
            mom[mu][ik].re*=p;
            mom[mu][ik].im*=p;
         }
      }
   }
   for (mu=0;mu<4;mu++)
      inv_dft4d(idp[mu],mom[mu],mom[mu]);   
   scatter_u1mom(mom,outmom);

   if (bc_type()==0)
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==0)
            outmom[8*(ix-VOLUME/2)+1]=0.0;
         if(global_time(ix)==L0*NPROC0-1)
            outmom[8*(ix-VOLUME/2)+0]=0.0;
      }
   }

   if ((bc_type()==0)||(bc_type()==2))
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==0)
         {
            outmom[8*(ix-VOLUME/2)+2]/=wgt;
            outmom[8*(ix-VOLUME/2)+3]/=wgt;
            outmom[8*(ix-VOLUME/2)+4]/=wgt;
            outmom[8*(ix-VOLUME/2)+5]/=wgt;
            outmom[8*(ix-VOLUME/2)+6]/=wgt;
            outmom[8*(ix-VOLUME/2)+7]/=wgt;
         }
      }
   }
      
   if (bc_type()==0)
   {
      for(ix=VOLUME/2;ix<VOLUME;ix++)
      {
         if(global_time(ix)==L0*NPROC0-1)
         {
            outmom[8*(ix-VOLUME/2)+2]/=wgt;
            outmom[8*(ix-VOLUME/2)+3]/=wgt;
            outmom[8*(ix-VOLUME/2)+4]/=wgt;
            outmom[8*(ix-VOLUME/2)+5]/=wgt;
            outmom[8*(ix-VOLUME/2)+6]/=wgt;
            outmom[8*(ix-VOLUME/2)+7]/=wgt;
         }
      }
   }

}
