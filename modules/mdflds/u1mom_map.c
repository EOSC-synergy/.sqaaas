
/*******************************************************************************
*
* File u1mom_map.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* Based on NSPT-1.4/modules/aflds/admap.c, NSPT-1.4/modules/aflds/adcom.c
* Copyright (C) 2015 Martin Luescher
*               
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reordering of the U(1) momentum variables in a form suitable for the FFT.
*
*   void gather_u1mom(double* inmom,complex_dble **mom)
*     Gathers the momemta inmom, stored in the "gauge friendly" format, and
*     assigns them to the complex array elements mom[mu][ix],
*     mu=0,..,3, ix=0,..,VOLUME, as described in the notes.
*
*   void scatter_u1mom(complex_dble **mom,double* outmom)
*     Inverse of the operation performed by gather_u1mom().
*
* Momenta are locally stored as the gauge potential in the "gauge friendly"
* format, see main/README.global. The program gather_u1mom() assigns the
* momentum variable inmom_{mu}(x) to the real part of the complex array element
* mom[mu][ix] and sets the imaginary part to zero, where
*
*  ix=x3+L3*x2+L2*L3*x1+L1*L2*L3*x0, 0<=xmu<Lmu.
*
* The array mom must be of size greater or equal to 4*VOLUME while the arrays
* inmom and outmom (that can be equal) must be allocated with size equal at
* least to 4*VOLUME+7*(BNDRY/4).
*
* The programs in this module act globally and must be called on all MPI
* processes simultaneously.
*
*******************************************************************************/

#define U1MOM_MAP_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "lattice.h"
#include "flags.h"
#include "u1flds.h"
#include "mdflds.h"
#include "global.h"

static int bc,np[4];
static double *buf=NULL,*sbuf,*rbuf;
static uidx_t *idx;

static void alloc_buf(void)
{
   int mu,nuk,n;

   bc=bc_type();
   np[0]=np[1]=np[2]=np[3]=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   if ((bc_cstar()>=2)&&(cpr[1]>=NPROC1/2))
      np[2]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+NPROC2+cpr[3])&0x1;
   if ((bc_cstar()>=3)&&(cpr[1]>=NPROC1/2))
      np[3]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+cpr[3]+NPROC3)&0x1;
   idx=uidx();
   n=0;

   for (mu=0;mu<4;mu++)
   {
      nuk=idx[mu].nuk;
      
      if (nuk>n)
         n=nuk;
   }
   
   buf=amalloc(2*n*sizeof(*buf),ALIGN);
   error(buf==NULL,1,"alloc_buf [u1mom_map.c]",
         "Unable to allocate send buffer");
}


static void pack_ad0(int mu,double *udb)
{
   int nu0,*iu,*ium;
   double *u;

   nu0=idx[mu].nu0;

   if (nu0>0)
   {
      u=buf;
      iu=idx[mu].iu0;
      ium=iu+nu0;

      for (;iu<ium;iu++)
      {
         (*u)=udb[*iu];
         u+=1;
      }
   }
}


static void unpack_ad0(int mu,double *udb)
{
   int nu0,*iu,*ium;
   double *u;

   nu0=idx[mu].nu0;

   if (nu0>0)
   {
      iu=idx[mu].iu0;
      ium=iu+nu0;
      u=buf;

      for (;iu<ium;iu++)
      {
         udb[*iu]=(*u);
         u+=1;
      }
   }
}


static void send_ad0(int ifc)
{
   int mu,nu0,nbf;
   int tag,saddr,raddr;
   MPI_Status stat;

   mu=ifc/2;
   nu0=idx[mu].nu0;

   if (nu0>0)
   {
      tag=mpi_tag();
      saddr=npr[ifc];
      raddr=npr[ifc^0x1];
      nbf=nu0;

      if (np[mu]==0)
      {
         MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
      }
      else
      {
         MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
      }

      rbuf+=nu0;
   }
}


static void gather_bnd_ad0(double *ad)
{
   int mu;

   if (NPROC>1)
   {
      if (buf==NULL)
         alloc_buf();

      sbuf=buf;
      rbuf=ad+4*VOLUME;

      for (mu=0;mu<4;mu++)
      {
         pack_ad0(mu,ad);
         send_ad0(2*mu);
      }
   }
}


static void scatter_bnd_ad0(double *ad)
{
   int mu,nu0;

   if (NPROC>1)
   {
      if (buf==NULL)
         alloc_buf();

      sbuf=ad+4*VOLUME;
      rbuf=buf;

      for (mu=0;mu<4;mu++)
      {
         nu0=idx[mu].nu0;
         send_ad0(2*mu+1);
         unpack_ad0(mu,ad);

         if (nu0>0)
         {
            rbuf-=nu0;
            sbuf+=nu0;
         }
      }
   }
}


void gather_u1mom(double* inmom,complex_dble **mom)
{
   int ix,iy,mu;

   gather_bnd_ad0(inmom);

   for (ix=0;ix<VOLUME;ix++)
   {
      iy=ipt[ix];

      for (mu=0;mu<4;mu++)
      {
         mom[mu][ix].re=inmom[link_offset(iy,mu)];
         mom[mu][ix].im=0.0;
      }
   }
}


void scatter_u1mom(complex_dble **mom,double* outmom)
{
   int ix,iy,mu;

   for (ix=0;ix<VOLUME;ix++)
   {
      iy=ipt[ix];

      for (mu=0;mu<4;mu++)
         outmom[link_offset(iy,mu)]=mom[mu][ix].re;
   }

   scatter_bnd_ad0(outmom);
}
