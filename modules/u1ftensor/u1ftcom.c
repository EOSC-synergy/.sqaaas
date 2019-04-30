
/*******************************************************************************
*
* File u1ftcom.c
*
* Copyright (C) 2016 Nazario Tantalo
*
* Based on openQCD-1.4/modules/tcharge/ftcom.c
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication of the U(1) field tensor components residing at the 
* boundaries of the local lattices.
*
* The externally accessible functions are
*
*   void copy_bnd_u1ft(int n,double *ft)
*     Fetches the boundary values the field ft from the neighbouring MPI
*     processes (see the notes). The boundary values at time NPROC0*L0
*     are fetched from the field at time 0 only in the case of periodic
*     boundary conditions.
*
*   void add_bnd_u1ft(int n,double *ft)
*     Adds the boundary values of the field ft to the field on the
*     neighbouring MPI processes. The boundary values at time NPROC0*L0
*     are added to the field at time 0 only in the case of periodic
*     boundary conditions.
*
* Notes:
*
* Both communication programs assume that the field ft has the same size as
* the n-th component of the symmetric field tensor F_{mu nu}, where n=0,..,5
* labels the (mu,nu)-planes (0,1),(0,2),(0,3),(2,3),(3,1),(1,2). For further
* explanations, see the files lattice/README.ftidx and u1tcharge/u1ftensor.c.
*
* The programs in this module perform global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define U1FTCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1ftensor.h"
#include "global.h"

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static double *u1ftbuf;
static ftidx_t *idx=NULL;


static void alloc_u1ftbuf(void)
{
   int n,nft,nbf;
   
   idx=ftidx();
   nbf=0;

   for (n=0;n<6;n++)
   {
      nft=idx[n].nft[0];
      if (nft>nbf)
         nbf=nft;
      
      nft=idx[n].nft[1];
      if (nft>nbf)
         nbf=nft;
   }

   u1ftbuf=amalloc(nbf*sizeof(*u1ftbuf),ALIGN);
   error(u1ftbuf==NULL,1,"alloc_u1ftbuf [u1ftcom.c]",
         "Unable to allocate communication buffers");
}


static void pack_buf(int n,int dir,double *ft)
{
   int bc,mu,nft;
   int *ift,*ifm;
   double *fb;

   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];

      if ((mu>0)||(cpr[0]>0)||(bc==3))
      {
         ift=idx[n].ift[dir];
         ifm=ift+nft;
         fb=u1ftbuf;

         for (;ift<ifm;ift++)
         {
            fb[0]=ft[*ift];
            fb+=1;
         }
      }
   }
}


static void fwd_send(int n,int dir,double *ft)
{
   int bc,mu,nft,nbf;
   int tag,saddr,raddr,np[4];
   double *sbuf,*rbuf;
   MPI_Status stat;

   np[0]=np[1]=np[2]=np[3]=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   if ((bc_cstar()>=2)&&(cpr[1]>=NPROC1/2))
      np[2]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+NPROC2+cpr[3])&0x1;
   if ((bc_cstar()>=3)&&(cpr[1]>=NPROC1/2))
      np[3]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+cpr[3]+NPROC3)&0x1;
   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];
      tag=mpi_tag();
      saddr=npr[2*mu];
      raddr=npr[2*mu+1];
      sbuf=u1ftbuf;
      rbuf=ft+VOLUME;
      if (dir==1)
         rbuf+=idx[n].nft[0];
      nbf=nft;

      if (np[mu]==0)
      {
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
      }
      else
      {
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))         
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
      }
   }
}


void copy_bnd_u1ft(int n,double *ft)
{
   if (NPROC>1)
   {
      if (idx==NULL)
	 alloc_u1ftbuf();
      
      pack_buf(n,1,ft);
      fwd_send(n,1,ft);
      pack_buf(n,0,ft);
      fwd_send(n,0,ft);
   }
}


static void bck_send(int n,int dir,double *ft)
{
   int bc,mu,nft,nbf;
   int tag,saddr,raddr,np[4];
   double *sbuf,*rbuf;
   MPI_Status stat;

   np[0]=np[1]=np[2]=np[3]=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   if ((bc_cstar()>=2)&&(cpr[1]>=NPROC1/2))
      np[2]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+NPROC2+cpr[3])&0x1;
   if ((bc_cstar()>=3)&&(cpr[1]>=NPROC1/2))
      np[3]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+cpr[3]+NPROC3)&0x1;
   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];
      tag=mpi_tag();
      saddr=npr[2*mu+1];
      raddr=npr[2*mu];
      sbuf=ft+VOLUME;
      if (dir==1)
         sbuf+=idx[n].nft[0];
      rbuf=u1ftbuf;
      nbf=nft;

      if (np[mu]==0)
      {
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         if ((mu>0)||(cpr[0]>0)||(bc==3))
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
      }
      else
      {
         if ((mu>0)||(cpr[0]>0)||(bc==3))         
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         if ((mu>0)||(cpr[0]<(NPROC0-1))||(bc==3))
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
      }
   }
}


static void unpack_buf(int n,int dir,double *ft)
{
   int bc,mu,nft;
   int *ift,*ifm;
   double *f,*fb;

   nft=idx[n].nft[dir];

   if (nft>0)
   {
      bc=bc_type();
      mu=plns[n][dir];   

      if ((mu>0)||(cpr[0]>0)||(bc==3))
      {
         ift=idx[n].ift[dir];
         ifm=ift+nft;
         fb=u1ftbuf;

         for (;ift<ifm;ift++)
         {
            f=ft+(*ift);

            (*f)+=(*fb);
      
            fb+=1;
         }
      }
   }
}


void add_bnd_u1ft(int n,double *ft)
{
   if (NPROC>1)
   {
      if (idx==NULL)
         alloc_u1ftbuf();         
   
      bck_send(n,0,ft);
      unpack_buf(n,0,ft);
      bck_send(n,1,ft);
      unpack_buf(n,1,ft);
   }
}
