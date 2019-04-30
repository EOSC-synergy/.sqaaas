
/*******************************************************************************
*
* File adcom.c
*
* Copyright (C) 2016 Marina Marinkovic
*               2017 Agostino Patella
*
* Based on openQCD-1.6/modules/udflds/udcom.c
* Copyright (C) 2005, 2009, 2010, 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication of the non-compact U(1) link variables residing at the
* exterior boundaries of the local lattices.
*
* The externally accessible function is
*
*   void copy_bnd_ad(void)
*     Copies the non-compact U(1) link variables from the neighbouring MPI
*     processes to the exterior boundaries of the local lattice. The field
*     variables on the spatial links at time NPROC0*L0 are fetched only in
*     the case of periodic boundary conditions.
*
* The program in this module performs global communications and must be
* called simultaneously on all MPI processes.
*
*******************************************************************************/

#define ADCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "lattice.h"
#include "u1flds.h"
#include "global.h"

static int bc,np[4];
static double *sbuf=NULL,*rbuf;
static uidx_t *idx;

static void alloc_sbuf(void)
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
   
   sbuf=amalloc(n*sizeof(*sbuf),ALIGN);
   error(sbuf==NULL,1,"alloc_sbuf [adcom.c]",
         "Unable to allocate send buffer");
}


static void pack_ad0(int mu)
{
   int nu0,*iu,*ium;
   double *u,*udb;

   udb=adfld();
   nu0=idx[mu].nu0;

   if (nu0>0)
   {
      u=sbuf;
      iu=idx[mu].iu0;
      ium=iu+nu0;

      for (;iu<ium;iu++)
      {
         (*u)=udb[*iu];
         u+=1;
      }
   }
}


static void pack_adk(int mu)
{
   int nuk,*iu,*ium;
   double *u,*udb;

   udb=adfld();
   nuk=idx[mu].nuk;

   if ((nuk>0)&&((mu>0)||(cpr[0]>0)||(bc==3)))
   {
      u=sbuf;
      iu=idx[mu].iuk;
      ium=iu+nuk;

      for (;iu<ium;iu++)
      {
         (*u)=udb[*iu];
         u+=1;
      }
   }
}


static void send_ad0(int mu)
{
   int nu0,nbf;
   int tag,saddr,raddr;
   MPI_Status stat;

   nu0=idx[mu].nu0;

   if (nu0>0)
   {
      tag=mpi_tag();
      saddr=npr[2*mu];
      raddr=npr[2*mu+1];
      nbf=nu0;   /*18 for SU(3) replaced with 1 for non-compact U(1)*/

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


static void send_adk(int mu)
{
   int nuk,nbf;
   int tag,saddr,raddr;
   MPI_Status stat;

   nuk=idx[mu].nuk;

   if (nuk>0)
   {
      tag=mpi_tag();
      saddr=npr[2*mu];
      raddr=npr[2*mu+1];
      nbf=nuk; /*18 for SU(3) replaced with 1 for non-compact U(1)*/

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

      rbuf+=nuk;
   }
}


void copy_bnd_ad(void)
{
   int mu;
   
   if (NPROC>1)
   {
      if (sbuf==NULL)
         alloc_sbuf();

      rbuf=adfld()+4*VOLUME;

      for (mu=0;mu<4;mu++)
      {
         pack_ad0(mu);
         send_ad0(mu);
      }

      for (mu=0;mu<4;mu++)
      {      
         pack_adk(mu);
         send_adk(mu);
      }
   }

   set_flags(COPIED_BND_AD);
}
