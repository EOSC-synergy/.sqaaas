
/*******************************************************************************
*
* File gflds_utils.c
*
* Copyright (C) 2016 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Initialization, translation and guage transformations of the gauge fields.
*
*******************************************************************************/

#define GFLDS_UTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "utils.h"
#include "lattice.h"
#include "flags.h"
#include "uflds.h"
#include "global.h"
#include "u1flds.h"
#include "gflds_utils.h"

#define FILENAME "gflds_utils.c"

#define N0 (NPROC0*L0)


static void unit_ud(void)
{
   size_t n;
   su3_dble unity,*u,*um,*udb;
   bc_parms_t bc;
   int ie;
   
   bc=bc_parms();
   
   ie=0;
   ie|=(bc.phi[1][0]!=0.0);
   ie|=(bc.phi[1][1]!=0.0);
   ie|=(bc.phi[1][2]!=0.0);
   error((bc.type==2)&ie,1,"unit_ud ["FILENAME"]","Nonzero boundary angles for open-SF boundary conditions");
   ie|=(bc.phi[0][0]!=0.0);
   ie|=(bc.phi[0][1]!=0.0);
   ie|=(bc.phi[0][2]!=0.0);
   error((bc.type==1)&ie,1,"unit_ud ["FILENAME"]","Nonzero boundary angles for SF boundary conditions");

   n=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc.type==1)||(bc.type==2)))
      n+=3;
   
   udb=udfld();

   cm3x3_unity(1,&unity);
   u=udb;
   um=udb+n;

   for (;u<um;u++)
      (*u)=unity;

   set_flags(UPDATED_UD);
   set_bc();
}


static void zero_ad(void)
{
   size_t n;
   double *a,*am,*adb;
   bc_parms_t bc;
   
   bc=bc_parms();
   
   n=4*VOLUME+7*(BNDRY/4);
   if ((cpr[0]==(NPROC0-1))&&((bc.type==1)||(bc.type==2)))
      n+=3;
   
   adb=adfld();

   a=adb;
   am=adb+n;

   for (;a<am;a++)
      (*a)=0.0;

   set_flags(UPDATED_AD);
   set_ad_bc();
}


void unit_gflds(void)
{
   unit_ud();
   zero_ad();
}


void random_gflds(void)
{
   if(gauge()==1)
   {
      random_ud();
      zero_ad();
   }
   else if(gauge()==3)
   {
      random_ud();
      random_ad();
   }
   else if(gauge()==2)
   {
      unit_ud();
      random_ad();
   }
}


int shift_gflds(int *s)
{
   int ret=0;
   
   if((gauge()&1)!=0)
      ret+=shift_ud(s);
   if((gauge()&2)!=0)
      ret+=shift_ad(s);
   
   return ret;
}


su3_dble *g3tr(void)
{
   su3_dble unity,*g,*gm;
   static su3_dble *gb=NULL;

   if(gb==NULL)
   {
      error(iup[0][0]==0,1,"g3tr ["FILENAME"]","Geometry arrays are not set");
      gb=amalloc(NSPIN*sizeof(*gb),ALIGN);
      error(gb==NULL,1,"g3tr ["FILENAME"]",
            "Unable to allocate memory space for the gauge transformation");
      
      cm3x3_unity(1,&unity);

      g=gb;
      gm=gb+NSPIN;
      for(;g<gm;g++)
         (*g)=unity;
   }

   return gb;
}


double *g1tr(void)
{
   double *g,*gm;
   static double *gb=NULL;

   if(gb==NULL)
   {
      error(iup[0][0]==0,1,"g1tr ["FILENAME"]","Geometry arrays are not set");
      gb=amalloc(NSPIN*sizeof(*gb),ALIGN);
      error(gb==NULL,1,"g1tr ["FILENAME"]",
            "Unable to allocate memory space for the gauge transformation");
   
      g=gb;
      gm=gb+NSPIN;
      for(;g<gm;g++)
         (*g)=0.;
   }

   return gb;
}


static void pack_and_send_g3buf(void)
{
   static su3_dble *gbuf=NULL;
   static int nfc[8],ofs[8];
   int ifc,ib,ix;
   int np,saddr,raddr;
   int nbf,tag;
   su3_dble *sbuf,*rbuf,*g;
   MPI_Status stat;
   
   g=g3tr();

   if (BNDRY!=0 && gbuf==NULL)
   {
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);
      error(gbuf==NULL,1,"pack_and_send_gbuf ["FILENAME"]",
            "Unable to allocate memory space for gbuf");
   }

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (ifc=0;ifc<8;ifc++)
   {
      for (ib=0;ib<nfc[ifc];ib++)
      {
         ix=map[ofs[ifc]+ib];
         gbuf[ofs[ifc]+ib]=g[ix];
      }
   }

   np=cpr[0]+cpr[1]+cpr[2]+cpr[3];

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=18*nfc[ifc];

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=npr[ifc^0x1];
         raddr=npr[ifc];
         sbuf=gbuf+ofs[ifc];
         rbuf=g+VOLUME+ofs[ifc];

         if (np&0x1)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


static void pack_and_send_g1buf(void)
{
   static double *gbuf=NULL;
   static int nfc[8],ofs[8];
   int ifc,ib,ix;
   int np,saddr,raddr;
   int nbf,tag;
   double *sbuf,*rbuf,*g;
   MPI_Status stat;
   
   g=g1tr();

   if (BNDRY!=0 && gbuf==NULL)
   {
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);
      error(gbuf==NULL,1,"pack_and_send_gbuf ["FILENAME"]",
            "Unable to allocate memory space for gbuf");
   }

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (ifc=0;ifc<8;ifc++)
   {
      for (ib=0;ib<nfc[ifc];ib++)
      {
         ix=map[ofs[ifc]+ib];
         gbuf[ofs[ifc]+ib]=g[ix];
      }
   }

   np=cpr[0]+cpr[1]+cpr[2]+cpr[3];

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=nfc[ifc];

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=npr[ifc^0x1];
         raddr=npr[ifc];
         sbuf=gbuf+ofs[ifc];
         rbuf=g+VOLUME+ofs[ifc];

         if (np&0x1)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


void random_g(void)
{
   int ix,t,bc;
   double *g1x;
   su3_dble unity,*g3x;
   static const su3_dble ud0={{0.0}};
   double twopi;

   bc=bc_type();
   
   if((gauge()&1)!=0)
   {
      unity=ud0;
      unity.c11.re=1.0;
      unity.c22.re=1.0;
      unity.c33.re=1.0;
      g3x=g3tr();

      for (ix=0;ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if ((t>0)||(bc!=1))
            random_su3_dble(g3x);
         else
            (*g3x)=unity;

         g3x+=1;
      }

      if (BNDRY>0)
         pack_and_send_g3buf();
   }
   
   if((gauge()&2)!=0)
   {
      g1x=g1tr();
      twopi=4*atan(1.);

      for (ix=0;ix<VOLUME;ix++)
      {
         t=global_time(ix);

         if ((t>0)||(bc!=1))
         {
            ranlxd(g1x,1);
            (*g1x)=((*g1x)-.5)*twopi;
         }
         else
            (*g1x)=0.;

         g1x+=1;
      }

      if (BNDRY>0)
         pack_and_send_g1buf();
   }

}


static void transform_ud(void)
{
   static su3_dble wd ALIGNED16;
   int ix,iy,t,ifc,bc;
   su3_dble *u,*g;
   
   bc=bc_type();

   u=udfld();
   g=g3tr();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         iy=iup[ix][0];
         su3xsu3dag(u,g+iy,&wd);
         su3xsu3(g+ix,&wd,u);
         u+=1;

         if (bc==3)
         {
            iy=idn[ix][0];
            su3xsu3dag(u,g+ix,&wd);
            su3xsu3(g+iy,&wd,u);
         }
         else if (bc!=0)
         {
            iy=idn[ix][0];
            su3xsu3(g+iy,u,&wd);
            (*u)=wd;
         }

         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               if (ifc&0x1)
               {
                  iy=idn[ix][ifc/2];
                  su3xsu3dag(u,g+ix,&wd);
                  su3xsu3(g+iy,&wd,u);
               }
               else
               {
                  iy=iup[ix][ifc/2];
                  su3xsu3dag(u,g+iy,&wd);
                  su3xsu3(g+ix,&wd,u);
               }
            }

            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc==3)
         {
            iy=iup[ix][0];
            su3xsu3dag(u,g+iy,&wd);
            su3xsu3(g+ix,&wd,u);
         }
         else if (bc!=0)
         {
            su3xsu3(g+ix,u,&wd);
            (*u)=wd;
         }

         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               su3xsu3dag(u,g+ix,&wd);
               su3xsu3(g+iy,&wd,u);
            }
            else
            {
               iy=iup[ix][ifc/2];
               su3xsu3dag(u,g+iy,&wd);
               su3xsu3(g+ix,&wd,u);
            }

            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               su3xsu3dag(u,g+ix,&wd);
               su3xsu3(g+iy,&wd,u);
            }
            else
            {
               iy=iup[ix][ifc/2];
               su3xsu3dag(u,g+iy,&wd);
               su3xsu3(g+ix,&wd,u);
            }

            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}



static void transform_ad(void)
{
   int ix,iy,t,ifc,bc;
   double *a,*g;
   
   bc=bc_type();

   a=adfld();
   g=g1tr();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         iy=iup[ix][0];
         (*a)=(*a)-g[iy]+g[ix];
         a+=1;

         if (bc==3)
         {
            iy=idn[ix][0];
            (*a)=(*a)-g[ix]+g[iy];
         }
         else if (bc!=0)
         {
            iy=idn[ix][0];
            (*a)=(*a)+g[iy];
         }

         a+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               if (ifc&0x1)
               {
                  iy=idn[ix][ifc/2];
                  (*a)=(*a)-g[ix]+g[iy];
               }
               else
               {
                  iy=iup[ix][ifc/2];
                  (*a)=(*a)-g[iy]+g[ix];
               }
            }

            a+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc==3)
         {
            iy=iup[ix][0];
            (*a)=(*a)-g[iy]+g[ix];
         }
         else if (bc!=0)
         {
            (*a)=(*a)+g[ix];
         }

         a+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               (*a)=(*a)-g[ix]+g[iy];
            }
            else
            {
               iy=iup[ix][ifc/2];
               (*a)=(*a)-g[iy]+g[ix];
            }

            a+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               (*a)=(*a)-g[ix]+g[iy];
            }
            else
            {
               iy=iup[ix][ifc/2];
               (*a)=(*a)-g[iy]+g[ix];
            }

            a+=1;
         }
      }
   }

   set_flags(UPDATED_AD);
}


void transform_gflds(void)
{
   if((gauge()&1)!=0)
      transform_ud();
   if((gauge()&2)!=0)
      transform_ad();
}
