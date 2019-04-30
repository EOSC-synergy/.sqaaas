/*******************************************************************************
*
* File force6.c
*
* Copyright (C) 2016, 2017 Alberto Ramos, Agostino Patella
*
* Based on openQCD-1.6/modules/forces/force0.c
* Copyright (C) 2005, 2009-2014, 2016 Martin Luescher, John Bulava
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Compact action of the double-precision U(1) gauge field and 
* associated force.
*
* The externally accessible functions are
*
*   void plaq_u1frc(void)
*     Computes the force deriving from the compact U(1) plaquette action,
*     omitting the prefactor 1/(qel*e0)^2, and assigns the result to the MD
*     force field. In the case of open, SF or open-SF boundary conditions,
*     the boundary improvement coefficients are set to their tree-level
*     value independently of the values stored in the parameter data base.
*
*   void force6(double c)
*     Computes the force deriving from the compact U(1) gauge action, 
*     including the prefactor beta_e/2, multiplies the calculated force 
*     by c and assigns the result to the MD force field. The coupling 
*     g0 and the other parameters of the gauge action are retrieved 
*     from the parameter data base.
*
*   double action6(int icom)
*     Computes the local part of the U(1) compact gauge action including 
*     the prefactor beta_e/2. The coupling g0 and the other parameters of 
*     the action are retrieved from the parameter data base. The 
*     program returns the sum of the local parts of the action over 
*     all MPI processes if icom=1 and otherwise just the local part.
*
*******************************************************************************/


#define FORCE6_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "flags.h"
#include "u1flds.h"
#include "utils.h"
#include "lattice.h"
#include "mdflds.h"
#include "forces.h"
#include "global.h"

#define N0 (NPROC0*L0)

static const int plns[6][2]={{0,1},{0,2},{0,3},{2,3},{3,1},{1,2}};
static int nfc[8],ofs[8],hofs[8],init=0,ism;
static complex_dble *u1db,*h1db;
static complex_dble wd[4] ALIGNED16;
static complex_dble vd[4] ALIGNED16;
static double X ALIGNED16;
static double Y[4] ALIGNED16;

static void set_ofs(void)
{
   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=VOLUME;
   ofs[1]=ofs[0]+(FACE0/2);
   ofs[2]=ofs[1]+(FACE0/2);
   ofs[3]=ofs[2]+(FACE1/2);
   ofs[4]=ofs[3]+(FACE1/2);
   ofs[5]=ofs[4]+(FACE2/2);
   ofs[6]=ofs[5]+(FACE2/2);
   ofs[7]=ofs[6]+(FACE3/2);

   hofs[0]=0;
   hofs[1]=hofs[0]+3*FACE0;
   hofs[2]=hofs[1]+3*FACE0;
   hofs[3]=hofs[2]+3*FACE1;
   hofs[4]=hofs[3]+3*FACE1;
   hofs[5]=hofs[4]+3*FACE2;
   hofs[6]=hofs[5]+3*FACE2;
   hofs[7]=hofs[6]+3*FACE3;

   init+=2;
}

static void set_staples(int n,int ix,int ia)
{
   int mu,nu,ifc;
   int iy,ib,ip[4];

   mu=plns[n][0];
   nu=plns[n][1];

   if (!ia)
   {
      iy=idn[ix][nu];

      if (iy<VOLUME)
      {
         plaq_uidx(n,iy,ip);

         u1xu1(u1db+ip[0],u1db+ip[1],wd+2);
         u1dagxu1(u1db+ip[2],wd+2,vd);
      }
      else
      {
         ifc=2*nu;

         if (iy<(VOLUME+(BNDRY/2)))
            ib=iy-ofs[ifc];
         else
            ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

         vd[0]=h1db[hofs[ifc]+3*ib+mu-(mu>nu)];
      }
   }

   iy=iup[ix][mu];

   if (iy<VOLUME)
   {
      plaq_uidx(n,iy,ip);

      u1xu1dag(u1db+ip[1],u1db+ip[3],wd+2);
      u1xu1(u1db+ip[0],wd+2,vd+1);
   }
   else
   {
      ifc=2*mu+1;

      if (iy<(VOLUME+(BNDRY/2)))
         ib=iy-ofs[ifc];
      else
         ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

      vd[1]=h1db[hofs[ifc]+3*ib+nu-(nu>mu)];
   }

   if (!ia)
   {
      iy=idn[ix][mu];

      if (iy<VOLUME)
      {
         plaq_uidx(n,iy,ip);

         u1xu1(u1db+ip[2],u1db+ip[3],wd+2);
         u1dagxu1(u1db+ip[0],wd+2,vd+2);
      }
      else
      {
         ifc=2*mu;

         if (iy<(VOLUME+(BNDRY/2)))
            ib=iy-ofs[ifc];
         else
            ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

         vd[2]=h1db[hofs[ifc]+3*ib+nu-(nu>mu)];
      }
   }

   iy=iup[ix][nu];

   if (iy<VOLUME)
   {
      plaq_uidx(n,iy,ip);

      u1xu1dag(u1db+ip[3],u1db+ip[1],wd+2);
      u1xu1(u1db+ip[2],wd+2,vd+3);
   }
   else
   {
      ifc=2*nu+1;

      if (iy<(VOLUME+(BNDRY/2)))
         ib=iy-ofs[ifc];
      else
         ib=iy-ofs[ifc]-(BNDRY/2)+nfc[ifc];

      vd[3]=h1db[hofs[ifc]+3*ib+mu-(mu>nu)];
   }
}


void plaq_u1frc(void)
{
   int bc,n,ix,t,ip[4];
   double r;
   double *fdb;
   mdflds_t *mdfs;

   error_root((gauge()&2)==0,1,
              "plaq_u1frc [force6.c]",
              "U(1) gauge field is not active");

   if (query_flags(ADBUF_UP2DATE)!=1)
      copy_bnd_ad();

   error_root(query_flags(ADBUF_UP2DATE)!=1,1,
              "plaq_u1frc [force6.c]",
              "Flags debug: Something is not working");

   bc=bc_type();
   u1db=u1dfld(EXT);
   mdfs=mdflds();
   fdb=(*mdfs).u1frc;
   set_u1frc2zero();

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t<(N0-1))||(bc!=0))
      {
         for (n=0;n<3;n++)
         {
            plaq_uidx(n,ix,ip);
            u1xu1dag(u1db+ip[1],u1db+ip[3],wd);
            u1dagxu1(u1db+ip[2],u1db+ip[0],wd+1);
            cm1x1_imtr(wd,wd+1,&X);

            if ((t<(N0-1))||(bc==3))
               *(fdb+ip[1])+=X;

            *(fdb+ip[3])-=X;
            *(fdb+ip[0])+=X;

            if ((t>0)||(bc!=1))
               *(fdb+ip[2])-=X;
         }
      }

      if ((t>0)||(bc!=1))
      {
         r=1.0;

         if (((t==0)&&(bc!=3))||((t==(N0-1))&&(bc==0)))
            r=0.5;

         for (n=3;n<6;n++)
         {
            plaq_uidx(n,ix,ip);
            u1xu1dag(u1db+ip[1],u1db+ip[3],wd);
            u1dagxu1(u1db+ip[2],u1db+ip[0],wd+1);
            cm1x1_imtr(wd,wd+1,&X);
            *(fdb+ip[1])+=r*X;
            *(fdb+ip[3])-=r*X;
            *(fdb+ip[0])+=r*X;
            *(fdb+ip[2])-=r*X;
         }
      }
   }

   add_bnd_u1frc();
}


void force6(double c)
{
   int bc,n,ix,t,ip[4],is_typeB;
   double c0,c1,*cG;
   double r0,r1,bf;
   double *fdb;
   mdflds_t *mdfs;
   u1lat_parms_t lat;
   bc_parms_t bcp;

   error_root((gauge()&2)==0,1,
              "force6 [force6.c]",
              "U(1) gauge field is not active");

   bcp=bc_parms();
   bc=bcp.type;

   lat=u1lat_parms();
   c*=lat.beta;
   if (bcp.cstar!=0) c*=0.5;
   c0=lat.c0;
   c1=lat.c1;
   cG=lat.cG;
   is_typeB=lat.SFtype;

   error_root(lat.type==1,1,
              "force6 [force6.c]",
              "Non-compact U(1) gauge action requested. I should not be here!");

   if (query_flags(ADBUF_UP2DATE)!=1)
      copy_bnd_ad();

   error_root(query_flags(ADBUF_UP2DATE)!=1,1,
              "force6 [force6.c]",
              "Flags debug: Something is not working");

   u1db=u1dfld(EXT);
   mdfs=mdflds();
   fdb=(*mdfs).u1frc;
   set_u1frc2zero();

   if (c0==1.0)
      h1db=NULL;
   else
   {
      if ((init&0x2)==0)
         set_ofs();

      if (query_flags(U1BSTAP_UP2DATE)!=1)
         set_u1_bstap();
      h1db=u1_bstap();
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t<(N0-1))||(bc!=0))
      {
         r0=c*c0;
         r1=c*c1;

         if ((t==0)&&(bc==1))
            r0*=cG[0];
         else if ((t==(N0-1))&&(bc!=3))
            r0*=cG[1];

         for (n=0;n<3;n++)
         {
            plaq_uidx(n,ix,ip);

            u1xu1dag(u1db+ip[1],u1db+ip[3],wd);
            u1xu1dag(wd,u1db+ip[2],wd+1);
            cm1x1_imtr(wd+1,u1db+ip[0],&X);
            *(fdb+ip[3]) -= r0*X;
            *(fdb+ip[0]) += r0*X;

            if ((t<(N0-1))||(bc==3))
            {
               *(fdb+ip[1]) += r0*X;
            }
            if ((t>0)||(bc!=1))
            {
               *(fdb+ip[2]) -= r0*X;
            }

            if (c0!=1.0)
            {
               if ((t==0)&&(bc==1)&&(!is_typeB))
               {
                  u1xu1(wd+1,u1db+ip[0],wd+2);
                  cm1x1_imtr(wd+2,wd+2,&X);
                  *(fdb+ip[1]) += r1*X;
                  *(fdb+ip[0]) += r1*X;
                  *(fdb+ip[3]) -= r1*X;
               }

               if ((t==(N0-1))&&(bc!=3)&&(!is_typeB))
               {
                  u1xu1(wd+1,u1db+ip[0],wd+2);
                  cm1x1_imtr(wd+2,wd+2,&X);
                  *(fdb+ip[0]) += r1*X;
                  *(fdb+ip[2]) -= r1*X;
                  *(fdb+ip[3]) -= r1*X;
               }

               set_staples(n,ix,0);

               if ((is_typeB)&&(t==0)&&(bc==1))
                  bf=1.5;
               else
                  bf=1.0;

               if ((is_typeB)&&(t==(N0-1))&&((bc==1)||(bc==2)))
                  bf*=1.5;
               else
                  bf*=1.0;

               cm1x1_imtr(wd+1,vd,&Y[0]);

               u1xu1dag(wd,vd+2,wd+2);
               cm1x1_imtr(wd+2,u1db+ip[0],&Y[2]);

               u1dagxu1(u1db+ip[2],u1db+ip[0],wd+2);
               u1dagxu1(vd+3,wd+2,wd+3);
               cm1x1_imtr(wd+3,u1db+ip[1],&Y[3]);

               u1dagxu1(u1db+ip[3],wd+2,wd+3);
               cm1x1_imtr(wd+3,vd+1,&Y[1]); 

               if ((t<(N0-1))||(bc==3))
               {
                  *(fdb+ip[1]) += bf*r1*Y[0];
               }

               if ((t>0)||(bc!=1))
               {
                  cm1x1_retr(vd,wd+1,&X);
                  *(fdb+ip[2]) -= bf*r1*Y[0];
               }

               *(fdb+ip[3]) -= bf*r1*Y[0];

               if ((t<(N0-2))||((t==(N0-2))&&(bc!=0))||(bc==3))
               {
                  *(fdb+ip[0]) += r1*Y[1];

                  if ((t>0)||(bc!=1))
                  {
                     *(fdb+ip[2]) -= r1*Y[1];
                  }

                  *(fdb+ip[3]) -= r1*Y[1];
               }

               if ((t>0)||(bc==3))
               {
                  *(fdb+ip[0]) += r1*Y[2];

                  if ((t<(N0-1))||(bc==3))
                  {
                     *(fdb+ip[1]) += r1*Y[2];
                  }

                  *(fdb+ip[3]) -= r1*Y[2];
               }

               *(fdb+ip[0]) += bf*r1*Y[3];

               if ((t>0)||(bc!=1))
               {
                  *(fdb+ip[2]) -= bf*r1*Y[3];
               }

               if ((t<(N0-1))||(bc==3))
               {
                  *(fdb+ip[1]) += bf*r1*Y[3];
               }
            }
         }
      }

      if ((t>0)||(bc!=1))
      {
         r0=c*c0;
         r1=c*c1;

         if ((t==0)&&(bc!=3))
         {
            r0*=(0.5*cG[0]);
            r1*=(0.5*cG[0]);
         }
         else if ((t==(N0-1))&&(bc==0))
         {
            r0*=(0.5*cG[1]);
            r1*=(0.5*cG[1]);
         }

         for (n=3;n<6;n++)
         {
            plaq_uidx(n,ix,ip);

            u1xu1dag(u1db+ip[1],u1db+ip[3],wd);
            u1xu1dag(wd,u1db+ip[2],wd+1);
            cm1x1_imtr(wd+1,u1db+ip[0],&X);
            *(fdb+ip[0]) += r0*X;
            *(fdb+ip[1]) += r0*X;
            *(fdb+ip[2]) -= r0*X;
            *(fdb+ip[3]) -= r0*X;

            if (c0!=1.0)
            {
               set_staples(n,ix,0);

               cm1x1_imtr(wd+1,vd,&Y[0]);

               u1xu1dag(wd,vd+2,wd+2);
               cm1x1_imtr(wd+2,u1db+ip[0],&Y[2]);

               u1dagxu1(u1db+ip[2],u1db+ip[0],wd+2);
               u1dagxu1(vd+3,wd+2,wd+3);
               cm1x1_imtr(wd+3,u1db+ip[1],&Y[3]);

               u1dagxu1(u1db+ip[3],wd+2,wd+3);
               cm1x1_imtr(wd+3,vd+1,&Y[1]); 

               *(fdb+ip[0]) += r1*(Y[1]+Y[2]+Y[3]);
               *(fdb+ip[1]) += r1*(Y[0]+Y[2]+Y[3]);
               *(fdb+ip[2]) -= r1*(Y[0]+Y[1]+Y[3]);
               *(fdb+ip[3]) -= r1*(Y[0]+Y[2]+Y[1]);
            }

         }
      }

   }

   add_bnd_u1frc();
}


static void wloops(int n,int ix,int t,double c0,double *trU)
{
   int bc,ip[4];

   bc=bc_type();
   plaq_uidx(n,ix,ip);

   trU[0]=0.0;
   trU[1]=0.0;
   trU[2]=0.0;
   trU[3]=0.0;

   if ((n>=3)||(t<(N0-1))||(bc!=0))
   {
      u1dagxu1(u1db+ip[2],u1db+ip[0],wd);
      u1xu1dag(u1db+ip[1],u1db+ip[3],wd+1);
      cm1x1_retr(wd,wd+1,trU);
      trU[0]=1.0-trU[0];
   }

   if (c0!=1.0)
   {
      set_staples(n,ix,1);

      if ((n<3)&&(((t==0)&&(bc==1))||
      ((t==(N0-1))&&((bc==1)||(bc==2)))))
      {
         u1xu1(wd,wd+1,wd+1);
         cm1x1_retr(wd+1,wd+1,trU+3);
         trU[3]=1.0-trU[3];
      }

      if ((n>=3)||(t<(N0-1))||(bc!=0))
      {
         u1xu1dag(u1db+ip[1],vd+3,wd+1);
         cm1x1_retr(wd,wd+1,trU+1);
         trU[1]=1.0-trU[1];
      }

      if ((n>=3)||(t<(N0-2))||((t==(N0-2))&&(bc!=0))||(bc==3))
      {
         u1xu1dag(vd+1,u1db+ip[3],wd+1);
         cm1x1_retr(wd,wd+1,trU+2);
         trU[2]=1.0-trU[2];
      }
   }
}


double action6(int icom)
{
   int bc,ix,t,n,is_typeB;
   double beta,c0,c1,*cG;
   double r0,r1,trU[4],act,bf1,bf3;
   u1lat_parms_t lat;
   bc_parms_t bcp;
   
   error_root((gauge()&2)==0,1,
              "action6 [force6.c]",
              "U(1) gauge field is not active");

   bcp=bc_parms();
   bc=bcp.type;

   lat=u1lat_parms();
   beta=lat.beta;
   if (bcp.cstar!=0) beta*=0.5;
   c0=lat.c0;
   c1=lat.c1;
   cG=lat.cG;
   is_typeB=lat.SFtype;

   error_root(lat.type==1,1,
              "action6 [force6.c]",
              "Non-compact U(1) gauge action requested. I should not be here!");

   u1db=u1dfld(EXT);

   if (c0==1.0)
      h1db=NULL;
   else
   {
      if ((init&0x2)==0)
         set_ofs();

      if (query_flags(U1BSTAP_UP2DATE)!=1)
         set_u1_bstap();
      h1db=u1_bstap();
   }

   if ((init&0x1)==0)
   {
      ism=init_hsum(1);
      init+=1;
   }

   reset_hsum(ism);

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      act=0.0;

      if ( ((is_typeB)&&(t==0)&&(bc==1))||
           ((is_typeB)&&(t==(N0-1))&&((bc==1)||(bc==2))) )
      {
         bf1=1.5;
         bf3=0.0;
      }
      else
      {
         bf1=1.0;
         bf3=0.5;
      }

      if ((t<(N0-1))||(bc!=0))
      {
         r0=c0;

         if ((t==0)&&(bc==1))
         r0*=cG[0];
         else if ((t==(N0-1))&&(bc!=3))
         r0*=cG[1];

         for (n=0;n<3;n++)
         {
            wloops(n,ix,t,c0,trU);
            act+=(r0*trU[0]+c1*(bf1*trU[1]+trU[2]+bf3*trU[3]));
         }
      }

      if ((t>0)||(bc!=1))
      {
         r0=c0;
         r1=c1;

         if ((t==0)&&(bc!=3))
         {
            r0*=(0.5*cG[0]);
            r1*=(0.5*cG[0]);
         }
         else if ((t==(N0-1))&&(bc==0))
         {
            r0*=(0.5*cG[1]);
            r1*=(0.5*cG[1]);
         }

         for (n=3;n<6;n++)
         {
            wloops(n,ix,t,c0,trU);
            act+=(r0*trU[0]+r1*(trU[1]+trU[2]));
         }
      }

      add_to_hsum(ism,&act);
   }

   if ((icom==1)&&(NPROC>1))
      global_hsum(ism,&act);
   else
      local_hsum(ism,&act);

   return beta*act;
}
