
/*******************************************************************************
*
* File cfactor.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of the gauge factors for electrically charged states.
*
* The externally accessible functions are
*
*   double* get_cfactor_phase(int coulomb,int mu)
*     Returns the pointer phi to the base address of a global real field
*     containing the phase Phi_mu(x), see the notes below. The real field is
*     organized in memory as a matter field so that phi[ix] is equal to Phi_mu(x)
*     with ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0]. The argument mu>=1 must be a C*
*     spatial direction. If coulomb=1 the gauge field used in the calculation of
*     the phase contains the Coulomb term A_mu^C(x).
*
*  void mul_cfactor(int iflag,int coulomb,int mu,spinor_dble *pk,spinor_dble *pl)
*     Multiplies the matter field pk[ix] for the gauge factor,
*        pl(x) = exp(i sign qhat Phi_mu(x)/2 ) pk(x),
*     and assigns the result to pl[ix]. For the definition of Phi_mu(x) see the
*     notes below. If iflag=0, sign=1.0 while sign=-1.0 if iflag=1. The argument
*     mu>=1 must be a C* spatial direction. If coulomb=1 the gauge field used in
*     the calculation of the phase contains the Coulomb term A_mu^C(x). qhat is
*     the electric charge in units of qel of the active flavour,
*     qhat=dirac_parms().qhat.
*
*  void mul_cfactor_muaverage(int iflag,int coulomb,spinor_dble *pk,spinor_dble *pl)
*     Multiplies the field pk[ix] for the average of the gauge factors over the
*     C* spatial directions,
*        pl(x) = sum_{k is C*}{ exp(i sign qhat Phi_k(x)/2 ) } pk(x),
*     and assigns the result to pl[ix]. For the definition of Phi_mu(x) see the
*     notes below. If iflag=0, sign=1.0 while sign=-1.0 if iflag=1. If coulomb=1
*     the gauge field used in the calculation of the phases contains the Coulomb
*     term A_mu^C(x). qhat is the electric charge in units of qel of the active
*     flavour, i.e. qhat=dirac_parms().qhat.
*
*
* Notes:
*
* The phase Phi_mu(ix) is defined in terms of the U(1) gauge field
*
*    A_mu(x) <-- adfld()[link_offset(ix,mu)]
*
* and the "Coulomb" gauge field
*
*    A_mu^C(x) = qel/nabla^2 sum_{k=1,2,3} nabla_k Ahat_{mu k}(x)
*
* where
*
*    nabla_k f(x) = 1/2 { f(x+k) - f(x-k) },
*
*    nabla^2 f(x) = sum_{k=1,2,3}{ f(x+k) + f(x-k) -2 f(x)},
*
* and Ahat_{mu nu}(x) is the U(1) field tensor computed by the u1ftensor()
* routine.
*
* If the flag coulomb=0, the phase is given by
*
*    Phi_k(x) = sum_{s=0}^{L_mu} A_mu(x+s mu),
*
* while, if coulomb=1, the phase is given by
*
*    Phi_k(x) = sum_{s=0}^{L_mu} { A_mu(x+s mu) + A_mu^C(x+s mu) }.
*
* The programs in this module act globally and must be called simultaneously on
* all MPI processes with the same translation vector.
*
*******************************************************************************/


#define CFACTOR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "lattice.h"
#include "flags.h"
#include "u1flds.h"
#include "cstates.h"
#include "global.h"

static int init=0,coulflag=-1;
static int lsize[4]={L0,L1,L2,L3};
static int psize[4]={-1,NPROC1/2,NPROC2,NPROC3};
static int faces[4]={FACE0,FACE1,FACE2,FACE3};
static int offset[4]={-1,0,(VOLUME+BNDRY),2*(VOLUME+BNDRY)};
static int *ix0b,np[4];
static double *lines;


static void set_lines(int coulomb)
{
   int mu,i,j,ix,iy,iz;
   int x[4],xm[4];
   double *ad,*ac=NULL;

   if (init==0)
   {
      np[0]=np[1]=np[2]=np[3]=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
      if ((bc_cstar()>=2)&&(cpr[1]>=NPROC1/2))
         np[2]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+NPROC2+cpr[3])&0x1;
      if ((bc_cstar()>=3)&&(cpr[1]>=NPROC1/2))
         np[3]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+cpr[3]+NPROC3)&0x1;

      lines=amalloc(3*(VOLUME+BNDRY)*sizeof(*lines),ALIGN);
      error(lines==NULL,1,"set_lines [cfactor.c]",
            "Unable to allocate service array");

      ix0b=amalloc(3*VOLUME*sizeof(*ix0b),ALIGN);
      error(ix0b==NULL,1,"set_lines [cfactor.c]",
            "Unable to allocate service array");

      init=1;
   }

   coulflag=coulomb;

   ad=adfld();
   copy_bnd_ad();

   if(coulflag)
      ac=get_coulomb_amu();

   for (mu=1;mu<=bc_cstar();++mu)
   {
      for(i=0;i<4;++i)
         xm[i]=lsize[i];
      xm[mu]=1;

      i=0;
      for (x[0]=0;x[0]<xm[0];++x[0])
      for (x[1]=0;x[1]<xm[1];++x[1])
      for (x[2]=0;x[2]<xm[2];++x[2])
      for (x[3]=0;x[3]<xm[3];++x[3])
      {
         ix=ipt[x[3]+L3*x[2]+L2*L3*x[1]+L1*L2*L3*x[0]];

         ix0b[ix+(mu-1)*VOLUME]=i;
         lines[ix+offset[mu]]=0.0;

         iz=ix;
         for(j=1;j<lsize[mu];++j)
         {
            iy=iup[iz][mu];

            ix0b[iy+(mu-1)*VOLUME]=i;
            lines[iy+offset[mu]]=lines[iz+offset[mu]]+2.0*ad[link_offset(iz,mu)];
            if(coulflag)
               lines[iy+offset[mu]]+=2.0*ac[iz+mu*(VOLUME+BNDRY)];

            iz=iy;
         }
         lines[i+VOLUME+offset[mu]]=0.5*lines[iz+offset[mu]]+ad[link_offset(iz,mu)];
         if(coulflag)
            lines[i+VOLUME+offset[mu]]+=ac[iz+mu*(VOLUME+BNDRY)];
         ++i;
      }
   }
}


static void add_next_long_lines(void)
{
   int il,mu,ix;
   int tag,saddr,raddr;
   double *sbuf,*rbuf;
   MPI_Status stat;

   for (mu=1;mu<=bc_cstar();++mu)
   {
      for (il=0;il<psize[mu];++il)
      {
         sbuf=lines+VOLUME+offset[mu];
         rbuf=sbuf+faces[mu];

         saddr=npr[2*mu+1];
         raddr=npr[2*mu];

         tag=mpi_tag();
         if (np[mu]==0)
         {
            MPI_Send(sbuf,faces[mu],MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,faces[mu],MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,faces[mu],MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,faces[mu],MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }

         for (ix=0;ix<faces[mu];++ix)
            sbuf[ix]=rbuf[ix];

         for (ix=0;ix<VOLUME;++ix)
            lines[ix+offset[mu]]+=lines[ix0b[ix+(mu-1)*VOLUME]+VOLUME+offset[mu]];
      }
   }
}


double* get_cfactor_phase(int coulomb,int mu)
{
   static double *cad=NULL;

   check_global_int("get_cfactor_phase",1,mu);
   check_global_int("get_cfactor_phase",1,coulomb);

   error(bc_cstar()<mu,1,"get_cfactor_phase [cfactor.c]",
         "invalid choice of mu: mu must be spatial and correspond to a C* direction");

   error(coulomb<0 || coulomb>1,1,"get_cfactor_phase [cfactor.c]",
         "invalid choice of the coulomb flag: must be 0 or 1");

   set_lines(coulomb);
   add_next_long_lines();
   cad=lines+offset[mu];

   return cad;
}

#define cmpl_mul(a,b,c) \
(a).re=(b).re*(c).re-(b).im*(c).im; \
(a).im=(b).re*(c).im+(b).im*(c).re


void mul_cfactor(int iflag,int coulomb,int mu,spinor_dble *pk,spinor_dble *pl)
{
   int ix,i,q;
   double *arg,d;
   complex_dble zx,sc,*prc,*psc;

   check_global_int("mul_cfactor",1,iflag);

   q=dirac_parms().qhat;

   error(q%2!=0,1,"mul_cfactor [cfactor.c]",
         "The electric charge of the spinor is not a multiple of two");

   error(iflag<0 || iflag>1,1,"mul_cfactor [cfactor.c]",
         "iflag must be 0 or 1");

   if (iflag)
      q=-q/2;
   else
      q=q/2;

   arg=get_cfactor_phase(coulomb,mu);

   for (ix=0;ix<VOLUME;ix++)
   {
      d=q*arg[ix];

      zx.re=cos(d);
      zx.im=sin(d);

      psc=(complex_dble*)(pk+ix);
      prc=(complex_dble*)(pl+ix);
      for(i=0;i<12;i++)
      {
         sc=psc[i];
         cmpl_mul(prc[i],zx,sc);
      }
   }
}


void mul_cfactor_muaverage(int iflag,int coulomb,spinor_dble *pk,spinor_dble *pl)
{
   int ix,i,q,mu;
   double *arg,phi;
   complex_dble zx,sc,*prc,*psc;

   check_global_int("mul_cfactor_muaverage",1,iflag);

   q=dirac_parms().qhat;

   error(q%2!=0,1,"mul_cfactor_muaverage [cfactor.c]",
         "The electric charge of the spinor is not a multiple of two");

   error(iflag<0 || iflag>1,1,"mul_cfactor [cfactor.c]",
         "iflag must be 0 or 1");

   if (iflag)
      q=-q/2;
   else
      q=q/2;

   arg=get_cfactor_phase(coulomb,1);

   for (ix=0;ix<VOLUME;ix++)
   {
      zx.re=0.0;
      zx.im=0.0;
      for(mu=1;mu<=bc_cstar();++mu)
      {
         phi=q*arg[ix+offset[mu]];
         zx.re+=cos(phi);
         zx.im+=sin(phi);
      }
      psc=(complex_dble*)(pk+ix);
      prc=(complex_dble*)(pl+ix);
      for(i=0;i<12;i++)
      {
         sc=psc[i];
         cmpl_mul(prc[i],zx,sc);
      }
   }
}

#undef cmpl_mul

