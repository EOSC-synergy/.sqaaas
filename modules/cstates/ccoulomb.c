
/*******************************************************************************
*
* File ccoulomb.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of the Coulomb U(1) gauge field A_mu^C(x)
*
* The externally accessible functions are
*
*   void nabla_sq_dvec(double *s, double *r)
*     Applies the spatial laplacian to a the double matter field s[ix]
*     and assigns the result to r[ix],
*        r(x) = sum_{k=1,2,3}[ s(x+k) + s(x-k) -2 s(x) ]
*     The variables r and s, that cannot be equal, must point to memory spaces
*     of (VOLUME+BNDRY) double numbers.
*
*   void div_sym_dvec(double *s1, double *s2, double *s3, double *r)
*     Applies the spatial symmetric divergence to the double matter fields
*     s1[ix], s2[ix], and s3[ix] and assigns the result to r[ix],
*        r(x) = sum_{k=1,2,3}[ s_k(x+k) + s_k(x-k) ]/2
*     The variables r and s1, s2 and s3 must point to memory spaces of
*     (VOLUME+BNDRY) double numbers. None of the variables s1, s2, s3 can be
*     equal to the result variable r.
*
*   double inv_nabla_sq(double *out, double *in)
*     Applies the inverse spatial laplacian to the double matter field in[ix]
*     and assigns the result to out[ix],
*        out(x) = 1/nabla_sq in(x)
*     by using the CG algorithm. The CG iteration is stopped when the relative
*     error is smaller than 10^(-8), and it throws an error if this target is
*     not achievable. The function returns an upper bound for the relative
*     residue.
*
*   double* get_coulomb_amu(void)
*     Returns the pointer amuc to the base address of a global real field that
*     contains the Coulomb A_mu^C(x), see the notes below. The real field is
*     organized in memory as a matter field so that amuc[ix] is equal to A_mu^C(x)
*     with ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0].
*
* Notes:
*
* The "Coulomb" gauge field is given by
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
* The programs in this module act globally and must be called simultaneously on
* all MPI processes with the same translation vector.
*
*******************************************************************************/


#define CCOULOMB_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "lattice.h"
#include "flags.h"
#include "u1flds.h"
#include "linalg.h"
#include "cstates.h"
#include "u1ftensor.h"
#include "global.h"


#ifdef M_PI
#undef M_PI
#endif
#define M_PI        3.14159265358979323846264338327950288


#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static double *camu=NULL;


static void pack_and_send_sfield(double *v)
{
   static int nfc[8],np[8];
   static int init=0;
   static double *vbuf=NULL;
   int ifc,ib,ix,eo,ofs[8];
   int saddr,raddr;
   int nbf,tag;
   double *sbuf,*rbuf;
   MPI_Status stat;

   if(init==0)
   {
      if (BNDRY!=0)
      {
         vbuf=amalloc(BNDRY*sizeof(*vbuf),4);
         error(vbuf==NULL,1,"pack_and_send_gbuf [ccoulomb.c]",
               "Unable to allocate memory space for vbuf");
      }

      nfc[0]=FACE0/2;
      nfc[1]=FACE0/2;
      nfc[2]=FACE1/2;
      nfc[3]=FACE1/2;
      nfc[4]=FACE2/2;
      nfc[5]=FACE2/2;
      nfc[6]=FACE3/2;
      nfc[7]=FACE3/2;

      np[0]=np[2]=np[4]=np[6]=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
      if ((bc_cstar()>=2)&&(cpr[1]>=NPROC1/2))
         np[4]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+NPROC2+cpr[3])&0x1;
      if ((bc_cstar()>=3)&&(cpr[1]>=NPROC1/2))
         np[6]=(cpr[0]+cpr[1]-NPROC1/2+cpr[2]+cpr[3]+NPROC3)&0x1;
      for(ifc=0;ifc<4;++ifc)
         np[2*ifc+1]=np[2*ifc];

      init=1;
   }

   for(eo=0;eo<2;++eo)
   {
      ofs[0]=eo*(BNDRY/2);
      ofs[1]=ofs[0]+nfc[0];
      ofs[2]=ofs[1]+nfc[1];
      ofs[3]=ofs[2]+nfc[2];
      ofs[4]=ofs[3]+nfc[3];
      ofs[5]=ofs[4]+nfc[4];
      ofs[6]=ofs[5]+nfc[5];
      ofs[7]=ofs[6]+nfc[6];

      for (ifc=2;ifc<8;ifc++)
      {
         for (ib=0;ib<nfc[ifc];ib++)
         {
            ix=map[ofs[ifc]+ib];
            vbuf[ofs[ifc]+ib]=v[ix];
         }
      }

      for (ifc=2;ifc<8;ifc++)
      {
         nbf=nfc[ifc];

         if (nbf>0)
         {
            tag=mpi_tag();
            saddr=npr[ifc^0x1];
            raddr=npr[ifc];
            sbuf=vbuf+ofs[ifc];
            rbuf=v+VOLUME+ofs[ifc];

            if (np[ifc]==0)
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
}


static void serv_nabla_sq_dvec(double *s, double *r)
{
   int mu,ix;

   pack_and_send_sfield(s);

   for (ix=0;ix<VOLUME;++ix)
   {
      r[ix]=-6.0*s[ix];
      for (mu=1;mu<4;++mu)
      {
         r[ix]+=s[iup[ix][mu]];
         r[ix]+=s[idn[ix][mu]];
      }
   }
}


void nabla_sq_dvec(double *s, double *r)
{
   error(s==r,1,"nabla_sq_dvec [ccoulomb.c]",
         "Input and output fields must be different");

   serv_nabla_sq_dvec(s,r);
}


static void serv_div_sym_dvec(double *s1, double *s2, double *s3, double *r)
{
   int ix;

   pack_and_send_sfield(s1);
   pack_and_send_sfield(s2);
   pack_and_send_sfield(s3);

   for (ix=0;ix<VOLUME;++ix)
   {
      r[ix] =s1[iup[ix][1]]-s1[idn[ix][1]];
      r[ix]+=s2[iup[ix][2]]-s2[idn[ix][2]];
      r[ix]+=s3[iup[ix][3]]-s3[idn[ix][3]];

      r[ix]*=0.5;
   }
}


void div_sym_dvec(double *s1, double *s2, double *s3, double *r)
{
   error(s1==r || s2==r || s3==r,1,"div_sym_dvec [ccoulomb.c]",
         "None of the input fields can be equal to the output field");

   serv_div_sym_dvec(s1,s2,s3,r);
}


static void lincomb_assign_dvec(int vol,double a,double *X,double b,double *Y)
{
   double *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      (*Y)=b*(*Y)+a*(*X);
      Y+=1;
   }
}


static void mul_assign_dvec(int vol,double a,double *X,double *Y)
{
   double *Xm;

   Xm=X+vol;

   for (;X<Xm;X++)
   {
      (*Y)=a*(*X);
      Y+=1;
   }
}


static void alloc_camu(void)
{
   if(camu==NULL)
   {
      camu=amalloc(9*(VOLUME+BNDRY)*sizeof(*camu),4);
      error(camu==NULL,1,"alloc_camu [ccoulomb.c]",
            "Unable to allocate camu array");
   }
}


double inv_nabla_sq(double *out, double *in)
{
   static double kappa,relerr,relres;
   static int niters, nmax=0, my_rank;
   
   double *p=NULL;
   double *r=NULL;
   double *Ap=NULL;
   double *tout=NULL;
   double rr, rrold, alpha, eps;
   int n;

   error(bc_cstar()==0,1,"inv_nabla_sq [ccoulomb.c]",
         "The spatial Laplace operator is not invertible without C* boundary conditions");

   if(camu==NULL)
      alloc_camu();
   
   if(nmax==0)
   {
      /*The asymmetry in the components comes from the orbifold construction*/
      kappa=sin(M_PI/(N1))*sin(M_PI/(N1));
      if(bc_cstar()>=2)
         kappa+=sin(M_PI/(2*N2))*sin(M_PI/(2*N2));
      if(bc_cstar()==3)
         kappa+=sin(M_PI/(2*N3))*sin(M_PI/(2*N3));
      
      kappa=3./kappa;
      
      relerr=1.e-8;
      
      rr=2./(1.+sqrt(kappa));
      if(rr<1.e-8)
      {
         niters=(int)( - log(relerr/(2.*kappa)) / rr ) + 1;
         nmax=(int)( - log(relerr/(2000.*kappa)) / rr )*2;
      }
      else
      {
         niters=(int)( log(relerr/(2.*kappa)) / log(1.-rr) ) + 1;
         nmax=(int)( log(relerr/(2000.*kappa)) / log(1.-rr) )*2;
      }
      
      relres=relerr/kappa;
      
      MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
      if(my_rank==0)
      {
         printf("\n"
                "The inversion of the spatial free Laplace operator is needed for the calculation\n"
                "of the Coulomb dressing factor. This is done at the moment with a CG (it would\n"
                "be better to use a Fourier representation).\n"
                "   Condition number = %.2e\n"
                "   Required relative error on solution (hardcoded) = %.2e\n"
                "   Required relative residue = %.2e\n"
                "   Expected number of iterations = %d\n"
                "   Maximum number of iterations = %d\n\n",
                kappa,relerr,relres,niters,nmax);
      }
      
      error(relres<1.e-13,1,"inv_nabla_sq [ccoulomb.c]",
            "The relative residues is smaller than 1e-13. Consider changing the required relative error (relerr variable in the code)");
   
      error(niters>100000,1,"inv_nabla_sq [ccoulomb.c]",
            "The expected number of iterations is larger than 1e+6. Consider changing the required relative error (relerr variable in the code)");
   }

   p   =camu+4*(VOLUME+BNDRY);
   r   =camu+5*(VOLUME+BNDRY);
   Ap  =camu+6*(VOLUME+BNDRY);
   tout=camu+7*(VOLUME+BNDRY);

   set_dvec2zero(VOLUME,tout);

   set_dvec2zero(VOLUME,r);
   muladd_assign_dvec(VOLUME,-1.0,in,r);
   assign_dvec2dvec(VOLUME,r,p);

   rr=norm_square_dvec(VOLUME,1,r);
   eps=relres*sqrt(rr);

   MPI_Bcast(&nmax,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   for(n=0;n<nmax;n++)
   {
      serv_nabla_sq_dvec(p,Ap);
      alpha=-rr/scalar_prod_dvec(VOLUME,1,p,Ap);
      muladd_assign_dvec(VOLUME,alpha,p,tout);
      muladd_assign_dvec(VOLUME,alpha,Ap,r);
      rrold=rr;
      rr=norm_square_dvec(VOLUME,1,r);
      if(sqrt(rr)<eps)
      {
         assign_dvec2dvec(VOLUME,tout,out);
         break;
      }

      lincomb_assign_dvec(VOLUME,1.0,r,rr/rrold,p);
   }

   error(n==nmax,1,"inv_nabla_sq [ccoulomb.c]",
         "Unable to invert nabla_sq in nmax iterations");


   return relres;
}


double* get_coulomb_amu(void)
{
   static double *ramu=NULL;
   double *div=NULL;
   double *s1=NULL;
   double *s2=NULL;
   double *s3=NULL;
   double qel,**ft;
   u1lat_parms_t lat;

   lat=u1lat_parms();
   qel=1.0/lat.invqel;

   if(camu==NULL)
      alloc_camu();

   div=camu+8*(VOLUME+BNDRY);
   s1 =camu+7*(VOLUME+BNDRY);
   s2 =camu+6*(VOLUME+BNDRY);
   s3 =camu+5*(VOLUME+BNDRY);

   ft=u1ftensor();

   /* mu=0 */

   mul_assign_dvec(VOLUME,qel,ft[0],s1);
   mul_assign_dvec(VOLUME,qel,ft[1],s2);
   mul_assign_dvec(VOLUME,qel,ft[2],s3);

   serv_div_sym_dvec(s1,s2,s3,div);
   inv_nabla_sq(camu,div);

   /* mu=1 */

   set_dvec2zero(VOLUME,s1);
   mul_assign_dvec(VOLUME,qel,ft[5],s2);
   mul_assign_dvec(VOLUME,-qel,ft[4],s3);

   serv_div_sym_dvec(s1,s2,s3,div);
   inv_nabla_sq(camu+(VOLUME+BNDRY),div);

   /* mu=2 */

   mul_assign_dvec(VOLUME,-qel,ft[5],s1);
   set_dvec2zero(VOLUME,s2);
   mul_assign_dvec(VOLUME,qel,ft[3],s3);

   serv_div_sym_dvec(s1,s2,s3,div);
   inv_nabla_sq(camu+2*(VOLUME+BNDRY),div);

   /* mu=3 */

   mul_assign_dvec(VOLUME,qel,ft[4],s1);
   mul_assign_dvec(VOLUME,-qel,ft[3],s2);
   set_dvec2zero(VOLUME,s3);

   serv_div_sym_dvec(s1,s2,s3,div);
   inv_nabla_sq(camu+3*(VOLUME+BNDRY),div);

   ramu=camu;

   return ramu;
}
