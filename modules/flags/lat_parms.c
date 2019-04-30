
/*******************************************************************************
*
* File lat_parms.c
*
* Copyright (C) 2009-2013, 2016 Martin Luescher, Isabel Campos
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Lattice parameters and boundary conditions.
*
* The externally accessible functions are
*
*   flds_parms_t set_flds_parms(int gauge,int nfl)
*     Sets the parameters which specify which fields need to be activated. The
*     parameters are
*
*       gauge          Active gauge fields.
*                      1 : SU(3)
*                      2 : U(1)
*                      3 : SU(3)xU(1)
*
*       nfl            Number of quark parameter sets (both valence and sea).
*
*     The return value is a structure that contains these parameters.
*
*   flds_parms_t flds_parms(void)
*     Returns the current field parameters in a structure that contains
*     the above parameters.
*
*   su3lat_parms_t set_su3lat_parms(double beta,double c0,
*                                   double cG,double cG_prime,int SFtype)
*     Sets the parameters for the SU(3) gauge action. The parameters are
*
*       beta           Inverse SU(3) bare coupling (beta=6/g0^2).
*
*       c0             Coefficient of the plaquette loops in the SU(3) gauge
*                      action (see doc/gauge_action.pdf).
*
*       cG,            SU(3) gauge action improvement coefficients at time 0
*       cG_prime       and T, respectively.
*
*       SFtype         Type of SF boundary for Luescher-Weisz actions
*                      (0: orbifold constuction - 1: Aoki-Frezzotti-Weisz)
*
*   su3lat_parms_t su3lat_parms(void)
*     Returns the current SU(3) gauge action parameters in a structure that
*     contains the above parameters.
*
*   u1lat_parms_t set_u1lat_parms(int type,double alpha,double invqel,
*                                 double lambda,double c0,
*                                 double cG,double cG_prime,int SFtype)
*     Sets the parameters for the U(1) gauge action. The parameters are
*
*       type           Type of U(1) gauge action (0: compact, 1: non-compact)
*
*       alpha          Bare fine-structure constant (alpha=1/(4 pi e0^2)).
*
*       invqel         Inverse of the elementary electric charge. Only matter
*                      with an electric charge equal to an integer value of
*                      1.0/invqel is allowed.
*
*       lambda         Gauge-fixing parameter (relevant only for the noncompact
*                      formulation).
*
*       c0             Coefficient of the plaquette loops in the compact U(1)
*                      gauge action (relevant only for the compact formulation).
*
*       cG,            U(1) gauge action improvement coefficients at time 0
*       cG_prime       and T, respectively (relevant only for the compact
*                      formulation).
*
*       SFtype         Type of SF boundary for Luescher-Weisz actions
*                      (0: orbifold constuction - 1: Aoki-Frezzotti-Weisz)
*
*   u1lat_parms_t u1lat_parms(void)
*     Returns the current U(1) gauge action parameters in a structure that
*     contains the above parameters.
*
*   dirac_parms_t set_qlat_parms(int ifl,double kappa,
*                                int qhat,double su3csw,
*                                double u1csw,double cF,double cF_prime,
*                                double th1,double th2,double th3)
*     Sets the Dirac-operator parameters for the inl-th quark flavour. The
*     parameters are
*
*       ifl            Number of quark flavour to be set.
*
*       kappa          Hopping parameter.
*
*       qhat           Integer electric charge (in units of the elementary one).
*
*       su3csw         Coefficients of the SU(3) Sheikholeslami-Wohlert term.
*
*       u1csw          Coefficients of the U(1) Sheikholeslami-Wohlert term.
*
*       cF,cF_prime    Fermion action improvement coefficients at time 0 and T,
*                      respectively.
*
*       th1,th2,th3    Angles specifying the phase-periodic boundary conditions
*                      for the quark fields in direction 1,2,3.       
*
*     The return value is a structure that contains the above parameters, but
*     kappa which is replaced by the bare quark mass m0.
*
*   dirac_parms_t qlat_parms(int ifl)
*     Returns the Dirac-operator parameters for the ifl-th quark flavour.
*
*   void print_lat_parms(void)
*     Prints the lattice parameters to stdout on MPI process 0.
*
*   void write_flds_bc_lat_parms(FILE *fdat)
*     Writes the global lattices sizes, fields, boundary and lattice parameters
*     to the file fdat on MPI process 0.
*
*   void check_flds_bc_lat_parms(FILE *fdat)
*     Compares the global lattice sizes, fields, boundary and the lattice
*     parameters with the values stored on the file fdat on MPI process 0,
*     assuming the latter were written to the file by the program
*     write_flds_bc_lat_parms().
*
*   bc_parms_t set_bc_parms(int type,int cstar,
*                           double *phi3,double *phi3_prime,
*                           double phi1,double phi1_prime)
*     Sets the boundary conditions and the associated parameters of the
*     action. The parameters are
*
*       type           Chosen type of time boundary condition (0: open, 1: SF,
*                      2: open-SF, 3: periodic).
*
*       cstar          Number of spatial directions with C-periodic boundary
*                      conditions.
*
*       phi3[0],       First two angles that define the boundary values of
*       phi3[1]        the SU(3) gauge field at time 0.
*
*       phi3_prime[0], First two angles that define the boundary values of
*       phi3_prime[1]  the SU(3) gauge field at time T.
*
*       phi1           Angle that defines the boundary values of the U(1)
*                      gauge field at time 0.
*
*       phi1_prime     Angle that defines the boundary values of the U(1)
*                      gauge field at time T.
*
*     The return value is a structure that contains these parameters plus
*     the third SU(3) angles. In this structure, the angles are stored in the
*     form of arrays phi3[2][3], where phi3[0][] and phi3[1][] are the
*     parameters at time 0 and T, respectively.
*
*   bc_parms_t bc_parms(void)
*     Returns a structure that contains the boundary parameters.
*
*   void print_bc_parms(void)
*     Prints the boundary parameters to stdout on MPI process 0.
*
*   int bc_type(void)
*     Returns the type of the chosen boundary conditions (0: open, 1: SF,
*     2: open-SF, 3: periodic).
*
*   void read_bc_parms(void)
*     On process 0, this program reads the following section from the stdin,
*     as explained in more details in doc/parms.pdf
*
*       [Boundary conditions]
*       type      <string|int>
*       cstar     <int>
*       su3phi    <double> <double>
*       su3phi'   <double> <double>
*       u1phi     <double>
*       u1phi'    <double>
*
*     This program assumes that 'set_flds_parms' has been already called.
*     Depending on the the gauge group selected with 'set_flds_parms', some
*     lines are not read and can be omitted in the input file. After reading
*     the stdin, the boundary conditions parameters are set with 'set_bc_parms'.
*
*   void read_glat_parms(void)
*     On process 0, this program reads a number of sections from the stdin
*     with the parameters of the gauge actions, as explained in more details
*     in doc/parms.pdf. This program assumes that 'set_flds_parms' and
*     'set_bc_parms' have been already called. If the SU(3) gauge field is
*     active, the program reads the following section
*
*       [SU(3) action]
*       beta         <double>
*       c0           <double>
*       SFtype       <string|int>
*       cG           <double>
*       cG'          <double>
*
*     Depending on the boundary conditions, some lines are not read and can be
*     omitted in the input file. After reading the stdin, the SU(3) action
*     parameters are set with 'set_su3lat_parms'.
*
*     If the U(1) gauge field is active, the program reads the following section
*    
*       [U(1) action]
*       type         <string|int>
*       alpha        <double>
*       invqel       <double>
*       c0           <double>
*       SFtype       <string|int>
*       cG           <double>
*       cG'          <double>
*
*     Depending on the boundary conditions, some lines are not read and can be
*     omitted in the input file. After reading the stdin, the U(1) action
*     parameters are set with 'set_u1lat_parms'.
*
*   void read_qlat_parms(void)
*     On process 0, this program reads a number of sections from the stdin
*     with the parameters of the quark actions, as explained in more details
*     in doc/parms.pdf. This program assumes that 'set_flds_parms' and
*     'set_bc_parms' have been already called.
*
*     If nfl is the number of quark flavours set with the 'set_flds_parms'
*     program, for each ifl=0,...,nfl-1, the program reads the following
*     sections
*       
*       [Flavour ifl]
*       qhat         <int>
*       kappa        <double>
*       su3csw       <double>
*       u1csw        <double>
*       cF           <double>
*       cF'          <double>
*       theta        <double>
*
*     Depending on the boundary conditions and the active gauge, some lines are
*     not read and can be omitted in the input file. After reading the stdin,
*     the quark flavour parameters are set with 'set_qlat_parms'.
*
* Notes:
*
* To ensure the consistency of the data base, the parameters must be set
* simultaneously on all processes. The data types flds_parms_t,su3lat_parms_t,
* u1lat_parms_t,dirac_parms_t,bc_parms_t are defined in the file flags.h.
*
* The programs set_flds_parms(), set_su3lat_parms(), set_u1lat_parms() and
* set_bc_parms() may be called at most once. set_bc_parms() must be called
* before set_su3lat_parms(), set_u1lat_parms() and set_qlat_parms(),
* Moreover, set_su3lat_parms(), set_u1lat_parms() and set_qlat_parms() may not
* be called after the geometry arrays are set up.
*
* The default values for the field parameters are
*       gauge          1 : SU(3)
*       nfl            0
* The default values for the SU(3) gauge action parameters are
*       beta           0.0
*       c0             1.0
*       cG=cG'         1.0
*       SFtype         0: orbifold construction
* The default values for the U(1) gauge action parameters are
*       type           0: compact
*       alpha          DBL_MAX
*       invqel         1.0
*       beta           0.0
*       lambda         0.0
*       c0             1.0
*       cG=cG'         1.0
*       SFtype         0: orbifold construction
* The default values for the quark parameters are
*       m0             DBL_MAX
*       qhat           0
*       su3csw         0.0
*       u1csw          0.0
*       cF=cF'         1.0
*       th1=th2=th3    0.0
* The default values for the boundary condition parameters are
*       type           0: open
*       cstar          0
*       su3phi=su3phi' 0.0
*       u1phi=u1phi'   0.0
*
* See the notes doc/gauge_action.pdf and doc/dirac.pdf for the detailed
* description of the lattice action and the boundary conditions.
*
*******************************************************************************/

#define LAT_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int flg_flds=0,flg_su3lat=0,flg_u1lat=0,flg_bc=0,*flg_qlat=NULL;
static flds_parms_t flds={1,0};
static su3lat_parms_t su3lat={0,0.0,1.0,0.0,{1.0,1.0}};
static u1lat_parms_t u1lat={0,0,DBL_MAX,1.0,0.0,0.0,1.0,0.0,{1.0,1.0}};
static double *qkappa=NULL;
static dirac_parms_t *qlat=NULL;
static bc_parms_t bc={0,0,{{0.0,0.0,0.0},{0.0,0.0,0.0}},{0.0,0.0}};


flds_parms_t set_flds_parms(int gauge,int nfl)
{
   int ifl;
   
   error(flg_flds!=0,1,"set_flds_parms [lat_parms.c]",
         "Attempt to reset the field parameters");

   error(flg_bc!=0,1,"set_flds_parms [lat_parms.c]",
         "Field parameters must be set before boundary conditions");

   error((flg_su3lat!=0)||(flg_u1lat!=0)||(flg_qlat!=NULL),1,
         "set_flds_parms [lat_parms.c]",
         "Field parameters must be set before lattice parameters");

   check_global_int("set_flds_parms",2,gauge,nfl);

   error_root((gauge<1)||(gauge>3),1,"set_flds_parms [lat_parms.c]",
              "Unknown gauge group");

   error_root(nfl<0,1,"set_flds_parms [lat_parms.c]",
              "Number of quark parameter sets must be non-negative");
   
   flds.gauge=gauge;
   flds.nfl=nfl;
   
   if (nfl!=0)
   {
      qkappa=malloc(nfl*sizeof(double));
      qlat=malloc(nfl*sizeof(dirac_parms_t));
      flg_qlat=malloc(nfl*sizeof(int));
      error((qkappa==NULL)||(qlat==NULL)||(flg_qlat==NULL),1,
            "set_flds_parms [lat_parms.c]",
            "Unable to allocate parameter array");
      for (ifl=0;ifl<nfl;ifl++)
      {
         flg_qlat[ifl]=0;
         qkappa[ifl]=0.0;
         qlat[ifl].qhat=0;
         qlat[ifl].m0=DBL_MAX;
         qlat[ifl].su3csw=0.0;
         qlat[ifl].u1csw=0.0;
         qlat[ifl].cF[0]=1.0;
         qlat[ifl].cF[1]=1.0;
         qlat[ifl].theta[0]=0.0;
         qlat[ifl].theta[1]=0.0;
         qlat[ifl].theta[2]=0.0;
      }
   }
   
   flg_flds=1;
   
   return flds;
}


bc_parms_t set_bc_parms(int type,int cstar,double *phi3,double *phi3_prime,double phi1, double phi1_prime)
{
   error(flg_bc!=0,1,"set_bc_parms [lat_parms.c]",
         "Attempt to reset the boundary conditions");

   error(iup[0][0]!=0,1,"set_bc_parms [lat_parms.c]",
         "Geometry arrays are already set");


   check_global_int("set_bc_parms",2,type,cstar);
   check_global_dble("set_bc_parms",6,phi3[0],phi3[1],phi3_prime[0],phi3_prime[1],
                                      phi1,phi1_prime);

   error_root((type<0)||(type>3),1,"set_bc_parms [lat_parms.c]",
              "Unknown type of time boundary condition");

   error_root((cstar<0)||(cstar>3),1,"set_bc_parms [lat_parms.c]",
              "Invalid number of directions with Cstar boundary conditions");

   bc.type=type;
   bc.cstar=cstar;

   if (type==1)
   {
      if (cstar==0)
      {
         if ((flds.gauge&1)!=0)
         {
            bc.phi3[0][0]=phi3[0];
            bc.phi3[0][1]=phi3[1];
            bc.phi3[0][2]=-phi3[0]-phi3[1];

            bc.phi3[1][0]=phi3_prime[0];
            bc.phi3[1][1]=phi3_prime[1];
            bc.phi3[1][2]=-phi3_prime[0]-phi3_prime[1];
         }
         if ((flds.gauge&2)!=0)
         {
            bc.phi1[0]=phi1;
            bc.phi1[1]=phi1_prime;
         }
      }
   }
   else if (type==2)
   {
      if(cstar==0)
      {
         if ((flds.gauge&1)!=0)
         {
            bc.phi3[1][0]=phi3_prime[0];
            bc.phi3[1][1]=phi3_prime[1];
            bc.phi3[1][2]=-phi3_prime[0]-phi3_prime[1];
         }
         if ((flds.gauge&2)!=0)
         {
            bc.phi1[1]=phi1_prime;
         }
      }
   }

   flg_bc=1;

   return bc;
}


su3lat_parms_t set_su3lat_parms(double beta,double c0,
                                double cG,double cG_prime,int SFtype)
{
   error(flg_flds==0,1,"set_su3lat_parms [lat_parms.c]",
         "Field parameters must be set first");

   error(flg_su3lat!=0,1,"set_su3lat_parms [lat_parms.c]",
         "Attempt to reset the SU(3) lattice parameters");

   error(iup[0][0]!=0,1,"set_su3lat_parms [lat_parms.c]",
         "Geometry arrays are already set");

   error(flg_bc==0,1,"set_su3lat_parms [lat_parms.c]",
         "Boundary conditions must be set first");

   check_global_int("set_bc_parms",1,SFtype);
   check_global_dble("set_su3lat_parms",4,beta,c0,cG,cG_prime);

   error_root(beta<=0.0,1,"set_su3lat_parms [lat_parms.c]",
              "Parameter beta must be positive");

   error_root(c0<=0.0,1,"set_su3lat_parms [lat_parms.c]",
              "Parameter c0 must be positive");

   error_root((SFtype<0)||(SFtype>1),1,"set_su3lat_parms [lat_parms.c]",
              "Unknown type of SF boundary condition");
   

   if ((flds.gauge&1)!=0)
   {
      su3lat.beta=beta;
      su3lat.c0=c0;
      su3lat.c1=0.125*(1.0-c0);
      if ((bc.type>=0)&&(bc.type<3))
      {
         su3lat.cG[0]=cG;
         if (bc.type==0)
            su3lat.cG[1]=cG;
         else if (bc.type==1)
            su3lat.cG[1]=cG;
         else if (bc.type==2)
            su3lat.cG[1]=cG_prime;
      }
      if ((bc.type>=1)&&(bc.type<=2)&&(c0!=1.0))
      {
         su3lat.SFtype=SFtype;
      }

      flg_su3lat=1;
   }

   return su3lat;
}


u1lat_parms_t set_u1lat_parms(int type,double alpha,double invqel,double lambda,
                              double c0,double cG,double cG_prime,int SFtype)
{
   error(flg_flds==0,1,"set_u1lat_parms [lat_parms.c]",
         "Field parameters must be set first");

   error(flg_u1lat!=0,1,"set_u1lat_parms [lat_parms.c]",
         "Attempt to reset the U(1) lattice parameters");

   error(flg_bc==0,1,"set_u1lat_parms [lat_parms.c]",
         "Boundary conditions must be set first");

   check_global_int("set_u1lat_parms",2,type,SFtype);
   check_global_dble("set_u1lat_parms",6,alpha,invqel,lambda,c0,cG,cG_prime);

   error_root((type<0)||(type>1),1,"set_u1lat_parms [lat_parms.c]",
              "Parameter type must be either 0 or 1");

   error_root(alpha<=0.0,1,"set_u1lat_parms [lat_parms.c]",
              "Parameter alpha must be positive");

   error_root(c0<=0.0,1,"set_u1lat_parms [lat_parms.c]",
              "Parameter c0 must be positive");

   error_root(invqel<=0.0,1,"set_u1lat_parms [lat_parms.c]",
              "Parameter invqel must be positive");

   error_root(lambda<0.0,1,"set_u1lat_parms [lat_parms.c]",
              "Parameter lambda must be positive");

   error_root((SFtype<0)||(SFtype>1),1,"set_u1lat_parms [lat_parms.c]",
              "Unknown type of SF boundary condition");
   

   if ((flds.gauge&2)!=0)
   {
      u1lat.type=type;
      u1lat.alpha=alpha;
      u1lat.invqel=invqel;
      u1lat.beta=invqel*invqel/(16.*atan(1)*alpha);
      if(type==1)
      {
         u1lat.lambda=lambda;
      }
      else
      {
         u1lat.c0=c0;
         u1lat.c1=0.125*(1.0-c0);
         if ((bc.type>=0)&&(bc.type<3))
         {
            u1lat.cG[0]=cG;
            if (bc.type==0)
               u1lat.cG[1]=cG;
            else if (bc.type==1)
               u1lat.cG[1]=cG;
            else if (bc.type==2)
               u1lat.cG[1]=cG_prime;
         }
      }
      if ((bc.type>=1)&&(bc.type<=2)&&(type==0)&&(c0!=1.0))
      {
         u1lat.SFtype=SFtype;
      }

      flg_u1lat=1;
   }

   return u1lat;
}



dirac_parms_t set_qlat_parms(int ifl,double kappa,int qhat,
                           double su3csw,double u1csw,double cF,double cF_prime,
                           double th1,double th2,double th3)
{
   error(flg_flds==0,1,"set_qlat_parms [lat_parms.c]",
         "Field parameters must be set first");

   error(iup[0][0]!=0,1,"set_qlat_parms [lat_parms.c]",
         "Geometry arrays are already set");

   error(flg_bc==0,1,"set_qlat_parms [lat_parms.c]",
         "Boundary conditions must be set first");

   error_root(flds.nfl==0,1,"set_qlat_parms [lat_parms.c]",
              "No fermion fields have been requested in field parameters");
   

   check_global_int("set_qlat_parms",2,ifl,qhat);

   error_root((ifl<0)||(ifl>=flds.nfl),1,"set_qlat_parms [lat_parms.c]",
              "Parameter ifl is out of range");

   error(flg_qlat[ifl]!=0,1,"set_qlat_parms [lat_parms.c]",
         "Attempt to reset already specified fermion lattice parameters");


   check_global_dble("set_qlat_parms",8,kappa,su3csw,u1csw,cF,cF_prime,
                                       th1,th2,th3);

   
   qkappa[ifl]=kappa;
   
   qlat[ifl].su3csw=0.0;
   if ((flds.gauge&1)!=0) qlat[ifl].su3csw=su3csw;
   
   qlat[ifl].qhat=0;
   qlat[ifl].u1csw=0.0;
   if ((flds.gauge&2)!=0)
   {
      qlat[ifl].qhat=qhat;
      qlat[ifl].u1csw=u1csw;
   }
   
   qlat[ifl].cF[0]=0.0;
   qlat[ifl].cF[1]=0.0;
   if ((bc.type>=0)&&(bc.type<3))
   {
      qlat[ifl].cF[0]=cF;
      if (bc.type==0)
         qlat[ifl].cF[1]=cF;
      else if (bc.type==1)
         qlat[ifl].cF[1]=cF;
      else if (bc.type==2)
         qlat[ifl].cF[1]=cF_prime;
   }
   
   qlat[ifl].theta[0]=0.0;
   qlat[ifl].theta[1]=0.0;
   qlat[ifl].theta[2]=0.0;
   if(bc.cstar==0)
   {
      qlat[ifl].theta[0]=th1;
      qlat[ifl].theta[1]=th2;
      qlat[ifl].theta[2]=th3;
   }

   qlat[ifl].m0=DBL_MAX;
   if (kappa!=0.0)
      qlat[ifl].m0=1.0/(2.0*kappa)-4.0;

   flg_qlat[ifl]=1;

   return qlat[ifl];
}




flds_parms_t flds_parms(void)
{
   return flds;
}


int gauge(void)
{
   return flds.gauge;
}


bc_parms_t bc_parms(void)
{
   return bc;
}


int bc_type(void)
{
   return bc.type;
}


int bc_cstar(void)
{
   return bc.cstar;
}


su3lat_parms_t su3lat_parms(void)
{
   return su3lat;
}


u1lat_parms_t u1lat_parms(void)
{
   return u1lat;
}


dirac_parms_t qlat_parms(int ifl)
{
   error_loc(flds.nfl==0,1,"qlat_parms [lat_parms.c]",
              "No fermion fields have been requested in field parameters");

   error_loc((ifl<0)||(ifl>=flds.nfl),1,"qlat_parms [lat_parms.c]",
              "Parameter ifl is out of range");

   return qlat[ifl];
}




void print_flds_parms(void)
{
   int my_rank;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      printf("Gauge group: ");
      if (flds.gauge==1)
         printf("SU(3)\n");
      else if (flds.gauge==2)
         printf("U(1)\n");
      else if (flds.gauge==3)
         printf("SU(3)xU(1)\n");
      
      printf("\n");
      fflush(stdout);
   }
}


void print_bc_parms(void)
{
   int my_rank,n[3],N[4];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   N[0]=N0;
   N[1]=N1;
   N[2]=N2;
   N[3]=N3;
   if (bc.type==0) N[0]-=1;
   if (bc.cstar!=0) N[1]/=2;

   if (my_rank==0)
   {
      if (bc.cstar==0)
      {
         printf("Periodic boundary conditions in space\n");
      }
      else
      {
         printf("C-periodic boundary conditions in %d spatial directions\n",bc.cstar);
      }
      
      if (bc.type==0)
      {
         printf("Open boundary conditions in time\n");
         printf("Physical lattice %dx%dx%dx%d\n",N[0],N[1],N[2],N[3]);
      }
      else if (bc.type==1)
      {
         printf("SF boundary conditions in time\n");
         printf("Physical lattice %dx%dx%dx%d\n",N[0],N[1],N[2],N[3]);

         if ((flds.gauge&1)!=0)
         {
            n[0]=fdigits(bc.phi3[0][0]);
            n[1]=fdigits(bc.phi3[0][1]);
            n[2]=fdigits(bc.phi3[0][2]);
            printf("su3phi = %.*f,%.*f,%.*f\n",IMAX(n[0],1),bc.phi3[0][0],
                   IMAX(n[1],1),bc.phi3[0][1],IMAX(n[2],1),bc.phi3[0][2]);

            n[0]=fdigits(bc.phi3[1][0]);
            n[1]=fdigits(bc.phi3[1][1]);
            n[2]=fdigits(bc.phi3[1][2]);
            printf("su3phi' = %.*f,%.*f,%.*f\n",IMAX(n[0],1),bc.phi3[1][0],
                   IMAX(n[1],1),bc.phi3[1][1],IMAX(n[2],1),bc.phi3[1][2]);
         }
         
         if ((flds.gauge&2)!=0)
         {
            n[0]=fdigits(bc.phi1[0]);
            printf("u1phi = %.*f\n",IMAX(n[0],1),bc.phi1[0]);
            
            n[0]=fdigits(bc.phi1[1]);
            printf("u1phi' = %.*f\n",IMAX(n[0],1),bc.phi1[1]);
         }
      }
      else if (bc.type==2)
      {
         printf("Open-SF boundary conditions in time\n");
         printf("Physical lattice %dx%dx%dx%d\n",N[0],N[1],N[2],N[3]);

         if ((flds.gauge&1)!=0)
         {
            n[0]=fdigits(bc.phi3[1][0]);
            n[1]=fdigits(bc.phi3[1][1]);
            n[2]=fdigits(bc.phi3[1][2]);
            printf("su3phi' = %.*f,%.*f,%.*f\n",IMAX(n[0],1),bc.phi3[1][0],
                   IMAX(n[1],1),bc.phi3[1][1],IMAX(n[2],1),bc.phi3[1][2]);
         }
         
         if ((flds.gauge&2)!=0)
         {
            n[0]=fdigits(bc.phi1[1]);
            printf("u1phi' = %.*f\n",IMAX(n[0],1),bc.phi1[1]);
         }
      }
      else
      {
         printf("Periodic boundary conditions in time\n");
         printf("Physical lattice %dx%dx%dx%d\n",N[0],N[1],N[2],N[3]);
      }

      printf("\n");
      fflush(stdout);
   }
}


void print_lat_parms(void)
{
   int my_rank,n,ifl;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      if ((flds.gauge&1)!=0) {
         printf("SU(3) action parameters:\n");
         n=fdigits(su3lat.beta);
         printf("beta = %.*f\n",IMAX(n,1),su3lat.beta);
         n=fdigits(su3lat.c0);
         printf("c0 = %.*f, ",IMAX(n,1),su3lat.c0);
         n=fdigits(su3lat.c1);
         printf("c1 = %.*f\n",IMAX(n,1),su3lat.c1);
         if (bc.type!=3) {
            n=fdigits(su3lat.cG[0]);
            printf("cG = %.*f\n",IMAX(n,1),su3lat.cG[0]);
         }
         if (bc.type==2) {
            n=fdigits(su3lat.cG[1]);
            printf("cG' = %.*f\n",IMAX(n,1),su3lat.cG[1]);
         }
         if (((bc.type==1)||(bc.type==2))&&(su3lat.c0!=1.0)) {
            if (su3lat.SFtype==0)
               printf("Orbifold-type SF boundary\n");
            else
               printf("AFW-typeB SF boundary\n");
         }
         printf("\n");
      }
      
      if ((flds.gauge&2)!=0) {
         printf("U(1) action parameters:\n");
         if(u1lat.type==1)
         {
            printf("Non-compact action\n");
         }
         else
         {
               printf("Compact action\n");
         }
         if (u1lat.alpha==DBL_MAX)
         {
            printf("alpha = infty\n");
         }
         else
         {
            n=fdigits(u1lat.alpha);
            printf("alpha = %.*f\n",IMAX(n,1),u1lat.alpha);
         }
         n=fdigits(u1lat.invqel);
         printf("qel = 1 / %.*f\n",IMAX(n,1),u1lat.invqel);
         if(u1lat.type==1)
         {
            n=fdigits(u1lat.lambda);
            printf("lambda = %.*f\n",IMAX(n,1),u1lat.lambda);
         }
         else
         {
            n=fdigits(u1lat.c0);
            printf("c0 = %.*f, ",IMAX(n,1),u1lat.c0);
            n=fdigits(u1lat.c1);
            printf("c1 = %.*f\n",IMAX(n,1),u1lat.c1);
            if (bc.type!=3) {
               n=fdigits(u1lat.cG[0]);
               printf("cG = %.*f\n",IMAX(n,1),u1lat.cG[0]);
            }
            if (bc.type==2) {
               n=fdigits(u1lat.cG[1]);
               printf("cG' = %.*f\n",IMAX(n,1),u1lat.cG[1]);
            }
            if (((bc.type==1)||(bc.type==2))&&(u1lat.c0!=1.0)) {
               if (u1lat.SFtype==0)
                  printf("Orbifold-type SF boundary\n");
               else
                  printf("AFW-typeB SF boundary\n");
            }
         }
         printf("\n");
      }
      
      if (flds.nfl>0)
         printf("Number of quark non-degenerate flavours: %d\n\n",flds.nfl);
      
      for (ifl=0;ifl<flds.nfl;ifl++)
      {
         printf("Flavour %d:\n",ifl);

         n=fdigits(qkappa[ifl]);
         printf("kappa = %.*f\n",IMAX(n,1),qkappa[ifl]);
         if ((flds.gauge&2)!=0)
         {
            printf("qhat = %d\n",qlat[ifl].qhat);
         }
         if ((flds.gauge&1)!=0)
         {
            n=fdigits(qlat[ifl].su3csw);
            printf("su3csw = %.*f\n",IMAX(n,1),qlat[ifl].su3csw);
         }
         if ((flds.gauge&2)!=0)
         {
            n=fdigits(qlat[ifl].u1csw);
            printf("u1csw = %.*f\n",IMAX(n,1),qlat[ifl].u1csw);
         }
         if (bc.type!=3)
         {
            n=fdigits(qlat[ifl].cF[0]);
            printf("cF = %.*f\n",IMAX(n,1),qlat[ifl].cF[0]);
         }
         if (bc.type==2)
         {
            n=fdigits(qlat[ifl].cF[1]);
            printf("cF' = %.*f\n",IMAX(n,1),qlat[ifl].cF[1]);
         }
         n=fdigits(qlat[ifl].theta[0]);
         printf("theta = %.*f,",IMAX(n,1),qlat[ifl].theta[0]);
         n=fdigits(qlat[ifl].theta[1]);
         printf("%.*f,",IMAX(n,1),qlat[ifl].theta[1]);
         n=fdigits(qlat[ifl].theta[2]);
         printf("%.*f\n",IMAX(n,1),qlat[ifl].theta[2]);

         printf("\n");
         
         fflush(stdout);
      }
   }
}


void print_flds_bc_lat_parms(void)
{
   print_flds_parms();
   print_bc_parms();
   print_lat_parms();
}




void write_flds_bc_lat_parms(FILE *fdat)
{
   int ifl;
   
   write_little_int(1,fdat,6,N0,N1,N2,N3,flds.gauge,flds.nfl);

   write_little_int(1,fdat,2,bc.type,bc.cstar);
   write_little_dble(1,fdat,8,bc.phi3[0][0],bc.phi3[0][1],bc.phi3[0][2],
                     bc.phi3[1][0],bc.phi3[1][1],bc.phi3[1][2],
                     bc.phi1[0],bc.phi1[1]);
   
   if ((flds.gauge&1)!=0)
   {
      write_little_int(1,fdat,1,su3lat.SFtype);
      write_little_dble(1,fdat,5,su3lat.beta,su3lat.c0,su3lat.c1,
                               su3lat.cG[0],su3lat.cG[1]);
   }

   if ((flds.gauge&2)!=0)
   {
      write_little_int(1,fdat,2,u1lat.type,u1lat.SFtype);
      write_little_dble(1,fdat,7,u1lat.alpha,u1lat.invqel,u1lat.lambda,
                               u1lat.c0,u1lat.c1,u1lat.cG[0],u1lat.cG[1]);
   }
   
   for (ifl=0;ifl<flds.nfl;ifl++)
   {
      write_little_int(1,fdat,1,qlat[ifl].qhat);
      write_little_dble(1,fdat,8,qkappa[ifl],
                      qlat[ifl].su3csw,qlat[ifl].u1csw,
                      qlat[ifl].cF[0],qlat[ifl].cF[1],
                      qlat[ifl].theta[0],qlat[ifl].theta[1],qlat[ifl].theta[2]);
   }
}


void check_flds_bc_lat_parms(FILE *fdat)
{
   int ifl;

   check_little_int("check_flds_bc_lat_parms",fdat,6,N0,N1,N2,N3,flds.gauge,flds.nfl);

   check_little_int("check_flds_bc_lat_parms",fdat,2,bc.type,bc.cstar);
   check_little_dble("check_flds_bc_lat_parms",fdat,8,
                     bc.phi3[0][0],bc.phi3[0][1],bc.phi3[0][2],
                     bc.phi3[1][0],bc.phi3[1][1],bc.phi3[1][2],
                     bc.phi1[0],bc.phi1[1]);
   
   if ((flds.gauge&1)!=0)
   {
      check_little_int("check_flds_bc_lat_parms",fdat,1,su3lat.SFtype);
      check_little_dble("check_flds_bc_lat_parms",fdat,5,su3lat.beta,su3lat.c0,su3lat.c1,
                               su3lat.cG[0],su3lat.cG[1]);
   }

   if ((flds.gauge&2)!=0)
   {
      check_little_int("check_flds_bc_lat_parms",fdat,2,u1lat.type,u1lat.SFtype);
      check_little_dble("check_flds_bc_lat_parms",fdat,7,u1lat.alpha,u1lat.invqel,
                               u1lat.lambda,u1lat.c0,u1lat.c1,
                               u1lat.cG[0],u1lat.cG[1]);
   }
   
   for (ifl=0;ifl<flds.nfl;ifl++)
   {
      check_little_int("check_flds_bc_lat_parms",fdat,1,qlat[ifl].qhat);
      check_little_dble("check_flds_bc_lat_parms",fdat,8,qkappa[ifl],
                      qlat[ifl].su3csw,qlat[ifl].u1csw,
                      qlat[ifl].cF[0],qlat[ifl].cF[1],
                      qlat[ifl].theta[0],qlat[ifl].theta[1],qlat[ifl].theta[2]);
   }
}




void read_bc_parms(void)
{
   int my_rank;
   int gg,bc,cs;
   double phi3[2],phi3_prime[2];
   double phi1,phi1_prime;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   error(flg_flds==0,1,"read_bc_parms [lat_parms.c]",
         "Field parameters must be set first");

   gg=gauge();
   
   if (my_rank==0)
   {
      find_section("Boundary conditions");
      read_line("type","%s",&line);
      bc=4;
      if ((strcmp(line,"open")==0)||(strcmp(line,"0")==0))
         bc=0;
      else if ((strcmp(line,"SF")==0)||(strcmp(line,"1")==0))
         bc=1;
      else if ((strcmp(line,"open-SF")==0)||(strcmp(line,"2")==0))
         bc=2;
      else if ((strcmp(line,"periodic")==0)||(strcmp(line,"3")==0))
         bc=3;
      else
         error_root(1,1,"read_bc_parms [main.c]",
                    "Unknown time boundary condition type %s",line);
      
      read_line("cstar","%d",&cs);

      phi3[0]=0.0;
      phi3[1]=0.0;
      phi3_prime[0]=0.0;
      phi3_prime[1]=0.0;
      phi1=0.0;
      phi1_prime=0.0;
      if ((cs==0)&&(((gg&1)!=0)||(flg_flds==0)))
      {
         if (bc==1)
            read_dprms("su3phi",2,phi3);
         if ((bc==1)||(bc==2))
            read_dprms("su3phi'",2,phi3_prime);
      }
      if ((cs==0)&&(((gg&2)!=0)||(flg_flds==0)))
      {
         if (bc==1)
            read_line("u1phi","%lf",&phi1);
         if ((bc==1)||(bc==2))
            read_line("u1phi'","%lf",&phi1_prime);
      }
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi3,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi3_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&phi1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&phi1_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_bc_parms(bc,cs,phi3,phi3_prime,phi1,phi1_prime);
}


void read_glat_parms(void)
{
   int my_rank;
   int gg,bc,sf;
   int type;
   double beta,c0,cG,cG_prime;
   double alpha,lambda,invqel;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   error(flg_flds==0,1,"read_glat_parms [lat_parms.c]",
         "Field parameters must be set first");

   error(flg_bc==0,1,"read_glat_parms [lat_parms.c]",
         "Boundary conditions must be set first");
   
   gg=gauge();
   bc=bc_type();
   
   if((gg&1)!=0)
   {
      if (my_rank==0)
      {
         find_section("SU(3) action");
         read_line("beta","%lf",&beta);
         read_line("c0","%lf",&c0);

         cG=1.0;
         cG_prime=1.0;
         if (bc!=3)
            read_line("cG","%lf",&cG);
         if (bc==2)
            read_line("cG'","%lf",&cG_prime);
         
         sf=0;
         if (((bc==1)||(bc==2))&&(c0!=1.0))
         {
            read_line("SFtype","%s",&line);
            if ((strcmp(line,"orbifold")==0)||
                (strcmp(line,"openQCD-1.4")==0)||(strcmp(line,"0")==0))
               sf=0;
            else if ((strcmp(line,"AFW-typeB")==0)||
                     (strcmp(line,"openQCD-1.2")==0)||(strcmp(line,"1")==0))
               sf=1;
            else
               error_root(1,1,"read_glat_parms [main.c]",
                          "Unknown SF type %s in SU(3) action",line);
         }
      }
      
      MPI_Bcast(&sf,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      set_su3lat_parms(beta,c0,cG,cG_prime,sf);
   }

   if((gg&2)!=0)
   {
      if (my_rank==0)
      {
         find_section("U(1) action");
         read_line("type","%s",line);

         type=0;
         if ((strcmp(line,"compact")==0)||(strcmp(line,"0")==0))
            type=0;
         else if ((strcmp(line,"non-compact")==0)||(strcmp(line,"1")==0))
            type=1;
         else
            error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                       "Unknown U(1) action type %s",line);
      
         read_line("alpha","%lf",&alpha);
         read_line("invqel","%lf",&invqel);

         lambda=c0=cG=cG_prime=0.0;
         if(type==0)
         {
            read_line("c0","%lf",&c0);
            if (bc!=3) read_line("cG","%lf",&cG);
            if (bc==2) read_line("cG'","%lf",&cG_prime);
         }
         else
         {
            read_line("lambda","%lf",&lambda);
         }

         sf=0;
         if (((bc==1)||(bc==2))&&(type==0)&&(c0!=1.0))
         {
            read_line("SFtype","%s",&line);
            if ((strcmp(line,"orbifold")==0)||
                (strcmp(line,"openQCD-1.4")==0)||(strcmp(line,"0")==0))
               sf=0;
            else if ((strcmp(line,"AFW-typeB")==0)||
                     (strcmp(line,"openQCD-1.2")==0)||(strcmp(line,"1")==0))
               sf=1;
            else
               error_root(1,1,"read_glat_parms [main.c]",
                          "Unknown SF type %s in U(1) action",line);
         }
      }
      
      MPI_Bcast(&type,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&sf,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&invqel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&lambda,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      set_u1lat_parms(type,alpha,invqel,lambda,c0,cG,cG_prime,sf);
   }
}


void read_qlat_parms(void)
{
   int my_rank;
   int gg,bc,cs;
   int ifl,qhat;
   double kappa,su3csw,u1csw,cF,cF_prime,th1,th2,th3;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   error(flg_flds==0,1,"read_qlat_parms [lat_parms.c]",
         "Field parameters must be set first");

   error(flg_bc==0,1,"read_qlat_parms [lat_parms.c]",
         "Boundary conditions must be set first");
   
   gg=gauge();
   bc=bc_type();
   cs=bc_cstar();
   
   for (ifl=0;ifl<flds.nfl;ifl++)
   {
      qhat=0;
      su3csw=u1csw=0.0;
      cF=cF_prime=0.0;
      th1=th2=th3=0.0;

      sprintf(line,"Flavour %d",ifl);
      if (my_rank==0)
      {
         find_section(line);
         read_line("kappa","%lf",&kappa);
         if ((gg&2)!=0)
         {
            read_line("qhat","%d",&qhat);
            read_line("u1csw","%lf",&u1csw);
         }
         if ((gg&1)!=0)
            read_line("su3csw","%lf",&su3csw);
         if (bc!=3) read_line("cF","%lf",&cF);
         if (bc==2) read_line("cF'","%lf",&cF_prime);
         if (cs==0) read_line("theta","%lf %lf %lf",&th1,&th2,&th3);
      }
      
      MPI_Bcast(&qhat,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&su3csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&u1csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&th1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&th2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&th3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      set_qlat_parms(ifl,kappa,qhat,su3csw,u1csw,cF,cF_prime,th1,th2,th3);
   }
}
