
/*******************************************************************************
*
* File rat_parms.c
*
* Copyright (C) 2012, 2013 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Rational function parameter data base
*
* The externally accessible functions are
*
*   rat_parms_t set_rat_parms(int irp,int num,int den,int degree,double *range,
*                             double mu0,double A,double *nu,double *mu,
*                             double delta,double *x)
*     Sets the parameters in the rational function parameter set number
*     irp and returns a structure containing them (see the notes).
*
*   rat_parms_t rat_parms(int irp)
*     Returns a structure containing the rational function parameter set
*     number irp (see the notes).
*
*   void read_rat_parms(int irp)
*     On process 0, this program scans stdin for a line starting with the
*     string "[Rational <int>]" (after any number of blanks), where <int> is
*     the integer value passed by the argument. An error occurs if no such
*     line or more than one is found. If a rational approximation for
*     (x+mu^2)^(-1/2) is needed, the input file section must look like
*
*       power   -1 2
*       degree  <int>
*       range   <double> <double>
*       mu      <double>
*
*     In this case the rational approzimation is calculated by the function
*     zolotarev() (see the ratfcts module directory), called by the
*     set_rat_parms function. If a rational approximation for a different
*     power is needed, the input file section must look like
*
*       power   <int> <int>
*       degree  <int>
*       range   <double> <double>
*       mu      <double>
*       A       <double>
*       delta   <double>
*       nu[0]    <double>
*       ...
*       nu[M]    <double>
*       mu[0]    <double>
*       ...
*       mu[M]    <double>
*       x[0]    <double>
*       ...
*       x[N]    <double>
*
*     where M=degree-1 and N=2*degree+1. A is the normalization of the
*     rational approximation, nu[i] are the square-roots of the (-1)*zeroes of
*     the numerator, mu[i] are the square-roots of the (-1)*zeroes of the
*     denominator, x[i] are the extrema of the error function, and delta is the
*     signed relative error calculated at the leftmost extremum.
*
*   void print_rat_parms(void)
*     Prints the defined rational function parameter sets to stdout on MPI
*     process 0.
*
*   void write_rat_parms(FILE *fdat)
*     Writes the defined rational function parameter sets to the file fdat 
*     on MPI process 0.
*
*   void check_rat_parms(FILE *fdat)
*     Compares the defined rational function parameter sets with those 
*     on the file fdat on MPI process 0, assuming the latter were written
*     to the file by the program write_rat_parms().
*
* Notes:
*
* Rational approximations of degree [n,n] for general rational powers are
* supported. If (x+mu0^2)^(-1/2) is requested, then the parameters of the rational
* approximation are calculated at run-time in terms of Zolotorev rational
* functions (see the modules ratfcts/zolotarev.c and ratfcts/ratfcts.c). If a
* different powers is requested, the parameters of the rational approximations
* are read from a file (see the module ratfcts/ratfcts.c), which can be
* generated with the MinMax program (<root>/minmax).
*
* The elements of a structure of type rat_parms_t are
*
*   power[2]     Numerator and denominator of the power to be approximated
*
*   degree       Degree of the rational function
*
*   range[2]     Lower and upper end of the spectrum of gamma_5*D (or
*                gamma5*Dhat). Notice that the range of the rational
*                approximation is [ range[0]^2 , range[1]^2 ].
*
*   A            Normalization of the rational approximation
*
*   mu0          shift parameter, (x+mu0^2)^(power[0]/power[1]) 
*
*   nu           Square-roots of the (-1)*zeroes of the numerator
*
*   mu           Square-roots of the (-1)*zeroes of the denominator
*
*   delta        Relative error of the rational approximation
*
* Up to 32 parameter sets, labeled by an index irp=0,1,..,31, can be
* specified. Once a set is defined, it cannot be changed by calling
* set_rat_parms() again. Rational function parameters must be globally
* the same.
*
* Except for rat_parms(), the programs in this module perform global
* operations and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define RAT_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "ratfcts.h"
#include "global.h"

#define IRPMAX 32

static int init=0;
static rat_parms_t rp[IRPMAX+1]={{{-1,2},0,{0.0,0.0},0.0,0.0,0.0,NULL,NULL}};


static void init_rp(void)
{
   int irp;
   
   for (irp=1;irp<=IRPMAX;irp++)
      rp[irp]=rp[0];

   init=1;
}


static double rat_delta(double n,double A,double *nu,double *mu,double p,double x,double mu0)
{
   int i;
   double delta;
   
   delta=A*pow(x+mu0*mu0,-p);
   for (i=0;i<n;i++)
      delta*=(x+nu[i]*nu[i])/(x+mu[i]*mu[i]);
   
   delta=1.0-delta;
   return delta;
}


rat_parms_t set_rat_parms(int irp,int num,int den,int degree,double *range,
                          double mu0,double A,double *nu,double *mu,
                          double delta,double *x)
{
   int ie,i;
   double d,dmax,eps,*a=NULL;
   
   if (init==0)
      init_rp();

   check_global_int("set_rat_parms",4,irp,num,den,degree);
   check_global_dble("set_rat_parms",3,range[0],range[1],mu0);
   if((num!=-1)||(den!=2))
   {
      check_global_dble("set_rat_parms",2,A,delta);
      check_global_dblearray("set_rat_parms",degree,nu);
      check_global_dblearray("set_rat_parms",degree,mu);
      check_global_dblearray("set_rat_parms",2*(degree+1),x);
   }

   ie=0;
   ie|=((num==0)||(den<=0));
   ie|=((irp<0)||(irp>=IRPMAX));
   ie|=(degree<1);
   ie|=(range[0]>=range[1]);
   ie|=(range[0]<=0.0);

   error_root(ie!=0,1,"set_rat_parms [rat_parms.c]",
              "Parameters are out of range");
   
   error_root(rp[irp].degree!=0,1,"set_rat_parms [rat_parms.c]",
              "Attempt to reset an already specified parameter set");

   rp[irp].nu=malloc(2*degree*sizeof(double));
   error_root(rp[irp].nu==NULL,1,"set_rat_parms [rat_parms.c]",
              "Unable to allocate memory");
   rp[irp].mu=rp[irp].nu+degree;


   rp[irp].power[0]=num;
   rp[irp].power[1]=den;
   rp[irp].degree=degree;
   rp[irp].range[0]=range[0];
   rp[irp].range[1]=range[1];
   rp[irp].mu0=mu0;
   
   if((num==-1)&&(den==2))
   {
      a=malloc(2*degree*sizeof(double));
      error_root(a==NULL,1,"set_rat_parms [rat_parms.c]",
                 "Unable to allocate memory");

      eps=(range[0]*range[0]+mu0*mu0)/(range[1]*range[1]+mu0*mu0);
      
      zolotarev(degree,eps,&A,a,&delta);

      rp[irp].delta=delta;
      rp[irp].A=A/range[1];
      for (i=0;i<degree;i++)
      {
         rp[irp].mu[i]=sqrt(range[1]*range[1]*a[2*i+1]+mu0*mu0);
         rp[irp].nu[i]=sqrt(range[1]*range[1]*a[2*i]+mu0*mu0);
      }

      free(a);
   }
   else
   {
      rp[irp].delta=fabs(delta);
      rp[irp].A=A;
      for (i=0;i<degree;i++)
      {
         rp[irp].mu[i]=mu[i];
         rp[irp].nu[i]=nu[i];
      }
      
      ie=0;
      dmax=0.0;
      for (i=0;i<2*(degree+1);i++)
      {
         ie|=(x[i]<range[0]*range[0]);
         ie|=(x[i]>range[1]*range[1]);
         d=fabs(1.0-(1-2*(i%2))*rat_delta(degree,A,nu,mu,(1.0*num)/den,x[i],mu0)/delta);
         if(d>dmax) dmax=d;
      }
      
      error_root(ie!=0,1,"set_rat_parms [rat_parms.c]",
         "The given extrema are outside of the expected range");
      error_root(dmax>1.e-2,1,"set_rat_parms [rat_parms.c]",
         "The given rational approximation does not reproduce the expected relative errors on extrema with enough precision");
   }

   return rp[irp];
}


rat_parms_t rat_parms(int irp)
{
   if (init==0)
      init_rp();

   if ((irp>=0)&&(irp<IRPMAX))
      return rp[irp];
   else
   {
      error_loc(1,1,"rat_parms [rat_parms.c]",
                "Rational function index is out of range");
      return rp[IRPMAX];
   }
}


void read_rat_parms(int irp)
{
   int my_rank;
   int degree=0,num,den,i;
   double range[2],delta,mu0;
   double A,*nu=NULL,*mu=NULL,*x=NULL;
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   if (my_rank==0)
   {
      sprintf(line,"Rational %d",irp);
      find_section(line);

      read_line("power","%d %d",&num,&den);
      read_line("degree","%d",&degree);
      read_line("range","%lf %lf",range,range+1);
      read_line("mu","%lf",&mu0);
   }

   MPI_Bcast(&degree,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&num,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&den,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(range,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   if ((num!=-1)||(den!=2))
   {
      nu=malloc(2*degree*sizeof(double));
      mu=nu+degree;
      error_root(nu==NULL,1,"read_rat_parms [rat_parms.c]",
                 "Unable to allocate memory");

      x=malloc(2*(degree+1)*sizeof(double));
      error_root(x==NULL,1,"read_rat_parms [rat_parms.c]",
                 "Unable to allocate memory");

      if (my_rank==0)
      {
         read_line("delta","%lf",&delta);
         read_line("A","%lf",&A);
         for (i=0;i<degree;i++)
         {
            sprintf(line,"nu[%d]",i);
            read_line(line,"%lf",nu+i);
            sprintf(line,"mu[%d]",i);
            read_line(line,"%lf",mu+i);
         }
         for (i=0;i<2*(degree+1);i++)
         {
            sprintf(line,"x[%d]",i);
            read_line(line,"%lf",x+i);
         }
      }
   }

   MPI_Bcast(&delta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&A,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if ((num!=-1)||(den!=2))
   {
      MPI_Bcast(nu,degree,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(mu,degree,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(x,2*(degree+1),MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   
   set_rat_parms(irp,num,den,degree,range,mu0,A,nu,mu,delta,x);
   
   if ((num!=-1)||(den!=2))
   {
      free(nu);
      free(x);
   }
}


void print_rat_parms(void)
{
   int my_rank,irp,n[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   
   if ((my_rank==0)&&(init==1))
   {
      for (irp=0;irp<IRPMAX;irp++)
      {
         if (rp[irp].degree!=0)
         {
            printf("Rational %d:\n",irp);
            printf("power = %d/%d\n",rp[irp].power[0],rp[irp].power[1]);
            printf("degree = %d\n",rp[irp].degree);
            n[0]=fdigits(rp[irp].range[0]);
            n[1]=fdigits(rp[irp].range[1]);            
            printf("range = [%.*f,%.*f]\n",IMAX(n[0],1),rp[irp].range[0],
                   IMAX(n[1],1),rp[irp].range[1]);
            n[0]=fdigits(rp[irp].mu0);
            printf("mu = %.*f\n",IMAX(n[0],1),rp[irp].mu0);
            printf("delta = %.2e\n\n",rp[irp].delta);
         }
      }
   }
}


void write_rat_parms(FILE *fdat)
{
   int irp;
   
   if (init==1)
   {
      for (irp=0;irp<IRPMAX;irp++)
      {
         if (rp[irp].degree!=0)
         {
            write_little_int(1,fdat,4,irp,rp[irp].power[0],rp[irp].power[1],
                                    rp[irp].degree);
            write_little_dble(1,fdat,4,rp[irp].range[0],rp[irp].range[1],
                                     rp[irp].A,rp[irp].delta);
            write_little_dblearray(1,fdat,rp[irp].degree,rp[irp].nu);
            write_little_dblearray(1,fdat,rp[irp].degree,rp[irp].mu);
         }
      }
   }
}


void check_rat_parms(FILE *fdat)
{
   int irp;
      
   if (init==1)
   {
      for (irp=0;irp<IRPMAX;irp++)
      {
         if (rp[irp].degree!=0)
         {
            check_little_int("check_rat_parms",fdat,4,
                           irp,rp[irp].power[0],rp[irp].power[1],
                           rp[irp].degree);
            check_little_dble("check_rat_parms",fdat,4,
                            rp[irp].range[0],rp[irp].range[1],
                            rp[irp].A,rp[irp].delta);
            check_little_dblearray("check_rat_parms",fdat,
                                 rp[irp].degree,rp[irp].nu);
            check_little_dblearray("check_rat_parms",fdat,
                                 rp[irp].degree,rp[irp].mu);
         }
      }
   }
}
