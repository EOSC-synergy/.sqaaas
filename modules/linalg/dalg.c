
/*******************************************************************************
*
* File dalg.c
*
* Copyright (C) 2016 Agostino Patella
* Copyright (C) 2015 Marina Marinkovic
* Copyright (C) 2005, 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic functions for double-precision fields
*
* The externally accessible functions are
*
*   void random_dvec(int vol,double *X)
*     Initializes the elements X to random values
*     with distribution proportional to exp{1/2 X^2}.
*
*   double norm_square_dvec(int vol,int icom,double *X)
*     Computes the square of the norm of the field X.
*
*   double scalar_prod_dvec(int vol,int icom,double *X,double *Y)
*     Computes the scalar product of the fields X and Y.
*
*   void set_dvec2zero(int vol,double *X)
*     Sets the elements of X to zero.
*
*   void assign_dvec2dvec(int vol,double *X,double *Y)
*     Assigns the field X to the field Y.
*
*   void swap_dvec(int vol,double *X,double *Y)
*     Swaps the fields X and Y.
*
*   void muladd_assign_dvec(int vol,double r,double *X,double *Y)
*     Adds r*X to Y.
*
* All programs in this module operate on arrays of doubles whose
* base address is passed through the arguments. The length of the array is
* specified by the parameter vol. Scalar products etc. are globally summed if
* the parameter icom is equal to 1. In this case the calculated values are
* guaranteed to be exactly the same on all processes.
*
*******************************************************************************/

#define DALG_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "random.h"
#include "linalg.h"
#include "global.h"

#define MAX_LEVELS 12
#define BLK_LENGTH 8

static int cnt[MAX_LEVELS];
static double smx[MAX_LEVELS];
static double c1=0.0,rb,sm;


void random_dvec(int vol,double *X)
{
   double *Xm;
   if (c1==0.0)
   {
      c1=sqrt(2.0);
   }
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      gauss_dble(&rb,1);
      (*X)=c1*rb;
   }
}


double norm_square_dvec(int vol,int icom,double *X)
{
   int n;
   double *Xm;
   
   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }
   
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      sm=((*X)*(*X));
      cnt[0]+=1;
      smx[0]+=sm;
      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];
         cnt[n-1]=0;
         smx[n-1]=0.0;
      }
   }
   
   for (n=1;n<MAX_LEVELS;n++)
      smx[0]+=smx[n];
   
   if ((icom==1)&&(NPROC>1))
   {
      sm=smx[0];
      MPI_Reduce(&sm,smx,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(smx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   
   return smx[0];
}


double scalar_prod_dvec(int vol,int icom,double *X,double *Y)
{
   int n;
   double *Xm;
   
   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }
   
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      sm=((*X)*(*Y));
      
      Y+=1;
      cnt[0]+=1;
      smx[0]+=sm;
      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];
         cnt[n-1]=0;
         smx[n-1]=0.0;
      }
   }
   
   for (n=1;n<MAX_LEVELS;n++)
      smx[0]+=smx[n];
      
   if ((icom==1)&&(NPROC>1))
   {
      sm=smx[0];
      MPI_Reduce(&sm,smx,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(smx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   
   return smx[0];
}


void set_dvec2zero(int vol,double *X)
{
   double *Xm;
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      (*X)=0.0;
   }
}


void assign_dvec2dvec(int vol,double *X,double *Y)
{
   double *Xm;
   
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      (*Y)=(*X);
      Y+=1;
   }
}


void swap_dvec(int vol,double *X,double *Y)
{
   double r;
   double *Xm;
   
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      r=(*Y);
      (*Y)=(*X);
      (*X)=r;
      
      Y+=1;
   }
}

void muladd_assign_dvec(int vol,double r,double *X,double *Y)
{
   double *Xm;
   
   Xm=X+vol;
   
   for (;X<Xm;X++)
   {
      (*Y)+=r*(*X);
      Y+=1;
   }
}
