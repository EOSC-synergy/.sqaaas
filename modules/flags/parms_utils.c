/*******************************************************************************
*
* File lat_parms.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/modules/flags/lat_parms.c
* Copyright (C) 2009-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Utilities for parameter management.
*
* The externally accessible functions are
*
*   void check_global_dble(const char fnm[256],int nargs,...)
*     Checks if all variables passed as [...] take the same value in all MPI
*     processes. A number of nargs double arguments is expected instead of
*     [...]. The string fnm is expected to contain the name of the function
*     that calls the test (used in the error message).
*
*   void check_global_dblearray(const char fnm[256],int size,double *data)
*     Checks if the array data[size] contains the same values in all MPI
*     processes. The string fnm is expected to contain the name of the function
*     that calls the test (used in the error message).
*
*   void check_global_int(const char fnm[256],int nargs,...)
*     Checks if all variables passed as [...] take the same value in all MPI
*     processes. A number of nargs int arguments is expected instead of [...].
*     The string fnm is expected to contain the name of the function that calls
*     the test (used in the error message).
*
*   void check_global_intarray(const char fnm[256],int size,int *data)
*     Checks if the array data[size] contains the same values in all MPI
*     processes. The string fnm is expected to contain the name of the function
*     that calls the test (used in the error message).
*
*   void write_little_dble(FILE *fdat,int nargs,...)
*     Writes the variables passed as [...] in the fdat file in little endian.
*     A number of nargs double arguments is expected instead of [...].
*
*   void write_little_int(FILE *fdat,int nargs,...)
*     Writes the variables passed as [...] in the fdat file in little endian.
*     A number of nargs int arguments is expected instead of [...].
*
*   void write_little_dblearray(FILE *fdat,int size,double *data)
*     Writes the array data[size] in the fdat file in little endian.
*
*   void write_little_intarray(FILE *fdat,int size,int *data)
*     Writes the array data[size] in the fdat file in little endian.
*
*   void check_fpar_dble(const char fnm[256],FILE *fdat,int nargs,...)
*     Reads the next nargs double values from the file fdat and checks that
*     they match the values of the variables passed as [...]. A number of nargs
*     double arguments is expected instead of [...]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*   void check_fpar_int(const char fnm[256],FILE *fdat,int nargs,...)
*     Reads the next nargs int values from the file fdat and checks that
*     they match the values of the variables passed as [...]. A number of nargs
*     double arguments is expected instead of [...]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*   void check_fpar_dblearray(const char fnm[256],FILE *fdat,int size,
*                             double *data)
*     Reads the next size double values from the file fdat and checks that
*     they match the values of the array data[size]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*   void check_fpar_intarray(const char fnm[256],FILE *fdat,int size,int *data)
*     Reads the next size int values from the file fdat and checks that
*     they match the values of the array data[size]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*******************************************************************************/

#define PARMS_UTILS_C

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "mpi.h"
#include "flags.h"
#include "global.h"


char msg[1024];


void check_global_dble(const char fnm[256],int nargs,...)
{
   double *dprms;
   int i,ie;
   va_list args;
   
   if (NPROC>1)
   {
      dprms=malloc(nargs*sizeof(double));
      
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         dprms[i]=va_arg(args,double);
      }
      va_end(args);
      
      MPI_Bcast(dprms,nargs,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
      ie=0;
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         ie|=(dprms[i]!=va_arg(args,double));
      }
      va_end(args);
      
      sprintf(msg,"Parameters are not global -- Called by %s",fnm);
      error(ie!=0,1,"check_global_dble [parms_utils.c]",msg);

      free(dprms);
   }
}


void check_global_dblearray(const char fnm[256],int size,double *data)
{
   double *dprms;
   int i,ie;
   
   if (NPROC>1)
   {
      dprms=malloc(size*sizeof(double));
      
      for(i=0;i<size;i++)
      {
         dprms[i]=data[i];
      }
      
      MPI_Bcast(dprms,size,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
      ie=0;
      for(i=0;i<size;i++)
      {
         ie|=(dprms[i]!=data[i]);
      }
      
      sprintf(msg,"Parameters are not global -- Called by %s",fnm);
      error(ie!=0,1,"check_global_dblearray [parms_utils.c]",msg);
      
      free(dprms);
   }
}


void check_global_int(const char fnm[256],int nargs,...)
{
   int *iprms;
   int i,ie;
   va_list args;
   
   if (NPROC>1)
   {
      iprms=malloc(nargs*sizeof(int));
      
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         iprms[i]=va_arg(args,int);
      }
      va_end(args);
      
      MPI_Bcast(iprms,nargs,MPI_INT,0,MPI_COMM_WORLD);
   
      ie=0;
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         ie|=(iprms[i]!=va_arg(args,int));
      }
      va_end(args);
      
      sprintf(msg,"Parameters are not global -- Called by %s",fnm);
      error(ie!=0,1,"check_global_int [parms_utils.c]",msg);
      
      free(iprms);
   }
}


void check_global_intarray(const char fnm[256],int size,int *data)
{
   int *iprms;
   int i,ie;
   
   if (NPROC>1)
   {
      iprms=malloc(size*sizeof(int));
      
      for(i=0;i<size;i++)
      {
         iprms[i]=data[i];
      }
      
      MPI_Bcast(iprms,size,MPI_INT,0,MPI_COMM_WORLD);
   
      ie=0;
      for(i=0;i<size;i++)
      {
         ie|=(iprms[i]!=data[i]);
      }
      
      sprintf(msg,"Parameters are not global -- Called by %s",fnm);
      error(ie!=0,1,"check_global_intarray [parms_utils.c]",msg);
   
      free(iprms);
   }
}


void write_little_dble(FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int iw,i;
   double *dprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      dprms=malloc(nargs*sizeof(double));
      
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         dprms[i]=va_arg(args,double);
      }
      va_end(args);

      if (endian==BIG_ENDIAN)
         bswap_double(nargs,dprms);

      iw=fwrite(dprms,sizeof(double),nargs,fdat);

      error_root(iw!=nargs,1,"write_little_dble [parms_utils.c]",
                 "Incorrect write count");
      
      free(dprms);
   }
}


void write_little_int(FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int iw,i;
   stdint_t *iprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      iprms=malloc(nargs*sizeof(stdint_t));
      
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         iprms[i]=(stdint_t)(va_arg(args,int));
      }
      va_end(args);

      if (endian==BIG_ENDIAN)
         bswap_int(nargs,iprms);

      iw=fwrite(iprms,sizeof(stdint_t),nargs,fdat);

      error_root(iw!=nargs,1,"write_little_int [parms_utils.c]",
                 "Incorrect write count");
      
      free(iprms);
   }
}


void write_little_dblearray(FILE *fdat,int size,double *data)
{
   int my_rank,endian;
   int iw,i;
   double *dprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      dprms=malloc(size*sizeof(double));
      
      for(i=0;i<size;i++)
      {
         dprms[i]=data[i];
      }

      if (endian==BIG_ENDIAN)
         bswap_double(size,dprms);

      iw=fwrite(dprms,sizeof(double),size,fdat);

      error_root(iw!=size,1,"write_little_dblearray [parms_utils.c]",
                 "Incorrect write count");
      
      free(dprms);
   }
}


void write_little_intarray(FILE *fdat,int size,int *data)
{
   int my_rank,endian;
   int iw,i;
   stdint_t *iprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      iprms=malloc(size*sizeof(stdint_t));
      
      for(i=0;i<size;i++)
      {
         iprms[i]=(stdint_t)(data[i]);
      }

      if (endian==BIG_ENDIAN)
         bswap_int(size,iprms);

      iw=fwrite(iprms,sizeof(stdint_t),size,fdat);

      error_root(iw!=size,1,"write_little_intarray [parms_utils.c]",
                 "Incorrect write count");
      
      free(iprms);
   }
}


void check_fpar_dble(const char fnm[256],FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int ir,ie,i;
   double *dprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      dprms=malloc(nargs*sizeof(double));

      ir=fread(dprms,sizeof(double),nargs,fdat);

      sprintf(msg,"Incorrect read count -- Called by %s",fnm);
      error_root(ir!=nargs,1,"check_fpar_dble [parms_utils.c]",msg);

      if (endian==BIG_ENDIAN)
         bswap_double(nargs,dprms);

      ie=0;
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         ie|=(dprms[i]!=(va_arg(args,double)));
      }
      va_end(args);

      sprintf(msg,"Parameters do not match -- Called by %s",fnm);
      error_root(ie!=0,1,"check_fpar_dble [parms_utils.c]",msg);

      free(dprms);
   }
}


void check_fpar_int(const char fnm[256],FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int ir,ie,i;
   stdint_t *iprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      iprms=malloc(nargs*sizeof(stdint_t));

      ir=fread(iprms,sizeof(stdint_t),nargs,fdat);

      sprintf(msg,"Incorrect read count -- Called by %s",fnm);
      error_root(ir!=nargs,1,"check_fpar_int [parms_utils.c]",msg);

      if (endian==BIG_ENDIAN)
         bswap_int(nargs,iprms);

      ie=0;
      va_start(args,nargs);
      for(i=0;i<nargs;i++)
      {
         ie|=(iprms[i]!=(stdint_t)(va_arg(args,int)));
      }
      va_end(args);

      sprintf(msg,"Parameters do not match -- Called by %s",fnm);
      error_root(ie!=0,1,"check_fpar_int [parms_utils.c]",msg);

      free(iprms);
   }
}


void check_fpar_dblearray(const char fnm[256],FILE *fdat,int size,double *data)
{
   int my_rank,endian;
   int ir,ie,i;
   double *dprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      dprms=malloc(size*sizeof(double));

      ir=fread(dprms,sizeof(double),size,fdat);

      sprintf(msg,"Incorrect read count -- Called by %s",fnm);
      error_root(ir!=size,1,"check_fpar_dblearray [parms_utils.c]",msg);

      if (endian==BIG_ENDIAN)
         bswap_double(size,dprms);

      ie=0;
      for(i=0;i<size;i++)
      {
         ie|=(dprms[i]!=data[i]);
      }

      sprintf(msg,"Parameters do not match -- Called by %s",fnm);
      error_root(ie!=0,1,"check_fpar_dblearray [parms_utils.c]",msg);

      free(dprms);
   }
}


void check_fpar_intarray(const char fnm[256],FILE *fdat,int size,int *data)
{
   int my_rank,endian;
   int ir,ie,i;
   stdint_t *iprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();

   if (my_rank==0)
   {
      iprms=malloc(size*sizeof(stdint_t));

      ir=fread(iprms,sizeof(stdint_t),size,fdat);

      sprintf(msg,"Incorrect read count -- Called by %s",fnm);
      error_root(ir!=size,1,"check_fpar_intarray [parms_utils.c]",msg);

      if (endian==BIG_ENDIAN)
         bswap_int(size,iprms);

      ie=0;
      for(i=0;i<size;i++)
      {
         ie|=(iprms[i]!=(stdint_t)(data[i]));
      }

      sprintf(msg,"Parameters do not match -- Called by %s",fnm);
      error_root(ie!=0,1,"check_fpar_intarray [parms_utils.c]",msg);

      free(iprms);
   }
}
