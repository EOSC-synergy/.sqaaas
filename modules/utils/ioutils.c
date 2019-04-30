/*******************************************************************************
*
* File ioutils.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on various code in openQCD-1.6
* Copyright (C) 2011-2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Basic I/O utility functions.
*
* The externally accessible functions are
*
*   int write_little_dble(int blocking,FILE *fdat,int nargs,...)
*     Attempts to write the variables passed as [...] in the fdat file in
*     little endian. A number of nargs double arguments is expected instead
*     of [...].
*     If iw==narg, where iw is the number of double values successfully written,
*     the program returns iw and the value of blocking is irrelevant.
*     If iw!=narg, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns iw, otherwise it
*     terminates with an error.
*
*   int write_little_int(int blocking,FILE *fdat,int nargs,...)
*     Attempts to write the variables passed as [...] in the fdat file in
*     little endian. A number of nargs int arguments is expected instead
*     of [...].
*     If iw==narg, where iw is the number of int values successfully written,
*     the program returns iw and the value of blocking is irrelevant.
*     If iw!=narg, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns iw, otherwise it
*     terminates with an error.
*
*   int write_little_dblearray(int blocking,FILE *fdat,int size,double *data)
*     Attempts to write to write the array data[size] in the fdat file in
*     little endian.
*     If iw==size, where iw is the number of double values successfully written,
*     the program returns iw and the value of blocking is irrelevant.
*     If iw!=size, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns iw, otherwise it
*     terminates with an error.
*
*   int write_little_intarray(int blocking,FILE *fdat,int size,int *data)
*     Attempts to write to write the array data[size] in the fdat file in
*     little endian.
*     If iw==size, where iw is the number of int values successfully written,
*     the program returns iw and the value of blocking is irrelevant.
*     If iw!=size, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns iw, otherwise it
*     terminates with an error.
*
*   int read_little_dble(int blocking,FILE *fdat,int nargs,...)
*     Attempts to read the next nargs double values from the file fdat, and
*     stores them in the variables pointed by the pointers passed as [...].
*     A number of nargs double* arguments is expected instead of [...].
*     If ir==narg, where ir is the number of double values successfully read,
*     the program returns ir and the value of blocking is irrelevant.
*     If ir!=narg, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns ir, otherwise it
*     terminates with an error.
*
*   int read_little_int(int blocking,FILE *fdat,int nargs,...)
*     Attempts to read the next nargs int values from the file fdat, and
*     stores them in the variables pointed by the pointers passed as [...].
*     A number of nargs int* arguments is expected instead of [...].
*     If ir==narg, where ir is the number of int values successfully read,
*     the program returns ir and the value of blocking is irrelevant.
*     If ir!=narg, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns ir, otherwise it
*     terminates with an error.
*
*   int read_little_dblearray(int blocking,FILE *fdat,int size,double *data)
*     Attempts to read the next size double values from the file fdat, and
*     stores them in the array data[size].
*     If ir==size, where ir is the number of double values successfully read,
*     the program returns ir and the value of blocking is irrelevant.
*     If ir!=size, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns ir, otherwise it
*     terminates with an error.
*
*   int read_little_intarray(int blocking,FILE *fdat,int size,int *data)
*     Attempts to read the next size int values from the file fdat, and
*     stores them in the array data[size].
*     If ir==size, where ir is the number of int values successfully read,
*     the program returns ir and the value of blocking is irrelevant.
*     If ir!=size, the behaviour of the program depends on the value of
*     blocking. If blocking==0 then the program returns ir, otherwise it
*     terminates with an error.
*
*   void check_little_dble(const char fnm[256],FILE *fdat,int nargs,...)
*     Reads the next nargs double values from the file fdat and checks that
*     they match the values of the variables passed as [...]. A number of nargs
*     double arguments is expected instead of [...]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*   void check_little_int(const char fnm[256],FILE *fdat,int nargs,...)
*     Reads the next nargs int values from the file fdat and checks that
*     they match the values of the variables passed as [...]. A number of nargs
*     double arguments is expected instead of [...]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*   void check_little_dblearray(const char fnm[256],FILE *fdat,int size,
*                             double *data)
*     Reads the next size double values from the file fdat and checks that
*     they match the values of the array data[size]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*   void check_little_intarray(const char fnm[256],FILE *fdat,int size,int *data)
*     Reads the next size int values from the file fdat and checks that
*     they match the values of the array data[size]. The string fnm is expected
*     to contain the name of the function that calls the test (used in the error
*     message).
*
*******************************************************************************/

#define IOUTILS_C

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include "mpi.h"
#include "flags.h"
#include "global.h"


char msg[1024];


int write_little_dble(int blocking,FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int iw,i;
   double *dprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   iw=0;

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

      error_root(blocking&&(iw!=nargs),1,"write_little_dble [ioutils.c]",
                 "Incorrect write count");
      
      free(dprms);
   }
   
   return iw;
}


int write_little_int(int blocking,FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int iw,i;
   stdint_t *iprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   iw=0;

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

      error_root(blocking&&(iw!=nargs),1,"write_little_int [ioutils.c]",
                 "Incorrect write count");
      
      free(iprms);
   }
   
   return iw;
}


int write_little_dblearray(int blocking,FILE *fdat,int size,double *data)
{
   int my_rank,endian;
   int iw,i;
   double *dprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   iw=0;

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

      error_root(blocking&&(iw!=size),1,"write_little_dblearray [ioutils.c]",
                 "Incorrect write count");
      
      free(dprms);
   }
   
   return iw;
}


int write_little_intarray(int blocking,FILE *fdat,int size,int *data)
{
   int my_rank,endian;
   int iw,i;
   stdint_t *iprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   iw=0;

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

      error_root(blocking&&(iw!=size),1,"write_little_intarray [ioutils.c]",
                 "Incorrect write count");
      
      free(iprms);
   }
   
   return iw;
}


int read_little_dble(int blocking,FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int ir,i;
   double *dprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   ir=0;

   if (my_rank==0)
   {
      dprms=malloc(nargs*sizeof(double));

      ir=fread(dprms,sizeof(double),nargs,fdat);

      error_root(blocking&&(ir!=nargs),1,"read_little_dble [ioutils.c]",
                 "Incorrect read count");
      
      if(ir==nargs)
      {
         if (endian==BIG_ENDIAN)
            bswap_double(nargs,dprms);

         va_start(args,nargs);
         for(i=0;i<nargs;i++)
         {
            *(va_arg(args,double*))=dprms[i];
         }
         va_end(args);
      }

      free(dprms);
   }
   
   return ir;
}


int read_little_int(int blocking,FILE *fdat,int nargs,...)
{
   int my_rank,endian;
   int ir,i;
   stdint_t *iprms;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   ir=0;

   if (my_rank==0)
   {
      iprms=malloc(nargs*sizeof(stdint_t));

      ir=fread(iprms,sizeof(stdint_t),nargs,fdat);

      error_root(blocking&&(ir!=nargs),1,"read_little_int [ioutils.c]",
                 "Incorrect read count");
      
      if(ir==nargs)
      {
         if (endian==BIG_ENDIAN)
            bswap_int(nargs,iprms);

         va_start(args,nargs);
         for(i=0;i<nargs;i++)
         {
            *(va_arg(args,int*))=iprms[i];
         }
         va_end(args);
      }

      free(iprms);
   }
   
   return ir;
}


int read_little_dblearray(int blocking,FILE *fdat,int size,double *data)
{
   int my_rank,endian;
   int ir,i;
   double *dprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   ir=0;

   if (my_rank==0)
   {
      dprms=malloc(size*sizeof(double));

      ir=fread(dprms,sizeof(double),size,fdat);

      error_root(blocking&&(ir!=size),1,"read_little_dblearray [ioutils.c]",
                 "Incorrect read count");
      
      if(ir==size)
      {
         if (endian==BIG_ENDIAN)
            bswap_double(size,dprms);

         for(i=0;i<size;i++)
         {
            data[i]=dprms[i];
         }
      }

      free(dprms);
   }
   
   return ir;
}


int read_little_intarray(int blocking,FILE *fdat,int size,int *data)
{
   int my_rank,endian;
   int ir,i;
   stdint_t *iprms;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   ir=0;
   
   if (my_rank==0)
   {
      iprms=malloc(size*sizeof(stdint_t));

      ir=fread(iprms,sizeof(stdint_t),size,fdat);

      error_root(blocking&&(ir!=size),1,"read_little_intarray [ioutils.c]",
                 "Incorrect read count");
      
      if(ir==size)
      {
         if (endian==BIG_ENDIAN)
            bswap_int(size,iprms);

         for(i=0;i<size;i++)
         {
            data[i]=(int)(iprms[i]);
         }
      }

      free(iprms);
   }
   
   return ir;
}


void check_little_dble(const char fnm[256],FILE *fdat,int nargs,...)
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
      error_root(ir!=nargs,1,"check_little_dble [ioutils.c]",msg);

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
      error_root(ie!=0,1,"check_little_dble [ioutils.c]",msg);

      free(dprms);
   }
}


void check_little_int(const char fnm[256],FILE *fdat,int nargs,...)
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
      error_root(ir!=nargs,1,"check_little_int [ioutils.c]",msg);

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
      error_root(ie!=0,1,"check_little_int [ioutils.c]",msg);

      free(iprms);
   }
}


void check_little_dblearray(const char fnm[256],FILE *fdat,int size,double *data)
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
      error_root(ir!=size,1,"check_little_dblearray [ioutils.c]",msg);

      if (endian==BIG_ENDIAN)
         bswap_double(size,dprms);

      ie=0;
      for(i=0;i<size;i++)
      {
         ie|=(dprms[i]!=data[i]);
      }

      sprintf(msg,"Parameters do not match -- Called by %s",fnm);
      error_root(ie!=0,1,"check_little_dblearray [ioutils.c]",msg);

      free(dprms);
   }
}


void check_little_intarray(const char fnm[256],FILE *fdat,int size,int *data)
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
      error_root(ir!=size,1,"check_little_intarray [ioutils.c]",msg);

      if (endian==BIG_ENDIAN)
         bswap_int(size,iprms);

      ie=0;
      for(i=0;i<size;i++)
      {
         ie|=(iprms[i]!=(stdint_t)(data[i]));
      }

      sprintf(msg,"Parameters do not match -- Called by %s",fnm);
      error_root(ie!=0,1,"check_little_intarray [ioutils.c]",msg);

      free(iprms);
   }
}
