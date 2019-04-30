
/*******************************************************************************
*
* File utils.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* based on the openQCD-1.6 version of this file written by
* Copyright (C) 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)

* Collection of basic utility programs
*
* The externally accessible functions are
*
*   int find_opt(int argc,char *argv[],char *opt)
*     This program compares the string opt with the arguments
*     argv[1],..,argv[argc-1] and returns the position of the first argument
*     that matches the string. If there is no matching argument, or if the
*     program is called from another process, the return value is 0.
*
*   void error(int test,int no,char *name,char *format,...)
*     Checks whether "test"=0 and if not aborts the program gracefully
*     with error number "no" after printing the "name" of the calling
*     program and an error message to stdout. The message is formed using
*     the "format" string and any additional arguments, exactly as in a
*     printf statement
*
*   void message(char *format,...)
*     Same as printf(), provided for compatibility
*
*******************************************************************************/

#define UTILS_C

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "utils.h"


int find_opt(int argc,char *argv[],char *opt)
{
   int k;

   for (k=1;k<argc;k++)
   {
      if (strcmp(argv[k],opt)==0)
         return k;
   }
   
   return 0;
}


void error(int test,int no,char *name,char *format,...)
{
   va_list args;

   if (test!=0)
   {
      printf("\nError in %s:\n",name);
      va_start(args,format);
      vprintf(format,args);
      va_end(args);
      printf("\nProgram aborted\n\n");
      exit(no);
   }
}


void message(char *format,...)
{
   va_list args;

   va_start(args,format);
   vprintf(format,args);
   va_end(args);
}
