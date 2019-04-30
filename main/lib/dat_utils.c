
/*******************************************************************************
*
* File dat_utils.c
*
* Copyright (C) 2017 Agostino Patella
*
* Based on openQCD-1.6/main/qcd1.c
*      and openQCD-1.6/main/ym1.c
* Copyright (C) 2011-2013, 2016 Martin Luescher, Isabel Campos
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#include <stdio.h>
#include "flags.h"
#include "utils.h"
#include "global.h"
#include "main_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)


static struct
{
   int nt,iac;
   double dH,avpl3,avpl1;
} dat;


void write_dat(FILE *fdat)
{
   int ew,iw;

   ew=3;
   if ((gauge()&1)!=0) ew++;
   if ((gauge()&2)!=0) ew++;

   iw=write_little_int(0,fdat,2,dat.nt,dat.iac);
   iw+=write_little_dble(0,fdat,1,dat.dH);
   if ((gauge()&1)!=0)
      iw+=write_little_dble(0,fdat,1,dat.avpl3);
   if ((gauge()&2)!=0)
      iw+=write_little_dble(0,fdat,1,dat.avpl1);

   error_root(iw!=ew,1,"write_dat [dat_utils.c]",
              "Incorrect write count");
}


int read_dat(FILE *fdat,int *nt)
{
   int er,ir;
   
   er=3;
   if ((gauge()&1)!=0) er++;
   if ((gauge()&2)!=0) er++;

   ir=read_little_int(0,fdat,1,&(dat.nt));
   if(ir!=1) return 0;
   
   ir+=read_little_int(0,fdat,1,&(dat.iac));
   
   ir+=read_little_dble(0,fdat,1,&(dat.dH));
   
   if ((gauge()&1)!=0)
      ir+=read_little_dble(0,fdat,1,&(dat.avpl3));
   
   if ((gauge()&2)!=0)
      ir+=read_little_dble(0,fdat,1,&(dat.avpl1));

   error_root(ir!=er,1,"read_dat [dat_utils.c]",
              "Incorrect read count");

   (*nt)=dat.nt;
   
   return 1;
}


void set_dat(int nt,int iac,double dH,double avpl3,double avpl1)
{
   dat.nt=nt;
   dat.iac=iac;
   dat.dH=dH;
   dat.avpl3=avpl3;
   dat.avpl1=avpl1;
}


void print_dat(void)
{
   printf("Trajectory no %d\n",dat.nt);
   printf("dH = %+.1e, ",dat.dH);
   printf("iac = %d\n",dat.iac);
   if ((gauge()&1)!=0)
      printf("Average SU(3) plaquette = %.6f\n",dat.avpl3);
   if ((gauge()&2)!=0)
      printf("Average U(1) plaquette = %.6f\n",dat.avpl1);
}
