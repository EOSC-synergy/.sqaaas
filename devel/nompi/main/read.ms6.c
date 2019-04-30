
/*******************************************************************************
*
* File read.ms6.c
*
* Copyright (C) 2017 Agostino Patella
* 
* Based on openQCD-1.6/devel/nompi/main/read1.c
* Copyright (C) 2010-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reads and evaluates data from the *.ms6.dat files created by the programs ms6,
* qcd1 and ym1. The file to be read has to be specified on the command line.
*
* This program writes the history of the SU(3) Wilson flow observables to the
* file <run name>.plot.ms6.<obs>.dat in the plots directory, and their
* expectation values with error and integrated autocorrelation time to the file
* <run name>.stat.ms6.<obs>.dat in the stats directory.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include "extras.h"
#include "random.h"

static struct
{
   int tmax,x0,nsrc;
   int nfl,qhat[2];
   double kappa[2],mu[2];
   double su3csw[2],u1csw[2];
   double cF[2],cF_prime[2];
   double th1[2],th2[2],th3[2];
} file_head;

typedef struct
{
   int nt;
   double *P,*A0,*A1,*A2,*A3;
} dat_t;

static int range,all,noauto;
static int nms,nfirst,nlast,neff;
static dat_t *adat;


static void read_file(char *fin)
{
   long ipos;
   FILE *fdat;
   int ir,id,t,tmax,endian,ifl;
   stdint_t istd[4];
   double dstd[9];
   double *p;
   static double *dbuf=NULL;


   fdat=fopen(fin,"rb");
   error(fdat==NULL,1,"read_file [read.ms6.c]","Unable to open data file");

   printf("Read data from file %s\n\n",fin);
   

   /* Read header */

   endian=endianness();

   ir=fread(istd,sizeof(stdint_t),4,fdat);
   error_root(ir!=4,1,"read_file [read.ms6.c]",
              "Incorrect file header");

   if (endian==BIG_ENDIAN)
   {
      bswap_int(4,istd);
   }

   file_head.tmax=(int)istd[0];
   file_head.x0=(int)istd[1];
   file_head.nsrc=(int)istd[2];
   file_head.nfl=(int)istd[3];

   error_root((file_head.nfl!=1)&&(file_head.nfl!=2),1,
              "read_file [read.ms6.c]","Incorrect value for nfl");
   error_root(file_head.nsrc<=0,1,
              "read_file [read.ms6.c]","Incorrect value for nsrc");
   error_root(file_head.tmax<=0,1,
              "read_file [read.ms6.c]","Incorrect value for tmax");
   error_root(file_head.tmax<=file_head.x0,1,
              "read_file [read.ms6.c]","Incorrect value for x0");

   for (ifl=0;ifl<file_head.nfl;++ifl)
   {
      ir=fread(istd,sizeof(stdint_t),1,fdat);
      ir+=fread(dstd,sizeof(double),9,fdat);
      error_root(ir!=10,1,"read_file [read.ms6.c]",
                 "Incorrect file header");

      if (endian==BIG_ENDIAN)
      {
         bswap_int(1,istd);
         bswap_double(9,dstd);
      }

      file_head.qhat[ifl]=(int)istd[0];
      file_head.kappa[ifl]=dstd[0];
      file_head.mu[ifl]=dstd[1];
      file_head.su3csw[ifl]=dstd[2];
      file_head.u1csw[ifl]=dstd[3];
      file_head.cF[ifl]=dstd[4];
      file_head.cF_prime[ifl]=dstd[5];
      file_head.th1[ifl]=dstd[6];
      file_head.th2[ifl]=dstd[7];
      file_head.th3[ifl]=dstd[8];
   }


   /* Count measures */

   tmax=file_head.tmax;
   dbuf=amalloc(5*tmax*sizeof(double),4);
   error((dbuf==NULL),1,"read_file [read.ms6.c]",
         "Unable to allocate buffer");

   ipos=ftell(fdat);
   nms=0;
   
   while(1)
   {
      if (fread(istd,sizeof(stdint_t),1,fdat)!=1)
         break;
      if (fread(dbuf,sizeof(double),5*tmax,fdat)!=5*tmax)
         break;
      nms++;
   }

   error(nms==0,1,"read_file [read.ms6.c]",
         "Empty data file");


   /* Allocate data array */

   adat=amalloc(nms*sizeof(*adat),3);
   p=amalloc(5*tmax*nms*sizeof(*p),4);

   error((adat==NULL)||(p==NULL),1,"alloc_adat [read.ms6.c]",
         "Unable to allocate data arrays");
   
   for (id=0;id<nms;id++)
   {
      adat[id].P=p;
      adat[id].A0=p+tmax;
      adat[id].A1=p+2*tmax;
      adat[id].A2=p+3*tmax;
      adat[id].A3=p+4*tmax;

      p+=5*tmax;
   }

   
   /* Read data */
   
   fseek(fdat,ipos,SEEK_SET);

   for (id=0;id<nms;id++)
   {
      ir=fread(istd,sizeof(stdint_t),1,fdat);
      error(ir!=1,1,"read_file [read.ms6.c]",
            "Error while reading data file");

      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      
      adat[id].nt=(int)istd[0];
      
      ir=fread(dbuf,sizeof(double),5*tmax,fdat);
      error(ir!=5*tmax,1,"read_file [read.ms6.c]",
            "Error while reading data file");

      if (endian==BIG_ENDIAN)
         bswap_double(5*tmax,dbuf);
      
      p=dbuf;
      for (t=0;t<tmax;t++)
      {
         adat[id].P[t]=(*p);
         p++;
      }
      for (t=0;t<tmax;t++)
      {
         adat[id].A0[t]=(*p);
         p++;
      }
      for (t=0;t<tmax;t++)
      {
         adat[id].A1[t]=(*p);
         p++;
      }
      for (t=0;t<tmax;t++)
      {
         adat[id].A2[t]=(*p);
         p++;
      }
      for (t=0;t<tmax;t++)
      {
         adat[id].A3[t]=(*p);
         p++;
      }
   }

   fclose(fdat);
}


static void select_range(void)
{
   int n,no,nf,nl;
   int np,dn,ie;

   printf("There are %d measurements (trajectories no %d - %d).\n",
          nms,adat[0].nt,adat[nms-1].nt);
   printf("Range [nfirst,nlast] of trajectories to analyse: ");
   if (all!=0)
   {
      nfirst=adat[0].nt;
      nlast=adat[nms-1].nt;
      printf("%d %d\n",nfirst,nlast);
   }
   else if (range!=0)
   {
      printf("%d %d\n",nfirst,nlast);
   }
   else
      scanf("%d %d",&nfirst,&nlast);

   nf=0;
   nl=0;

   for (n=0;n<nms;n++)
   {
      no=adat[n].nt;

      if (no<nfirst)
         nf+=1;

      if (no<=nlast)
         nl+=1;
   }

   nfirst=nf;
   nlast=nl;
   neff=nlast-nfirst;

   printf("Keep %d measurements (trajectories no %d - %d).\n\n",
          neff,adat[nfirst].nt,adat[nlast-1].nt);

   error(neff<2,1,"select_range [read.ms6.c]",
         "Selected range contains less than 2 measurements");

   np=adat[nfirst].nt;
   dn=adat[nfirst+1].nt-adat[nfirst].nt;

   if (dn<=0)
      ie=1;
   else
   {
      ie=0;

      for (n=(nfirst+1);n<nlast;n++)
      {
         no=adat[n].nt;

         if ((no-np)!=dn)
         {
            ie=2;
            break;
         }

         np=no;
      }
   }

   error(ie!=0,1,"select_range [read.ms6.c]",
         "Varying trajectory number separation in selected range");
}


static void print_plots(char *fin)
{
   int n,ims,x0,t,tmax;
   char base[NAME_SIZE],plt_file[NAME_SIZE],*p;
   dat_t *ndat;
   FILE *fout;

   p=strstr(fin,".ms6.dat");
   error(p==NULL,1,"print_plots [read.ms6.c]","Unexpected data file name");
   n=p-fin;

   p=strrchr(fin,'/');
   if (p==NULL)
      p=fin;
   else
      p+=1;
   n-=(p-fin);

   error(n>=NAME_SIZE,1,"print_plots [read.ms6.c]","File name is too long");
   strncpy(base,p,n);
   base[n]='\0';

   error(name_size("plots/%s.plot.ms6.dA0.dat",base)>=NAME_SIZE,1,
         "print_plots [read.ms6.c]","File name is too long");

   x0=file_head.x0;
   tmax=file_head.tmax;

   printf("Data printed to files:\n");
   
   
   sprintf(plt_file,"plots/%s.plot.ms6.P.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms6.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable P\n");
   fprintf(fout,"# ---------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# P(t):     < P(%d) P(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          P(0) P(1) P(2)...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (t=0;t<tmax;t++)
         fprintf(fout," %.4e  ",(*ndat).P[t]);
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms6.P.dat\n",base);


   sprintf(plt_file,"plots/%s.plot.ms6.A0.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms6.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A0\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# A0(t):    < P(%d) A0(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          A0(0) A0(1) A0(2)...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (t=0;t<tmax;t++)
         fprintf(fout," %.4e  ",(*ndat).A0[t]);
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms6.A0.dat\n",base);


   sprintf(plt_file,"plots/%s.plot.ms6.A1.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms6.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A1\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# A1(t):    < P(%d) A1(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          A1(0) A1(1) A1(2)...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (t=0;t<tmax;t++)
         fprintf(fout," %.4e  ",(*ndat).A1[t]);
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms6.A1.dat\n",base);


   sprintf(plt_file,"plots/%s.plot.ms6.A2.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms6.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A2\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# A2(t):    < P(%d) A2(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          A2(0) A2(1) A2(2)...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (t=0;t<tmax;t++)
         fprintf(fout," %.4e  ",(*ndat).A2[t]);
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms6.A2.dat\n",base);


   sprintf(plt_file,"plots/%s.plot.ms6.A3.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms6.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A3\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# A3(t):    < P(%d) A3(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          A3(0) A3(1) A3(2)...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (t=0;t<tmax;t++)
         fprintf(fout," %.4e  ",(*ndat).A3[t]);
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms6.A3.dat\n",base);


   printf("\n");
}


static double acorr(int n,double *a,int *lambda,int *w,double *sigma)
{
   int tmax;
   double tau;
   
   if (noauto) return 0.5;
   
   (*lambda)=((n/10)>100)?100:(n/10);
   while(1)
   {
      tmax=(n+1-(*lambda))/2;
      error(tmax<=0,1,"acorr [read.ms6.c]","Unable to estimate autocorrelation");
      tau=tauint(n,a,tmax,(*lambda),w,sigma);
      if((*lambda)>(int)(3*tau)) break;
      (*lambda)=((int)(3*tau)>((*lambda)+10))?(int)(3*tau):((*lambda)+10);
   }
   
   return tau;
}

#define MIN_NBS 1000
static void print_stats(char *fin)
{
   int m,n,dt,t,t1,tmax,x0,lambda,w,fail;
   char base[NAME_SIZE],plt_file[NAME_SIZE],*p;
   FILE *fout;
   double *a,mean,err,tau,err_tau;
   double d[2];
   float *r;
   int nbs,bs;

   p=strstr(fin,".ms6.dat");
   error(p==NULL,1,"print_stats [read.ms6.c]","Unexpected data file name");
   n=p-fin;

   p=strrchr(fin,'/');
   if (p==NULL)
      p=fin;
   else
      p+=1;
   n-=(p-fin);

   error(n>=NAME_SIZE,1,"print_stats [read.ms6.c]","File name is too long");
   strncpy(base,p,n);
   base[n]='\0';

   error(name_size("stats/%s.stat.ms6.mpcac.dat",base)>=NAME_SIZE,1,
         "print_stats [read.ms6.c]","File name is too long");

   nbs=neff;
   if(nbs<MIN_NBS) nbs=MIN_NBS;
   a=malloc(nbs*sizeof(double));
   error(a==NULL,1,"print_stats [read.ms6.c]",
         "Unable to allocate data array");

   r=malloc(neff*sizeof(float));
   error(r==NULL,1,"print_stats [read.ms6.c]",
         "Unable to allocate random number array");

   x0=file_head.x0;
   tmax=file_head.tmax;

   printf("Analysis results printed to files:\n");
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.P.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable P\n");
   fprintf(fout,"# ---------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): < P(%d) P(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs            err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (t=0;t<tmax;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].P[t];
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",t,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.P.dat\n",base);
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.A0.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A1\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): < P(%d) A0(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs             err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (t=0;t<tmax;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].A0[t];
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e    %.1e    %6d   %6d\n",t,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.A0.dat\n",base);
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.A1.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A1\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): < P(%d) A1(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs             err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (t=0;t<tmax;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].A1[t];
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e    %.1e    %6d   %6d\n",t,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.A1.dat\n",base);
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.A2.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A2\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): < P(%d) A2(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs             err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (t=0;t<tmax;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].A2[t];
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e    %.1e    %6d   %6d\n",t,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.A2.dat\n",base);
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.A3.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable A3\n");
   fprintf(fout,"# ----------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): < P(%d) A3(t) >\n",x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs             err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (t=0;t<tmax;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].A3[t];
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e    %.1e    %6d   %6d\n",t,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.A3.dat\n",base);
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.mpcac.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable mpcac\n");
   fprintf(fout,"# -------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): < P(%d) [ A0(t+1)-A0(t-1) ] > / 4*< P(%d) P(t) >\n",x0,x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs             err        tau       neff\n");
   fprintf(fout,"#\n");

   for (t=1;t<tmax-1;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].P[t];
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      d[1]=4.0*average(neff,a);
      
      for (n=0;n<neff;n++)
         a[n]=-adat[nfirst+n].A0[t+1]+adat[nfirst+n].A0[t-1];
      d[0]=acorr(neff,a,&lambda,&w,&err_tau);
      if(d[0]>tau) tau=d[0];
      d[0]=average(neff,a);
      
      mean=d[0]/d[1];
   
      nbs=neff;
      if(nbs<MIN_NBS) nbs=MIN_NBS;
      for (bs=0;bs<nbs;bs++)
      {
         ranlxs(r,neff);
         d[0]=d[1]=0.0;
         for (m=0;m<neff;m++)
         {
            n=((int)(r[m]*(neff)));
            d[0]+=-adat[nfirst+n].A0[t+1]+adat[nfirst+n].A0[t-1];
            d[1]+=4.0*adat[nfirst+n].P[t];
         }
         a[bs]=d[0]/d[1];
      }
      
      err=sigma0(nbs,a)*sqrt(nbs-1.0)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e   %d\n",t,mean,err,tau,neff/(int)ceil(2.0*tau));
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.mpcac.dat\n",base);
   
   
   sprintf(plt_file,"stats/%s.stat.ms6.mpich.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable mpi with cosh^(-1)\n");
   fprintf(fout,"# --------------------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# IMPORTANT! This is a good choice only for periodic b.c.'s in time!\n");
   fprintf(fout,"# IMPORTANT! Exact zeros in the observable and/or the error mean that\n");
   fprintf(fout,"# IMPORTANT! the cosh could not be inverted!\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): cosh^(-1) < P(%d) [ P(t+1)+P(t-1) ] > / 2*< P(%d) P(t) >\n",x0,x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  dt   obs             err        tau       neff\n");
   fprintf(fout,"#\n");

   for (dt=0;dt<tmax/2;dt++)
   {
      t=(x0+dt)%tmax;
      t1=(x0-dt+tmax)%tmax;
      
      for (n=0;n<neff;n++)
         a[n]=(adat[nfirst+n].P[t]
               +adat[nfirst+n].P[t1]);
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      d[1]=average(neff,a);
      
      for (n=0;n<neff;n++)
         a[n]=0.5*(adat[nfirst+n].P[(t+1)%tmax]+adat[nfirst+n].P[(t-1+tmax)%tmax]
                   +adat[nfirst+n].P[(t1+1)%tmax]+adat[nfirst+n].P[(t1-1+tmax)%tmax]);
      d[0]=acorr(neff,a,&lambda,&w,&err_tau);
      if(d[0]>tau) tau=d[0];
      d[0]=average(neff,a);
      
      fail=0;
      mean=d[0]/d[1];
      if (mean<1.0)
         fail=1;
      else
         mean=log(mean+sqrt(mean*mean-1));

      nbs=neff;
      if(nbs<MIN_NBS) nbs=MIN_NBS;
      for (bs=0;(bs<nbs)&&(!fail);bs++)
      {
         ranlxs(r,neff);
         d[0]=d[1]=0.0;
         for (m=0;m<neff;m++)
         {
            n=((int)(r[m]*(neff)));
            d[0]+=0.5*(adat[nfirst+n].P[(t+1)%tmax]+adat[nfirst+n].P[(t-1+tmax)%tmax]
                       +adat[nfirst+n].P[(t1+1)%tmax]+adat[nfirst+n].P[(t1-1+tmax)%tmax]);
            d[1]+=(adat[nfirst+n].P[t]
                   +adat[nfirst+n].P[t1]);
         }
         a[bs]=d[0]/d[1];
         if (a[bs]<1.0)
            fail=1;
         else
            a[bs]=log(a[bs]+sqrt(a[bs]*a[bs]-1));
      }
      
      if (fail)
         err=0.0;
      else
         err=sigma0(nbs,a)*sqrt(nbs-1.0)*sqrt(2.0*tau);

      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e   %d\n",dt,mean,err,tau,neff/(int)ceil(2.0*tau));
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.mpich.dat\n",base);


   sprintf(plt_file,"stats/%s.stat.ms6.mpish.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms6.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms6, observable mpi with sinh^(-1)\n");
   fprintf(fout,"# --------------------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs(t): -sinh^(-1) < P(%d) [ P(t+1)-P(t-1) ] > / 2*< P(%d) P(t) >\n",x0,x0);
   fprintf(fout,"#\n");
   fprintf(fout,"#  t    obs             err        tau       neff\n");
   fprintf(fout,"#\n");

   for (t=1;t<tmax-1;t++)
   {
      for (n=0;n<neff;n++)
         a[n]=adat[nfirst+n].P[t];
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      d[1]=average(neff,a);
      
      for (n=0;n<neff;n++)
         a[n]=0.5*(adat[nfirst+n].P[t+1]-adat[nfirst+n].P[t-1]);
      d[0]=acorr(neff,a,&lambda,&w,&err_tau);
      if(d[0]>tau) tau=d[0];
      d[0]=average(neff,a);
      
      mean=-d[0]/d[1];
      mean=log(mean+sqrt(1.0+mean*mean));

      nbs=neff;
      if(nbs<MIN_NBS) nbs=MIN_NBS;
      for (bs=0;bs<nbs;bs++)
      {
         ranlxs(r,neff);
         d[0]=d[1]=0.0;
         for (m=0;m<neff;m++)
         {
            n=((int)(r[m]*(neff)));
            d[0]+=0.5*(adat[nfirst+n].P[t+1]-adat[nfirst+n].P[t-1]);
            d[1]+=adat[nfirst+n].P[t];
         }
         a[bs]=-d[0]/d[1];
         a[bs]=log(a[bs]+sqrt(1.0+a[bs]*a[bs]));
      }
      
      err=sigma0(nbs,a)*sqrt(nbs-1.0)*sqrt(2.0*tau);
      fprintf(fout,"  %3d   %+.6e   %.2e   %.1e   %d\n",t,mean,err,tau,neff/(int)ceil(2.0*tau));
   }

   fclose(fout);

   printf(" stats/%s.stat.ms6.mpish.dat\n",base);


   printf("\n");
   
   free(a);
}


int main(int argc,char *argv[])
{
   error(argc<2,1,"main [read.ms6.c]",
         "Syntax: read.ms6 [-r <first> <last> | -a] [-p] [-noauto] <filename>\n"
         "        -r      first and last configurations to be used for analysis\n"
         "        -a      all configurations are to be used for analysis\n"
         "        -noauto autocorrelation time is assumed to be 0.5\n");

   range=find_opt(argc,argv,"-r");

   if (range!=0)
   {
      error(sscanf(argv[range+1],"%d",&nfirst)!=1,1,"main [read.ms6.c]",
            "Syntax: read.ms6 [-r <first> <last> | -a] [-p] [-noauto] <filename>\n"
            "        -r      first and last configurations to be used for analysis\n"
            "        -a      all configurations are to be used for analysis\n"
            "        -noauto autocorrelation time is assumed to be 0.5\n");
      error(sscanf(argv[range+2],"%d",&nlast)!=1,1,"main [read.ms6.c]",
            "Syntax: read.ms6 [-r <first> <last> | -a] [-p] [-noauto] <filename>\n"
            "        -r      first and last configurations to be used for analysis\n"
            "        -a      all configurations are to be used for analysis\n"
            "        -noauto autocorrelation time is assumed to be 0.5\n");
   }
   
   all=find_opt(argc,argv,"-a");

   error((all!=0)&&(range!=0),1,"main [read.ms6.c]",
         "Syntax: read.ms6 [-r <first> <last> | -a] [-p] [-noauto] <filename>\n"
         "        -r      first and last configurations to be used for analysis\n"
         "        -a      all configurations are to be used for analysis\n"
         "        -noauto autocorrelation time is assumed to be 0.5\n");
   
   noauto=find_opt(argc,argv,"-noauto");

   rlxs_init(0,time(NULL));

   printf("\n");
   printf("Data written by the program ms6\n");
   printf("-------------------------------\n\n");

   read_file(argv[argc-1]);

   print_plots(argv[argc-1]);

   select_range();

   print_stats(argv[argc-1]);
   
   exit(0);
}

