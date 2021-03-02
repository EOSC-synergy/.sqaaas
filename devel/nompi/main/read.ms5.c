
/*******************************************************************************
*
* File read.ms5.c
*
* Copyright (C) 2017 Agostino Patella
* 
* Based on openQCD-1.6/devel/nompi/main/read1.c
* Copyright (C) 2010-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Reads and evaluates data from the *.ms5.dat files created by the programs ms5,
* qcd1 and ym1. The file to be read has to be specified on the command line.
*
* This program writes the history of the SU(3) Wilson flow observables to the
* file <run name>.plot.ms5.<obs>.dat in the plots directory, and their expectation
* values with error and integrated autocorrelation time to the file
* <run name>.stat.ms5.<obs>.dat in the stats directory.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "extras.h"

static struct
{
   int dn,nn,tmax;
   double eps;
} file_head;

typedef struct
{
   int nt;
   double **Usl,**Msl,**Fsl[6];
} dat_t;



static int nms,nfirst,nlast,neff;
static dat_t *adat;


static void read_file(char *fin)
{
   long ipos;
   FILE *fdat;
   int ir,in,id,t,nn,tmax,endian,pl;
   stdint_t istd[3];
   double dstd[1];
   double **pp,*p;
   static double *dbuf=NULL;


   fdat=fopen(fin,"rb");
   error(fdat==NULL,1,"read_file [read.ms5.c]","Unable to open data file");

   printf("Read data from file %s\n\n",fin);
   

   /* Read header */

   endian=endianness();

   ir=fread(istd,sizeof(stdint_t),3,fdat);
   ir+=fread(dstd,sizeof(double),1,fdat);
   error_root(ir!=4,1,"read_file [read.ms5.c]",
              "Incorrect file header");

   if (endian==BIG_ENDIAN)
   {
      bswap_int(3,istd);
      bswap_double(1,dstd);
   }

   file_head.dn=(int)istd[0];
   file_head.nn=(int)istd[1];
   file_head.tmax=(int)istd[2];
   file_head.eps=dstd[0];

   error_root(file_head.dn<=0,1,"read_file [read.ms5.c]","Incorrect value for dn");
   error_root(file_head.nn<0,1,"read_file [read.ms5.c]","Incorrect value for nn");
   error_root(file_head.tmax<=0,1,"read_file [read.ms5.c]","Incorrect value for tmax");
   error_root(file_head.eps<=0.0,1,"read_file [read.ms5.c]","Incorrect value for eps");


   /* Count measures */

   nn=file_head.nn;
   tmax=file_head.tmax;
   dbuf=amalloc(8*(nn+1)*tmax*sizeof(double),4);
   error((dbuf==NULL),1,"read_file [read.ms5.c]",
         "Unable to allocate buffer");

   ipos=ftell(fdat);
   nms=0;
   
   while(1)
   {
      if (fread(istd,sizeof(stdint_t),1,fdat)!=1)
         break;
      if (fread(dbuf,sizeof(double),8*(nn+1)*tmax,fdat)!=8*(nn+1)*tmax)
         break;
      nms++;
   }

   error(nms==0,1,"read_file [read.ms5.c]",
         "Empty data file");


   /* Allocate data array */

   adat=amalloc(nms*sizeof(*adat),3);
   pp=amalloc(8*(nn+1)*nms*sizeof(*pp),3);
   p=amalloc(8*(nn+1)*tmax*nms*sizeof(*p),4);

   error((adat==NULL)||(pp==NULL)||(p==NULL),1,"alloc_adat [read.ms5.c]",
         "Unable to allocate data arrays");
   
   for (id=0;id<nms;id++)
   {
      adat[id].Usl=pp;
      adat[id].Msl=pp+nn+1;
      adat[id].Fsl[0]=pp+2*(nn+1);
      adat[id].Fsl[1]=pp+3*(nn+1);
      adat[id].Fsl[2]=pp+4*(nn+1);
      adat[id].Fsl[3]=pp+5*(nn+1);
      adat[id].Fsl[4]=pp+6*(nn+1);
      adat[id].Fsl[5]=pp+7*(nn+1);

      for (in=0;in<(8*(nn+1));in++)
      {
         *pp=p;
         pp+=1;
         p+=tmax;
      }
   }

   
   /* Read data */
   
   fseek(fdat,ipos,SEEK_SET);

   for (id=0;id<nms;id++)
   {
      ir=fread(istd,sizeof(stdint_t),1,fdat);
      error(ir!=1,1,"read_file [read.ms5.c]",
            "Error while reading data file");

      if (endian==BIG_ENDIAN)
         bswap_int(1,istd);
      
      adat[id].nt=(int)istd[0];
      
      ir=fread(dbuf,sizeof(double),8*(nn+1)*tmax,fdat);
      error(ir!=8*(nn+1)*tmax,1,"read_file [read.ms5.c]",
            "Error while reading data file");

      if (endian==BIG_ENDIAN)
         bswap_double(8*(nn+1)*tmax,dbuf);
      
      p=dbuf;
      for (in=0;in<=nn;in++)
      {
         for (t=0;t<tmax;t++)
         {
            adat[id].Usl[in][t]=(*p);
            p++;
         }
      }

      for (in=0;in<=nn;in++)
      {
         for (t=0;t<tmax;t++)
         {
            adat[id].Msl[in][t]=(*p);
            p++;
         }
      }

      for (pl=0;pl<6;pl++)
      {
         for (in=0;in<=nn;in++)
         {
            for (t=0;t<tmax;t++)
            {
               adat[id].Fsl[pl][in][t]=(*p);
               p++;
            }
         }
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

   error(neff<2,1,"select_range [read.ms5.c]",
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

   error(ie!=0,1,"select_range [read.ms5.c]",
         "Varying trajectory number separation in selected range");
}


static void print_plots(char *fin)
{
   int n,ims,t,in,tmax,nn,dn,pl;
   char spl[6][3]={"01","02","03","23","31","12"};
   char base[NAME_SIZE],plt_file[NAME_SIZE],*p;
   double sum,eps;
   dat_t *ndat;
   FILE *fout;

   p=strstr(fin,".ms5.dat");
   error(p==NULL,1,"print_plots [read.ms5.c]","Unexpected data file name");
   n=p-fin;

   p=strrchr(fin,'/');
   if (p==NULL)
      p=fin;
   else
      p+=1;
   n-=(p-fin);

   error(n>=NAME_SIZE,1,"print_plots [read.ms5.c]","File name is too long");
   strncpy(base,p,n);
   base[n]='\0';

   error(name_size("plots/%s.plot.ms5.Fsl0.999.dat",base)>=NAME_SIZE,1,
         "print_plots [read.ms5.c]","File name is too long");

   nn=file_head.nn;
   dn=file_head.dn;
   eps=file_head.eps;
   tmax=file_head.tmax;

   printf("Data printed to files:\n");
   
   for (t=0;t<tmax;t++)
   {
      sprintf(plt_file,"plots/%s.plot.ms5.Usl.%03d.dat",base,t);
      fout=fopen(plt_file,"w");
      error(fout==NULL,1,"print_plots [read.ms5.c]",
            "Unable to open output file. Make sure that the 'plots' directory exists.");

      fprintf(fout,"#\n");
      fprintf(fout,"# Data written by the program ms5, observable Usl\n");
      fprintf(fout,"# -----------------------------------------------\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# Number of measurements = %d\n",nms);
      fprintf(fout,"#\n");
      fprintf(fout,"# nt:       trajectory number\n");
      fprintf(fout,"# Usl(x0):  U(1) Wilson action density, summed over time slice x0\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# column no        flowtime\n");
      for (in=0;in<=nn;in++)
         fprintf(fout,"# %4d            %.4e\n",in+2,in*dn*eps);
      fprintf(fout,"#\n");
      fprintf(fout,"#  nt          Usl(%d)...\n",t);
      fprintf(fout,"#\n");

      ndat=adat;

      for (ims=0;ims<nms;ims++)
      {
         fprintf(fout," %6d       ",(*ndat).nt);
         for (in=0;in<=nn;in++)
            fprintf(fout," %.4e  ",(*ndat).Usl[in][t]);
         fprintf(fout,"\n");
         ndat+=1;
      }

      fclose(fout);
   }

   printf(" plots/%s.plot.ms5.Usl.*.dat\n",base);
   
   for (t=0;t<tmax;t++)
   {
      sprintf(plt_file,"plots/%s.plot.ms5.Msl.%03d.dat",base,t);
      fout=fopen(plt_file,"w");
      error(fout==NULL,1,"print_plots [read.ms5.c]",
            "Unable to open output file. Make sure that the 'plots' directory exists.");

      fprintf(fout,"#\n");
      fprintf(fout,"# Data written by the program ms5, observable Msl\n");
      fprintf(fout,"# -----------------------------------------------\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# Number of measurements = %d\n",nms);
      fprintf(fout,"#\n");
      fprintf(fout,"# nt:       trajectory number\n");
      fprintf(fout,"# Msl(x0):  Maxwell action density, summed over time slice x0\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# column no        flowtime\n");
      for (in=0;in<=nn;in++)
         fprintf(fout,"# %4d            %.4e\n",in+2,in*dn*eps);
      fprintf(fout,"#\n");
      fprintf(fout,"#  nt          Msl(%d)...\n",t);
      fprintf(fout,"#\n");

      ndat=adat;

      for (ims=0;ims<nms;ims++)
      {
         fprintf(fout," %6d       ",(*ndat).nt);
         for (in=0;in<=nn;in++)
            fprintf(fout," %.4e  ",(*ndat).Msl[in][t]);
         fprintf(fout,"\n");
         ndat+=1;
      }

      fclose(fout);
   }
   
   printf(" plots/%s.plot.ms5.Msl.*.dat\n",base);
   
   for (pl=0;pl<6;pl++)
   {
      for (t=0;t<tmax;t++)
      {
         sprintf(plt_file,"plots/%s.plot.ms5.Fsl%d.%03d.dat",base,pl,t);
         fout=fopen(plt_file,"w");
         error(fout==NULL,1,"print_plots [read.ms5.c]",
               "Unable to open output file. Make sure that the 'plots' directory exists.");

         fprintf(fout,"#\n");
         fprintf(fout,"# Data written by the program ms5, observable Fsl[%d]\n",pl);
         fprintf(fout,"# ---------------------------------------------------\n");
         fprintf(fout,"#\n");
         fprintf(fout,"# Number of measurements = %d\n",nms);
         fprintf(fout,"#\n");
         fprintf(fout,"# nt:           trajectory number\n");
         fprintf(fout,"# Fsl[%d](x0):  Electromagnetic flux density across plane %s, summed over time slice x0\n",pl,spl[pl]);
         fprintf(fout,"#\n");
         fprintf(fout,"# column no        flowtime\n");
         for (in=0;in<=nn;in++)
            fprintf(fout,"# %4d            %.4e\n",in+2,in*dn*eps);
         fprintf(fout,"#\n");
         fprintf(fout,"#  nt          Fsl[%d](%d)...\n",pl,t);
         fprintf(fout,"#\n");

         ndat=adat;

         for (ims=0;ims<nms;ims++)
         {
            fprintf(fout," %6d       ",(*ndat).nt);
            for (in=0;in<=nn;in++)
               fprintf(fout," %.4e  ",(*ndat).Fsl[pl][in][t]);
            fprintf(fout,"\n");
            ndat+=1;
         }

         fclose(fout);
      }
   
      printf(" plots/%s.plot.ms5.Fsl%d.*.dat\n",base,pl);
   }

   sprintf(plt_file,"plots/%s.plot.ms5.Utot.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms5.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms5, observable Utot\n");
   fprintf(fout,"# ------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# Utot:     U(1) Wilson action\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# column no        flowtime\n");
   for (in=0;in<=nn;in++)
      fprintf(fout,"# %4d            %.4e\n",in+2,in*dn*eps);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          Utot...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (in=0;in<=nn;in++)
      {
         sum=0.0;
         for (t=0;t<tmax;t++)
            sum+=(*ndat).Usl[in][t];
         fprintf(fout," %.4e  ",sum);
      }
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms5.Utot.dat\n",base);

   sprintf(plt_file,"plots/%s.plot.ms5.Mtot.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_plots [read.ms5.c]",
         "Unable to open output file. Make sure that the 'plots' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms5, observable Mtot\n");
   fprintf(fout,"# ------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# Number of measurements = %d\n",nms);
   fprintf(fout,"#\n");
   fprintf(fout,"# nt:       trajectory number\n");
   fprintf(fout,"# Mtot:     Maxwell action\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# column no        flowtime\n");
   for (in=0;in<=nn;in++)
      fprintf(fout,"# %4d            %.4e\n",in+2,in*dn*eps);
   fprintf(fout,"#\n");
   fprintf(fout,"#  nt          Mtot...\n");
   fprintf(fout,"#\n");

   ndat=adat;

   for (ims=0;ims<nms;ims++)
   {
      fprintf(fout," %6d       ",(*ndat).nt);
      for (in=0;in<=nn;in++)
      {
         sum=0.0;
         for (t=0;t<tmax;t++)
            sum+=(*ndat).Msl[in][t];
         fprintf(fout," %.4e  ",sum);
      }
      fprintf(fout,"\n");
      ndat+=1;
   }

   fclose(fout);

   printf(" plots/%s.plot.ms5.Mtot.dat\n",base);

   for (pl=0;pl<6;pl++)
   {
      sprintf(plt_file,"plots/%s.plot.ms5.Flux%d.dat",base,pl);
      fout=fopen(plt_file,"w");
      error(fout==NULL,1,"print_plots [read.ms5.c]",
            "Unable to open output file. Make sure that the 'plots' directory exists.");

      fprintf(fout,"#\n");
      fprintf(fout,"# Data written by the program ms5, observable Flux[%d]\n",pl);
      fprintf(fout,"# ----------------------------------------------------\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# Number of measurements = %d\n",nms);
      fprintf(fout,"#\n");
      fprintf(fout,"# nt:       trajectory number\n");
      fprintf(fout,"# Flux[%d]: Electromagnetic flux density across plane %s\n",pl,spl[pl]);
      fprintf(fout,"#\n");
      fprintf(fout,"# column no        flowtime\n");
      for (in=0;in<=nn;in++)
         fprintf(fout,"# %4d            %.4e\n",in+2,in*dn*eps);
      fprintf(fout,"#\n");
      fprintf(fout,"#  nt          Flux[%d]...\n",pl);
      fprintf(fout,"#\n");

      ndat=adat;

      for (ims=0;ims<nms;ims++)
      {
         fprintf(fout," %6d       ",(*ndat).nt);
         for (in=0;in<=nn;in++)
         {
            sum=0.0;
            for (t=0;t<tmax;t++)
               sum+=(*ndat).Fsl[pl][in][t];
            fprintf(fout," %.4e  ",sum);
         }
         fprintf(fout,"\n");
         ndat+=1;
      }

      fclose(fout);

      printf(" plots/%s.plot.ms5.Flux[%d].dat\n",base,pl);
   }

   printf("\n");
}


static double acorr(int n,double *a,int *lambda,int *w,double *sigma)
{
   int tmax;
   double tau;
   
   (*lambda)=((n/10)>100)?100:(n/10);
   while(1)
   {
      tmax=(n+1-(*lambda))/2;
      error(tmax<=0,1,"acorr [read.ms5.c]","Unable to estimate autocorrelation");
      tau=tauint(n,a,tmax,(*lambda),w,sigma);
      if((*lambda)>(int)(3*tau)) break;
      (*lambda)=((int)(3*tau)>((*lambda)+10))?(int)(3*tau):((*lambda)+10);
   }
   
   return tau;
}


static void print_stats(char *fin)
{
   int n,t,in,tmax,nn,dn,lambda,w,pl;
   char spl[6][3]={"01","02","03","23","31","12"};
   char base[NAME_SIZE],plt_file[NAME_SIZE],*p;
   FILE *fout;
   double *a,mean,err,tau,err_tau,eps;

   p=strstr(fin,".ms5.dat");
   error(p==NULL,1,"print_stats [read.ms5.c]","Unexpected data file name");
   n=p-fin;

   p=strrchr(fin,'/');
   if (p==NULL)
      p=fin;
   else
      p+=1;
   n-=(p-fin);

   error(n>=NAME_SIZE,1,"print_stats [read.ms5.c]","File name is too long");
   strncpy(base,p,n);
   base[n]='\0';

   error(name_size("stats/%s.stat.ms5.Flux0sq.999.dat",base)>=NAME_SIZE,1,
         "print_stats [read.ms5.c]","File name is too long");

   a=malloc(neff*sizeof(double));
   error(a==NULL,1,"print_stats [read.ms5.c]",
         "Unable to allocate data array");

   nn=file_head.nn;
   dn=file_head.dn;
   eps=file_head.eps;
   tmax=file_head.tmax;

   printf("Analysis results printed to files:\n");
   
   for (t=0;t<tmax;t++)
   {
      sprintf(plt_file,"stats/%s.stat.ms5.Usl.%03d.dat",base,t);
      fout=fopen(plt_file,"w");
      error(fout==NULL,1,"print_stats [read.ms5.c]",
            "Unable to open output file. Make sure that the 'stats' directory exists.");

      fprintf(fout,"#\n");
      fprintf(fout,"# Data written by the program ms5, observable Usl\n");
      fprintf(fout,"# -----------------------------------------------\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
      fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
      fprintf(fout,"# Number of measurements = %d\n",neff);
      fprintf(fout,"#\n");
      fprintf(fout,"# obs: U(1) Wilson action density, summed over time slice %d\n",t);
      fprintf(fout,"#\n");
      fprintf(fout,"# flowtime     obs            err        tau        err(tau)   window   lambda\n");
      fprintf(fout,"#\n");

      for (in=0;in<=nn;in++)
      {
         for (n=0;n<neff;n++)
            a[n]=adat[nfirst+n].Usl[in][t];
         
         tau=acorr(neff,a,&lambda,&w,&err_tau);
         mean=average(neff,a);
         err=sigma0(neff,a)*sqrt(2.0*tau);
         fprintf(fout,"  %.4e   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",in*dn*eps,mean,err,tau,err_tau,w,lambda);
      }

      fclose(fout);
   }

   printf(" stats/%s.stat.ms5.Usl.*.dat\n",base);
   
   for (t=0;t<tmax;t++)
   {
      sprintf(plt_file,"stats/%s.stat.ms5.Msl.%03d.dat",base,t);
      fout=fopen(plt_file,"w");
      error(fout==NULL,1,"print_stats [read.ms5.c]",
            "Unable to open output file. Make sure that the 'stats' directory exists.");

      fprintf(fout,"#\n");
      fprintf(fout,"# Data written by the program ms5, observable Msl\n");
      fprintf(fout,"# -----------------------------------------------\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
      fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
      fprintf(fout,"# Number of measurements = %d\n",neff);
      fprintf(fout,"#\n");
      fprintf(fout,"# obs: Maxwell action density, summed over time slice %d\n",t);
      fprintf(fout,"#\n");
      fprintf(fout,"# flowtime     obs            err        tau        err(tau)   window   lambda\n");
      fprintf(fout,"#\n");

      for (in=0;in<=nn;in++)
      {
         for (n=0;n<neff;n++)
            a[n]=adat[nfirst+n].Msl[in][t];
         
         tau=acorr(neff,a,&lambda,&w,&err_tau);
         mean=average(neff,a);
         err=sigma0(neff,a)*sqrt(2.0*tau);
         fprintf(fout,"  %.4e   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",in*dn*eps,mean,err,tau,err_tau,w,lambda);
      }

      fclose(fout);
   }

   printf(" stats/%s.stat.ms5.Msl.*.dat\n",base);
   
   for (pl=0;pl<6;pl++)
   {
      for (t=0;t<tmax;t++)
      {
         sprintf(plt_file,"stats/%s.stat.ms5.Fsl%dsq.%03d.dat",base,pl,t);
         fout=fopen(plt_file,"w");
         error(fout==NULL,1,"print_stats [read.ms5.c]",
               "Unable to open output file. Make sure that the 'stats' directory exists.");

         fprintf(fout,"#\n");
         fprintf(fout,"# Data written by the program ms5, observable Fsl[%d]^2\n",pl);
         fprintf(fout,"# ------------------------------------------------\n");
         fprintf(fout,"#\n");
         fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
         fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
         fprintf(fout,"# Number of measurements = %d\n",neff);
         fprintf(fout,"#\n");
         fprintf(fout,"# obs: Electromagnetic flux density across plane %s, summed over time slice %d, then squared\n",spl[pl],t);
         fprintf(fout,"#\n");
         fprintf(fout,"# flowtime     obs            err        tau        err(tau)   window   lambda\n");
         fprintf(fout,"#\n");

         for (in=0;in<=nn;in++)
         {
            for (n=0;n<neff;n++)
            {
               a[n]=adat[nfirst+n].Fsl[pl][in][t];
               a[n]=a[n]*a[n];
            }
            
            tau=acorr(neff,a,&lambda,&w,&err_tau);
            mean=average(neff,a);
            err=sigma0(neff,a)*sqrt(2.0*tau);
            fprintf(fout,"  %.4e   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",in*dn*eps,mean,err,tau,err_tau,w,lambda);
         }

         fclose(fout);
      }

      printf(" stats/%s.stat.ms5.Fsl%dsq.*.dat\n",base,pl);
   }
   
   sprintf(plt_file,"stats/%s.stat.ms5.Utot.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms5.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms5, observable Utot\n");
   fprintf(fout,"# ------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs: U(1) Wilson action\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# flowtime     obs            err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (in=0;in<=nn;in++)
   {
      for (n=0;n<neff;n++)
      {
         a[n]=0.0;
         for (t=0;t<tmax;t++)
            a[n]+=adat[nfirst+n].Usl[in][t];
      }
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %.4e   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",in*dn*eps,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms5.Mtot.dat\n",base);
   
   sprintf(plt_file,"stats/%s.stat.ms5.Mtot.dat",base);
   fout=fopen(plt_file,"w");
   error(fout==NULL,1,"print_stats [read.ms5.c]",
         "Unable to open output file. Make sure that the 'stats' directory exists.");

   fprintf(fout,"#\n");
   fprintf(fout,"# Data written by the program ms5, observable Mtot\n");
   fprintf(fout,"# ------------------------------------------------\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
   fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
   fprintf(fout,"# Number of measurements = %d\n",neff);
   fprintf(fout,"#\n");
   fprintf(fout,"# obs: Maxwell action\n");
   fprintf(fout,"#\n");
   fprintf(fout,"# flowtime     obs            err        tau        err(tau)   window   lambda\n");
   fprintf(fout,"#\n");

   for (in=0;in<=nn;in++)
   {
      for (n=0;n<neff;n++)
      {
         a[n]=0.0;
         for (t=0;t<tmax;t++)
            a[n]+=adat[nfirst+n].Msl[in][t];
      }
      
      tau=acorr(neff,a,&lambda,&w,&err_tau);
      mean=average(neff,a);
      err=sigma0(neff,a)*sqrt(2.0*tau);
      fprintf(fout,"  %.4e   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",in*dn*eps,mean,err,tau,err_tau,w,lambda);
   }

   fclose(fout);

   printf(" stats/%s.stat.ms5.Mtot.dat\n",base);
   
   for (pl=0;pl<6;pl++)
   {
      sprintf(plt_file,"stats/%s.stat.ms5.Flux%dsq.dat",base,pl);
      fout=fopen(plt_file,"w");
      error(fout==NULL,1,"print_stats [read.ms5.c]",
            "Unable to open output file. Make sure that the 'stats' directory exists.");

      fprintf(fout,"#\n");
      fprintf(fout,"# Data written by the program ms5, observable Flux[%d]\n",pl);
      fprintf(fout,"# ------------------------------------------------\n");
      fprintf(fout,"#\n");
      fprintf(fout,"# First trajectory       = %d\n",adat[nfirst].nt);
      fprintf(fout,"# Last trajectory        = %d\n",adat[nlast-1].nt);
      fprintf(fout,"# Number of measurements = %d\n",neff);
      fprintf(fout,"#\n");
      fprintf(fout,"# obs: Squared electromagnetic flux across plane %s\n",spl[pl]);
      fprintf(fout,"#\n");
      fprintf(fout,"# flowtime     obs            err        tau        err(tau)   window   lambda\n");
      fprintf(fout,"#\n");

      for (in=0;in<=nn;in++)
      {
         for (n=0;n<neff;n++)
         {
            a[n]=0.0;
            for (t=0;t<tmax;t++)
               a[n]+=adat[nfirst+n].Fsl[pl][in][t];
            a[n]=a[n]*a[n];
         }
         
         tau=acorr(neff,a,&lambda,&w,&err_tau);
         mean=average(neff,a);
         err=sigma0(neff,a)*sqrt(2.0*tau);
         fprintf(fout,"  %.4e   %.6e   %.2e   %.1e    %.1e    %6d   %6d\n",in*dn*eps,mean,err,tau,err_tau,w,lambda);
      }

      fclose(fout);

      printf(" stats/%s.stat.ms5.Flux%dsq.dat\n",base,pl);
   }

   printf("\n");
   
   free(a);
}


int main(int argc,char *argv[])
{

   error(argc!=2,1,"main [read.ms5.c]","Syntax: read.ms5 <filename>");

   printf("\n");
   printf("Data written by the program ms5\n");
   printf("-------------------------------\n\n");

   read_file(argv[1]);

   print_plots(argv[1]);

   select_range();

   print_stats(argv[1]);
   
   exit(0);
}
