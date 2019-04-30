
/*******************************************************************************
*
* File iso1.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher, Isabel Campos
*               2017            Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* HMC simulation program for QCD+QED with Wilson quarks.
*
* Syntax: iso1 -i <filename> [-noloc] [-noexp] [-rmold] [-noms]
*                            [-c <filename> [-a [-norng]]]
*
* For usage instructions see the file README.iso1.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "flags.h"
#include "random.h"
#include "su3fcts.h"
#include "linalg.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "archive.h"
#include "forces.h"
#include "update.h"
#include "wflow.h"
#include "tcharge.h"
#include "u1ftensor.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

typedef struct
{
   int nt,iac;
   double dH,avpl3,avpl1;
} dat_t;

static struct
{
   int dn,nn,tmax;
   double eps;
} file_head;

static struct
{
   int nt;
   double **Wsl,**Ysl,**Qsl;
} data3;

static struct
{
   int nt;
   double **Usl,**Msl,**Fsl[6];
} data1;

static int my_rank,noloc,noexp,rmold,noms,norng;
static int scnfg,append,endian;
static int level,seed;
static int nth,ntr,dtr_log,dtr_ms,dtr_cnfg;
static int ipgrd[2],flint;
static double *Wact,*Yact,*Qtop;
static double *Uact,*Mact,*Flux[6];

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char ms3dat_file[NAME_SIZE],ms3dat_save[NAME_SIZE];
static char ms1dat_file[NAME_SIZE],ms1dat_save[NAME_SIZE];
static char rng_file[NAME_SIZE],rng_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],end_file[NAME_SIZE];
static char nbase[NAME_SIZE],cnfg[NAME_SIZE];
static FILE *flog;

static hmc_parms_t hmc;


static int write_dat(FILE *fdat,int n,dat_t *ndat)
{
   int i,iw,ic;
   stdint_t istd[2];
   double dstd[3];

   ic=0;

   for (i=0;i<n;i++)
   {
      istd[0]=(stdint_t)((*ndat).nt);
      istd[1]=(stdint_t)((*ndat).iac);

      dstd[0]=(*ndat).dH;
      dstd[1]=(*ndat).avpl3;
      dstd[2]=(*ndat).avpl1;

      if (endian==BIG_ENDIAN)
      {
         bswap_int(2,istd);
         bswap_double(3,dstd);
      }

      iw=fwrite(istd,sizeof(stdint_t),2,fdat);
      iw+=fwrite(dstd,sizeof(double),3,fdat);

      if (iw!=5)
         return ic;

      ic+=1;
      ndat+=1;
   }

   return ic;
}


static int read_dat(FILE *fdat,int n,dat_t *ndat)
{
   int i,ir,ic;
   stdint_t istd[2];
   double dstd[3];

   ic=0;

   for (i=0;i<n;i++)
   {
      ir=fread(istd,sizeof(stdint_t),2,fdat);
      ir+=fread(dstd,sizeof(double),3,fdat);

      if (ir!=5)
         return ic;

      if (endian==BIG_ENDIAN)
      {
         bswap_int(2,istd);
         bswap_double(3,dstd);
      }

      (*ndat).nt=(int)(istd[0]);
      (*ndat).iac=(int)(istd[1]);

      (*ndat).dH=dstd[0];
      (*ndat).avpl3=dstd[1];
      (*ndat).avpl1=dstd[2];

      ic+=1;
      ndat+=1;
   }

   return ic;
}


static void alloc_data(void)
{
   int nn,tmax;
   int in;
   double **pp,*p;

   nn=file_head.nn;
   tmax=file_head.tmax;

   pp=amalloc(11*(nn+1)*sizeof(*pp),3);
   p=amalloc(11*(nn+1)*(tmax+1)*sizeof(*p),4);

   error((pp==NULL)||(p==NULL),1,"alloc_data [iso1.c]",
         "Unable to allocate data arrays");

   data3.Wsl=pp;
   data3.Ysl=pp+nn+1;
   data3.Qsl=pp+2*(nn+1);
   data1.Usl=pp+3*(nn+1);
   data1.Msl=pp+4*(nn+1);
   data1.Fsl[0]=pp+5*(nn+1);
   data1.Fsl[1]=pp+6*(nn+1);
   data1.Fsl[2]=pp+7*(nn+1);
   data1.Fsl[3]=pp+8*(nn+1);
   data1.Fsl[4]=pp+9*(nn+1);
   data1.Fsl[5]=pp+10*(nn+1);

   for (in=0;in<(11*(nn+1));in++)
   {
      *pp=p;
      pp+=1;
      p+=tmax;
   }

   Wact=p;
   p+=nn+1;
   Yact=p;
   p+=nn+1;
   Qtop=p;
   p+=nn+1;
   Uact=p;
   p+=nn+1;
   Mact=p;
   p+=nn+1;
   Flux[0]=p;
   p+=nn+1;
   Flux[1]=p;
   p+=nn+1;
   Flux[2]=p;
   p+=nn+1;
   Flux[3]=p;
   p+=nn+1;
   Flux[4]=p;
   p+=nn+1;
   Flux[5]=p;
}


static void write_file_head(FILE *fdat)
{
   int iw;
   stdint_t istd[3];
   double dstd[1];

   istd[0]=(stdint_t)(file_head.dn);
   istd[1]=(stdint_t)(file_head.nn);
   istd[2]=(stdint_t)(file_head.tmax);
   dstd[0]=file_head.eps;

   if (endian==BIG_ENDIAN)
   {
      bswap_int(3,istd);
      bswap_double(1,dstd);
   }

   iw=fwrite(istd,sizeof(stdint_t),3,fdat);
   iw+=fwrite(dstd,sizeof(double),1,fdat);

   error_root(iw!=4,1,"write_file_head [iso1.c]",
              "Incorrect write count");
}


static void check_file_head(FILE *fdat)
{
   int ir;
   stdint_t istd[3];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),3,fdat);
   ir+=fread(dstd,sizeof(double),1,fdat);

   error_root(ir!=4,1,"check_file_head [iso1.c]",
              "Incorrect read count");

   if (endian==BIG_ENDIAN)
   {
      bswap_int(3,istd);
      bswap_double(1,dstd);
   }

   error_root(((int)(istd[0])!=file_head.dn)||
              ((int)(istd[1])!=file_head.nn)||
              ((int)(istd[2])!=file_head.tmax)||
              (dstd[0]!=file_head.eps),1,"check_file_head [iso1.c]",
              "Unexpected value of dn,nn,tmax or eps");
}


static void write_data3(FILE *fdat)
{
   int iw,nn,tmax;
   int in,t;
   stdint_t istd[1];
   double dstd[1];

   istd[0]=(stdint_t)(data3.nt);

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   iw=fwrite(istd,sizeof(stdint_t),1,fdat);

   nn=file_head.nn;
   tmax=file_head.tmax;

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data3.Wsl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data3.Ysl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data3.Qsl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   error_root(iw!=(1+3*(nn+1)*tmax),1,"write_data3 [iso1.c]",
              "Incorrect write count");
}


static void write_data1(FILE *fdat)
{
   int iw,nn,tmax;
   int in,t,k;
   stdint_t istd[1];
   double dstd[1];

   istd[0]=(stdint_t)(data1.nt);

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   iw=fwrite(istd,sizeof(stdint_t),1,fdat);

   nn=file_head.nn;
   tmax=file_head.tmax;

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data1.Usl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         dstd[0]=data1.Msl[in][t];

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         iw+=fwrite(dstd,sizeof(double),1,fdat);
      }
   }

   for (k=0;k<6;k++)
   {
      for (in=0;in<=nn;in++)
      {
         for (t=0;t<tmax;t++)
         {
            dstd[0]=data1.Fsl[k][in][t];

            if (endian==BIG_ENDIAN)
               bswap_double(1,dstd);

            iw+=fwrite(dstd,sizeof(double),1,fdat);
         }
      }
   }

   error_root(iw!=(1+8*(nn+1)*tmax),1,"write_data1 [iso1.c]",
              "Incorrect write count");
}


static int read_data3(FILE *fdat)
{
   int ir,nn,tmax;
   int in,t;
   stdint_t istd[1];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);

   if (ir!=1)
      return 0;

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   data3.nt=(int)(istd[0]);

   nn=file_head.nn;
   tmax=file_head.tmax;

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data3.Wsl[in][t]=dstd[0];
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data3.Ysl[in][t]=dstd[0];
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data3.Qsl[in][t]=dstd[0];
      }
   }

   error_root(ir!=(1+3*(nn+1)*tmax),1,"read_data3 [iso1.c]",
              "Read error or incomplete data record");

   return 1;
}


static int read_data1(FILE *fdat)
{
   int ir,nn,tmax;
   int in,t,k;
   stdint_t istd[1];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);

   if (ir!=1)
      return 0;

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   data1.nt=(int)(istd[0]);

   nn=file_head.nn;
   tmax=file_head.tmax;

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data1.Usl[in][t]=dstd[0];
      }
   }

   for (in=0;in<=nn;in++)
   {
      for (t=0;t<tmax;t++)
      {
         ir+=fread(dstd,sizeof(double),1,fdat);

         if (endian==BIG_ENDIAN)
            bswap_double(1,dstd);

         data1.Msl[in][t]=dstd[0];
      }
   }

   for (k=0;k<6;k++)
   {
      for (in=0;in<=nn;in++)
      {
         for (t=0;t<tmax;t++)
         {
            ir+=fread(dstd,sizeof(double),1,fdat);

            if (endian==BIG_ENDIAN)
               bswap_double(1,dstd);

            data1.Fsl[k][in][t]=dstd[0];
         }
      }
   }

   error_root(ir!=(1+8*(nn+1)*tmax),1,"read_data1 [iso1.c]",
              "Read error or incomplete data record");

   return 1;
}


static void read_dirs(void)
{
   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);
      read_line("dat_dir","%s",dat_dir);
      if (noloc==0)
         read_line("loc_dir","%s",loc_dir);
      else
         loc_dir[0]='\0';
      if ((noexp==0)||((scnfg)&&(cnfg[strlen(cnfg)-1]!='*')))
         read_line("cnfg_dir","%s",cnfg_dir);
      else
         cnfg_dir[0]='\0';
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   if (noloc==0)
      check_dir(loc_dir);

   if (noexp==0)
      check_dir_root(cnfg_dir);

   check_dir_root(log_dir);
   check_dir_root(dat_dir);
   error_root(name_size("%s/%s.log~",log_dir,nbase)>=NAME_SIZE,1,
              "setup_files [iso1.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms3.dat~",dat_dir,nbase)>=NAME_SIZE,1,
              "setup_files [iso1.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.par",dat_dir,nbase);
   sprintf(dat_file,"%s/%s.dat",dat_dir,nbase);
   sprintf(ms3dat_file,"%s/%s.ms3.dat",dat_dir,nbase);
   sprintf(ms1dat_file,"%s/%s.ms5.dat",dat_dir,nbase);
   sprintf(rng_file,"%s/%s.rng",dat_dir,nbase);
   sprintf(end_file,"%s/%s.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
   sprintf(ms3dat_save,"%s~",ms3dat_file);
   sprintf(ms1dat_save,"%s~",ms1dat_file);
   sprintf(rng_save,"%s~",rng_file);
}


static void read_flds_bc_lat_parms(FILE *fpar)
{
   int nfl,ifl,bc,sf,cs,type,qhat;
   double phi[2],phi_prime[2];
   double beta,c0,cG,cG_prime;
   double alpha,lambda,invqel;
   double kappa,su3csw,u1csw,cF,cF_prime,th1,th2,th3;

   if (my_rank==0)
   {
      find_section("Quark action");
      read_line("nfl","%d",&nfl);
   }
   MPI_Bcast(&nfl,1,MPI_INT,0,MPI_COMM_WORLD);

   set_flds_parms(3,nfl);

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
         error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                    "Unknown time boundary condition type %s",line);
      
      read_line("cstar","%d",&cs);

      sf=0;
      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
      if ((cs==0)&&(bc==1))
         read_dprms("phi",2,phi);
      if ((bc==1)||(bc==2))
      {
         read_line("SFtype","%s",&line);
         if ((strcmp(line,"orbifold")==0)||
             (strcmp(line,"openQCD-1.4")==0)||(strcmp(line,"0")==0))
            sf=0;
         else if ((strcmp(line,"AFW-typeB")==0)||
                  (strcmp(line,"openQCD-1.2")==0)||(strcmp(line,"1")==0))
            sf=1;
         else
            error_root(1,1,"read_flds_bc_lat_parms [main.c]",
                       "Unknown SF type %s",line);
         
         if (cs==0)
            read_dprms("phi'",2,phi_prime);
      }
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&sf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_bc_parms(bc,sf,cs,phi,phi_prime);

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
   }
   
   MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_su3lat_parms(beta,c0,cG,cG_prime);

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
   }
   
   MPI_Bcast(&type,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&invqel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&lambda,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&c0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cG,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cG_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_u1lat_parms(type,alpha,invqel,lambda,c0,cG,cG_prime);

   for (ifl=0;ifl<nfl;ifl++)
   {
      cF=cF_prime=0.0;
      th1=th2=th3=0.0;

      sprintf(line,"Flavour %d",ifl);
      if (my_rank==0)
      {
         find_section(line);
         read_line("kappa","%lf",&kappa);
         read_line("qhat","%d",&qhat);
         read_line("su3csw","%lf",&su3csw);
         read_line("u1csw","%lf",&u1csw);
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

   if (append)
      check_flds_bc_lat_parms(fpar);
   else
      write_flds_bc_lat_parms(fpar);
}


static void read_schedule(FILE *fpar)
{
   int ie,ir,iw;
   stdint_t istd[3];

   if (my_rank==0)
   {
      find_section("MD trajectories");
      read_line("nth","%d",&nth);
      read_line("ntr","%d",&ntr);
      read_line("dtr_log","%d",&dtr_log);
      if (noms==0)
         read_line("dtr_ms","%d",&dtr_ms);
      else
         dtr_ms=0;
      read_line("dtr_cnfg","%d",&dtr_cnfg);

      error_root((append!=0)&&(nth!=0),1,"read_schedule [iso1.c]",
                 "Continuation run: nth must be equal to zero");

      ie=0;
      ie|=(nth<0);
      ie|=(ntr<1);
      ie|=(dtr_log<1);
      ie|=(dtr_log>dtr_cnfg);
      ie|=((dtr_cnfg%dtr_log)!=0);
      ie|=((nth%dtr_cnfg)!=0);
      ie|=((ntr%dtr_cnfg)!=0);

      if (noms==0)
      {
         ie|=(dtr_ms<dtr_log);
         ie|=(dtr_ms>dtr_cnfg);
         ie|=((dtr_ms%dtr_log)!=0);
         ie|=((dtr_cnfg%dtr_ms)!=0);
      }

      error_root(ie!=0,1,"read_schedule [iso1.c]",
                 "Improper value of nth,ntr,dtr_log,dtr_ms or dtr_cnfg");
   }

   MPI_Bcast(&nth,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ntr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dtr_log,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dtr_ms,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dtr_cnfg,1,MPI_INT,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      if (append)
      {
         ir=fread(istd,sizeof(stdint_t),3,fpar);
         error_root(ir!=3,1,"read_schedule [iso1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
            bswap_int(3,istd);

         ie=0;
         ie|=(istd[0]!=(stdint_t)(dtr_log));
         ie|=(istd[1]!=(stdint_t)(dtr_ms));
         ie|=(istd[2]!=(stdint_t)(dtr_cnfg));

         error_root(ie!=0,1,"read_schedule [iso1.c]",
                    "Parameters do not match previous run");
      }
      else
      {
         istd[0]=(stdint_t)(dtr_log);
         istd[1]=(stdint_t)(dtr_ms);
         istd[2]=(stdint_t)(dtr_cnfg);

         if (endian==BIG_ENDIAN)
            bswap_int(3,istd);

         iw=fwrite(istd,sizeof(stdint_t),3,fpar);
         error_root(iw!=3,1,"read_schedule [iso1.c]",
                    "Incorrect write count");
      }
   }
}


static void read_actions(FILE *fpar)
{
   int i,k,l,nact,*iact;
   int npf,nlv,nmu;
   double tau,*mu;
   action_parms_t ap;
   rat_parms_t rp;

   if (my_rank==0)
   {
      find_section("HMC parameters");
      nact=count_tokens("actions");
      read_line("npf","%d",&npf);
      read_line("nlv","%d",&nlv);
      read_line("tau","%lf",&tau);
   }

   MPI_Bcast(&nact,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&npf,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nlv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&tau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   if (nact>0)
   {
      iact=malloc(nact*sizeof(*iact));
      error(iact==NULL,1,"read_actions [iso1.c]",
            "Unable to allocate temporary array");
      if (my_rank==0)
         read_iprms("actions",nact,iact);
      MPI_Bcast(iact,nact,MPI_INT,0,MPI_COMM_WORLD);
   }
   else
      iact=NULL;

   nmu=0;

   for (i=0;i<nact;i++)
   {
      k=iact[i];
      ap=action_parms(k);

      if (ap.action==ACTIONS)
         read_action_parms(k);

      ap=action_parms(k);

      if ((ap.action==ACF_RAT)||(ap.action==ACF_RAT_SDET))
      {
         l=ap.irat[0];
         rp=rat_parms(l);

         if (rp.degree==0)
            read_rat_parms(l);
      }
      else if ((nmu==0)&&((ap.action==ACF_TM1)||
                          (ap.action==ACF_TM1_EO)||
                          (ap.action==ACF_TM1_EO_SDET)||
                          (ap.action==ACF_TM2)||
                          (ap.action==ACF_TM2_EO)))
      {
         if (my_rank==0)
         {
            find_section("HMC parameters");
            nmu=count_tokens("mu");
         }

         MPI_Bcast(&nmu,1,MPI_INT,0,MPI_COMM_WORLD);
      }
   }

   if (nmu>0)
   {
      mu=malloc(nmu*sizeof(*mu));
      error(mu==NULL,1,"read_actions [iso1.c]",
            "Unable to allocate temporary array");

      if (my_rank==0)
      {
         find_section("HMC parameters");
         read_dprms("mu",nmu,mu);
      }

      MPI_Bcast(mu,nmu,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
      mu=NULL;

   hmc=set_hmc_parms(nact,iact,npf,nmu,mu,nlv,tau);

   if (nact>0)
      free(iact);
   if (nmu>0)
      free(mu);

   if (append)
   {
      check_hmc_parms(fpar);
      check_action_parms(fpar);
   }
   else
   {
      write_hmc_parms(fpar);
      write_action_parms(fpar);
   }
}


static void read_integrator(FILE *fpar)
{
   int nlv,i,j,k,l;
   mdint_parms_t mdp;
   force_parms_t fp;
   rat_parms_t rp;

   nlv=hmc.nlv;

   for (i=0;i<nlv;i++)
   {
      read_mdint_parms(i);
      mdp=mdint_parms(i);

      for (j=0;j<mdp.nfr;j++)
      {
         k=mdp.ifr[j];
         fp=force_parms(k);

         if (fp.force==FORCES)
            read_force_parms2(k);

         fp=force_parms(k);

         if ((fp.force==FRF_RAT)||(fp.force==FRF_RAT_SDET))
         {
            l=fp.irat[0];
            rp=rat_parms(l);

            if (rp.degree==0)
               read_rat_parms(l);
         }
      }
   }

   if (append)
   {
      check_rat_parms(fpar);
      check_mdint_parms(fpar);
      check_force_parms(fpar);
   }
   else
   {
      write_rat_parms(fpar);
      write_mdint_parms(fpar);
      write_force_parms(fpar);
   }
}


static void read_sap_parms(FILE *fpar)
{
   int bs[4];

   if (my_rank==0)
   {
      find_section("SAP");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   set_sap_parms(bs,1,4,5);

   if (append)
      check_sap_parms(fpar);
   else
      write_sap_parms(fpar);
}


static void read_dfl_parms(FILE *fpar)
{
   int bs[4],Ns;
   int ninv,nmr,ncy,nkv,nmx,nsm,qhat;
   double kappa,mu,su3csw,u1csw,cF,cF_prime,th1,th2,th3,res,dtau;

   if (my_rank==0)
   {
      find_section("Deflation subspace");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
      read_line("Ns","%d",&Ns);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_parms(bs,Ns);

   if (my_rank==0)
   {
      find_section("Deflation subspace generation");
      read_line("kappa","%lf",&kappa);
      read_line("qhat","%d",&qhat);
      read_line("mu","%lf",&mu);
      read_line("su3csw","%lf",&su3csw);
      read_line("u1csw","%lf",&u1csw);
      if (bc_type()!=3) read_line("cF","%lf",&cF);
      if (bc_type()==2) read_line("cF'","%lf",&cF_prime);
      if (bc_cstar()==0) read_line("theta","%lf %lf %lf",&th1,&th2,&th3);
      read_line("ninv","%d",&ninv);
      read_line("nmr","%d",&nmr);
      read_line("ncy","%d",&ncy);
   }

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&qhat,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&su3csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&u1csw,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&cF_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&th1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&th2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&th3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_gen_parms(kappa,mu,qhat,su3csw,u1csw,cF,cF_prime,th1,th2,th3,ninv,nmr,ncy);

   if (my_rank==0)
   {
      find_section("Deflation projection");
      read_line("nkv","%d",&nkv);
      read_line("nmx","%d",&nmx);
      read_line("res","%lf",&res);
   }

   MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   set_dfl_pro_parms(nkv,nmx,res);

   if (my_rank==0)
   {
      find_section("Deflation update scheme");
      read_line("dtau","%lf",&dtau);
      read_line("nsm","%d",&nsm);
   }

   MPI_Bcast(&dtau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nsm,1,MPI_INT,0,MPI_COMM_WORLD);
   set_dfl_upd_parms(dtau,nsm);

   if (append)
      check_dfl_parms(fpar);
   else
      write_dfl_parms(fpar);
}


static void read_solvers(FILE *fpar)
{
   int nact,*iact,nlv;
   int nfr,*ifr,i,j,k;
   int nsp,isap,idfl;
   mdint_parms_t mdp;
   action_parms_t ap;
   force_parms_t fp;
   solver_parms_t sp;

   nact=hmc.nact;
   iact=hmc.iact;
   nlv=hmc.nlv;

   isap=0;
   idfl=0;

   for (i=0;i<nact;i++)
   {
      ap=action_parms(iact[i]);

      if ((ap.action==ACF_TM1)||
          (ap.action==ACF_TM1_EO)||
          (ap.action==ACF_TM1_EO_SDET)||
          (ap.action==ACF_TM2)||
          (ap.action==ACF_TM2_EO)||
          (ap.action==ACF_RAT)||
          (ap.action==ACF_RAT_SDET))
      {
         if ((ap.action==ACF_TM2)||
             (ap.action==ACF_TM2_EO))
            nsp=2;
         else
            nsp=1;

         for (k=0;k<nsp;k++)
         {
            j=ap.isp[k];
            sp=solver_parms(j);

            if (sp.solver==SOLVERS)
            {
               read_solver_parms(j);
               sp=solver_parms(j);

               if (sp.solver==SAP_GCR)
                  isap=1;
               else if (sp.solver==DFL_SAP_GCR)
               {
                  isap=1;
                  idfl=1;
               }
            }
         }
      }
   }

   for (i=0;i<nlv;i++)
   {
      mdp=mdint_parms(i);
      nfr=mdp.nfr;
      ifr=mdp.ifr;

      for (j=0;j<nfr;j++)
      {
         fp=force_parms(ifr[j]);

         if ((fp.force==FRF_TM1)||
             (fp.force==FRF_TM1_EO)||
             (fp.force==FRF_TM1_EO_SDET)||
             (fp.force==FRF_TM2)||
             (fp.force==FRF_TM2_EO)||
             (fp.force==FRF_RAT)||
             (fp.force==FRF_RAT_SDET))
         {
            k=fp.isp[0];
            sp=solver_parms(k);

            if (sp.solver==SOLVERS)
            {
               read_solver_parms(k);
               sp=solver_parms(k);

               if (sp.solver==SAP_GCR)
                  isap=1;
               else if (sp.solver==DFL_SAP_GCR)
               {
                  isap=1;
                  idfl=1;
               }
            }
         }
      }
   }

   if (append)
      check_solver_parms(fpar);
   else
      write_solver_parms(fpar);

   if (isap)
      read_sap_parms(fpar);

   if (idfl)
      read_dfl_parms(fpar);
}


static void read_wflow_parms(FILE *fpar)
{
   int nstep,dnms,ie,ir,iw;
   stdint_t istd[3];
   double eps,dstd[1];

   if (my_rank==0)
   {
      if (append)
      {
         ir=fread(istd,sizeof(stdint_t),1,fpar);
         error_root(ir!=1,1,"read_wflow_parms [iso1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
            bswap_int(1,istd);

         error_root(istd[0]!=(stdint_t)(noms==0),1,"read_wflow_parms [iso1.c]",
                    "Attempt to mix measurement with other runs");
      }
      else
      {
         istd[0]=(stdint_t)(noms==0);

         if (endian==BIG_ENDIAN)
            bswap_int(1,istd);

         iw=fwrite(istd,sizeof(stdint_t),1,fpar);
         error_root(iw!=1,1,"read_wflow_parms [iso1.c]",
                    "Incorrect write count");
      }

      if (noms==0)
      {
         find_section("Wilson flow");
         read_line("integrator","%s",line);
         read_line("eps","%lf",&eps);
         read_line("nstep","%d",&nstep);
         read_line("dnms","%d",&dnms);

         if (strcmp(line,"EULER")==0)
            flint=0;
         else if (strcmp(line,"RK2")==0)
            flint=1;
         else if (strcmp(line,"RK3")==0)
            flint=2;
         else
            error_root(1,1,"read_wflow_parms [iso1.c]","Unknown integrator");

         error_root((dnms<1)||(nstep<dnms)||((nstep%dnms)!=0),1,
                    "read_wflow_parms [iso1.c]",
                    "nstep must be a multiple of dnms");
      }
      else
      {
         flint=0;
         eps=0.0;
         nstep=1;
         dnms=1;
      }
   }

   MPI_Bcast(&flint,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&eps,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dnms,1,MPI_INT,0,MPI_COMM_WORLD);

   file_head.dn=dnms;
   file_head.nn=nstep/dnms;
   file_head.tmax=N0;
   file_head.eps=eps;

   if ((my_rank==0)&&(noms==0))
   {
      if (append)
      {
         ir=fread(istd,sizeof(stdint_t),3,fpar);
         ir+=fread(dstd,sizeof(double),1,fpar);
         error_root(ir!=4,1,"read_wflow_parms [iso1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
         {
            bswap_int(3,istd);
            bswap_double(1,dstd);
         }

         ie=0;
         ie|=(istd[0]!=(stdint_t)(flint));
         ie|=(istd[1]!=(stdint_t)(nstep));
         ie|=(istd[2]!=(stdint_t)(dnms));
         ie|=(dstd[0]!=eps);

         error_root(ie!=0,1,"read_wflow_parms [iso1.c]",
                    "Parameters do not match previous run");
      }
      else
      {
         istd[0]=(stdint_t)(flint);
         istd[1]=(stdint_t)(nstep);
         istd[2]=(stdint_t)(dnms);
         dstd[0]=eps;

         if (endian==BIG_ENDIAN)
         {
            bswap_int(3,istd);
            bswap_double(1,dstd);
         }

         iw=fwrite(istd,sizeof(stdint_t),3,fpar);
         iw+=fwrite(dstd,sizeof(double),1,fpar);
         error_root(iw!=4,1,"read_wflow_parms [iso1.c]",
                    "Incorrect write count");
      }
   }
}


static void read_infile(int argc,char *argv[])
{
   int ifile;
   FILE *fin=NULL,*fpar=NULL;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      noloc=find_opt(argc,argv,"-noloc");
      noexp=find_opt(argc,argv,"-noexp");
      rmold=find_opt(argc,argv,"-rmold");
      noms=find_opt(argc,argv,"-noms");
      scnfg=find_opt(argc,argv,"-c");
      append=find_opt(argc,argv,"-a");
      norng=find_opt(argc,argv,"-norng");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1))||(scnfg==(argc-1))||
                 ((append!=0)&&(scnfg==0)),1,"read_infile [iso1.c]",
                 "Syntax: iso1 -i <filename> [-noloc] [-noexp] "
                 "[-rmold] [-noms] [-c <filename> [-a [-norng]]]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [iso1.c]",
                 "Machine has unknown endianness");

      error_root((noexp)&&(noloc),1,"read_infile [iso1.c]",
		 "The concurrent use of -noloc and -noexp is not permitted");

      if (scnfg)
      {
         strncpy(cnfg,argv[scnfg+1],NAME_SIZE-1);
         cnfg[NAME_SIZE-1]='\0';
      }
      else
         cnfg[0]='\0';

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [iso1.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&noloc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&rmold,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noms,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&scnfg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&norng,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);
   }

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();

   if (my_rank==0)
   {
      if (append)
         fpar=fopen(par_file,"rb");
      else
         fpar=fopen(par_file,"wb");

      error_root(fpar==NULL,1,"read_infile [iso1.c]",
                 "Unable to open parameter file");
   }

   read_flds_bc_lat_parms(fpar);
   read_schedule(fpar);
   read_actions(fpar);
   read_integrator(fpar);
   read_solvers(fpar);
   read_wflow_parms(fpar);

   if (my_rank==0)
   {
      fclose(fin);
      fclose(fpar);

      if (append==0)
         copy_file(par_file,par_save);
   }
}


static void check_old_log(int ic,int *nl,int *icnfg)
{
   int ir,isv;
   int np[4],bp[4];
   FILE *fold;

   fold=fopen(log_file,"r");
   error_root(fold==NULL,1,"check_old_log [iso1.c]",
              "Unable to open log file");
   (*nl)=0;
   (*icnfg)=0;
   ir=1;
   isv=0;

   while (fgets(line,NAME_SIZE,fold)!=NULL)
   {
      if (strstr(line,"process grid")!=NULL)
      {
         ir&=(sscanf(line,"%dx%dx%dx%d process grid, %dx%dx%dx%d",
                     np,np+1,np+2,np+3,bp,bp+1,bp+2,bp+3)==8);

         ipgrd[0]=((np[0]!=NPROC0)||(np[1]!=NPROC1)||
                   (np[2]!=NPROC2)||(np[3]!=NPROC3));
         ipgrd[1]=((bp[0]!=NPROC0_BLK)||(bp[1]!=NPROC1_BLK)||
                   (bp[2]!=NPROC2_BLK)||(bp[3]!=NPROC3_BLK));
      }
      else if (strstr(line,"Trajectory no")!=NULL)
      {
         ir&=(sscanf(line,"Trajectory no %d",nl)==1);
         isv=0;
      }
      else if (strstr(line,"Configuration no")!=NULL)
      {
         ir&=(sscanf(line,"Configuration no %d",icnfg)==1);
         isv=1;
      }
   }

   fclose(fold);

   error_root(ir!=1,1,"check_old_log [iso1.c]","Incorrect read count");

   error_root(ic!=(*icnfg),1,"check_old_log [iso1.c]",
              "Continuation run:\n"
              "Initial configuration is not the last one of the previous run");

   error_root(isv==0,1,"check_old_log [iso1.c]",
              "Continuation run:\n"
              "The log file extends beyond the last configuration save");
}


static void check_old_dat(int nl)
{
   int nt;
   dat_t ndat;
   FILE *fdat;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [iso1.c]",
              "Unable to open data file");
   nt=0;

   while (read_dat(fdat,1,&ndat)==1)
      nt=ndat.nt;

   fclose(fdat);

   error_root(nt!=nl,1,"check_old_dat [iso1.c]",
              "Continuation run: Incomplete or too many data records");
}


static void check_old_ms3dat(int nl)
{
   int ic,ir,nt,pnt,dnt;
   FILE *fdat;

   fdat=fopen(ms3dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_ms3dat [iso1.c]",
              "Unable to open data file");

   check_file_head(fdat);

   nt=0;
   dnt=0;
   pnt=0;

   for (ic=0;;ic++)
   {
      ir=read_data3(fdat);

      if (ir==0)
      {
         error_root(ic==0,1,"check_old_ms3dat [iso1.c]",
                    "No data records found");
         break;
      }

      nt=data3.nt;

      if (ic==1)
      {
         dnt=nt-pnt;
         error_root(dnt<1,1,"check_old_ms3dat [iso1.c]",
                    "Incorrect trajectory separation");
      }
      else if (ic>1)
         error_root(nt!=(pnt+dnt),1,"check_old_ms3dat [iso1.c]",
                    "Trajectory sequence is not equally spaced");

      pnt=nt;
   }

   fclose(fdat);

   error_root((nt!=nl)||((ic>1)&&(dnt!=dtr_ms)),1,
              "check_old_ms3dat [iso1.c]","Last trajectory numbers "
              "or the trajectory separations do not match");
}


static void check_old_ms1dat(int nl)
{
   int ic,ir,nt,pnt,dnt;
   FILE *fdat;

   fdat=fopen(ms1dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_ms1dat [iso1.c]",
              "Unable to open data file");

   check_file_head(fdat);

   nt=0;
   dnt=0;
   pnt=0;

   for (ic=0;;ic++)
   {
      ir=read_data1(fdat);

      if (ir==0)
      {
         error_root(ic==0,1,"check_old_ms1dat [iso1.c]",
                    "No data records found");
         break;
      }

      nt=data1.nt;

      if (ic==1)
      {
         dnt=nt-pnt;
         error_root(dnt<1,1,"check_old_ms1dat [iso1.c]",
                    "Incorrect trajectory separation");
      }
      else if (ic>1)
         error_root(nt!=(pnt+dnt),1,"check_old_ms1dat [iso1.c]",
                    "Trajectory sequence is not equally spaced");

      pnt=nt;
   }

   fclose(fdat);

   error_root((nt!=nl)||((ic>1)&&(dnt!=dtr_ms)),1,
              "check_old_ms1dat [iso1.c]","Last trajectory numbers "
              "or the trajectory separations do not match");
}


static void check_files(int *nl,int *icnfg)
{
   int icmax,ic;
   FILE *fdat;

   ipgrd[0]=0;
   ipgrd[1]=0;

   if (my_rank==0)
   {
      if (noloc)
         error_root(cnfg[strlen(cnfg)-1]=='*',1,
                    "check_files [iso1.c]","Attempt to read an "
                    "imported configuration when -noloc is set");

      if (append)
      {
         error_root(strstr(cnfg,nbase)!=cnfg,1,"check_files [iso1.c]",
                    "Continuation run:\n"
                    "Run name does not match the previous one");
         error_root(sscanf(cnfg+strlen(nbase),"n%d",&ic)!=1,1,
                    "check_files [iso1.c]","Continuation run:\n"
                    "Unable to read configuration number from file name");

         check_old_log(ic,nl,icnfg);
         check_old_dat(*nl);
         if (noms==0)
         {
            check_old_ms3dat(*nl);
            check_old_ms1dat(*nl);
         }

         (*icnfg)+=1;
      }
      else
      {
         error_root(fopen(log_file,"r")!=NULL,1,
                    "check_files [iso1.c]",
                    "Attempt to overwrite old *.log file");

         error_root(fopen(dat_file,"rb")!=NULL,1,
                    "check_files [iso1.c]",
                    "Attempt to overwrite old *.dat file");
         
         if (noms==0)
         {
            error_root((fopen(ms3dat_file,"rb")!=NULL)||
                       (fopen(ms1dat_file,"rb")!=NULL),1,
                       "check_files [iso1.c]",
                       "Attempt to overwrite old *.ms*.dat file");

            fdat=fopen(ms1dat_file,"wb");
            error_root(fdat==NULL,1,"check_files [iso1.c]",
                       "Unable to open measurement data file (1)");
            write_file_head(fdat);
            fclose(fdat);

            fdat=fopen(ms3dat_file,"wb");
            error_root(fdat==NULL,1,"check_files [iso1.c]",
                       "Unable to open measurement data file (2)");
            write_file_head(fdat);
            fclose(fdat);
         }

         (*nl)=0;
         (*icnfg)=1;
      }

      icmax=(*icnfg)+(ntr-nth)/dtr_cnfg;

      if (noloc==0)
         error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,icmax,NPROC-1)>=
                    NAME_SIZE,1,"check_files [iso1.c]",
                    "loc_dir name is too long");

      if (noexp==0)
         error_root(name_size("%s/%sn%d",cnfg_dir,nbase,icmax)>=NAME_SIZE,1,
                    "check_files [iso1.c]","cnfg_dir name is too long");

      if (scnfg)
      {
         if (cnfg[strlen(cnfg)-1]=='*')
            error_root(name_size("%s/%s%d",loc_dir,cnfg,NPROC-1)>=NAME_SIZE,1,
                       "check_files [iso1.c]","loc_dir name is too long");
         else
            error_root(name_size("%s/%s",cnfg_dir,cnfg)>=NAME_SIZE,1,
                       "check_files [iso1.c]","cnfg_dir name is too long");
      }
   }

   MPI_Bcast(nl,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(icnfg,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void init_rng(int icnfg)
{
   int ic;

   if (append)
   {
      if (cnfg[strlen(cnfg)-1]!='*')
      {
         if (norng)
            start_ranlux(level,seed^(icnfg-1));
         else
         {
            ic=import_ranlux(rng_file);
            error_root(ic!=(icnfg-1),1,"init_rng [iso1.c]",
                       "Configuration number mismatch (*.rng file)");
         }
      }
   }
   else
      start_ranlux(level,seed);
}


static void init_gflds(void)
{
   int cnfg_type;
   char *p;

   if (scnfg)
   {
      if (cnfg[strlen(cnfg)-1]!='*')
      {
         sprintf(cnfg_file,"%s/%s",cnfg_dir,cnfg);
         cnfg_type=import_cnfg(cnfg_file);
         error(cnfg_type!=gauge(),1,"init_gflds [iso1.c]",
               "Starting configuration is not SU(3)xU(1)");
      }
      else
      {
         sprintf(line,"%s/%s",loc_dir,cnfg);
         p=line+strlen(line)-1;
         p[0]='\0';
         sprintf(cnfg_file,"%s_%d",line,my_rank);
         cnfg_type=read_cnfg(cnfg_file);
         error(cnfg_type!=gauge(),1,"init_gflds [iso1.c]",
               "Starting configuration is not SU(3)xU(1)");
      }
   }
   else
   {
      random_ud();
   }
}


static void store_gflds(su3_dble *usv,double *asv)
{
   cm3x3_assign(4*VOLUME,udfld(),usv);
   assign_dvec2dvec(4*VOLUME,adfld(),asv);
}


static void recall_gflds(su3_dble *usv,double *asv)
{
   cm3x3_assign(4*VOLUME,usv,udfld());
   set_flags(UPDATED_UD);
   assign_dvec2dvec(4*VOLUME,asv,adfld());
   set_flags(UPDATED_AD);
}


static void set_data(int nt)
{
   int in,dn,nn,x0;
   double eps;

   data1.nt=nt;
   data3.nt=nt;
   dn=file_head.dn;
   nn=file_head.nn;
   eps=file_head.eps;

   for (in=0;in<nn;in++)
   {
      Wact[in]=plaq_action_slices(data3.Wsl[in]);
      Yact[in]=ym_action_slices(data3.Ysl[in]);
      Qtop[in]=tcharge_slices(data3.Qsl[in]);
      if (bc_cstar()!=0)
      {
         for(x0=0;x0<N0;x0++)
         {
            data3.Wsl[in][x0]*=0.5;
            data3.Ysl[in][x0]*=0.5;
            data3.Qsl[in][x0]*=0.5;
         }
         Wact[in]*=0.5;
         Yact[in]*=0.5;
         Qtop[in]*=0.5;
      }

      if (flint==0)
         fwd_su3_euler(dn,eps);
      else if (flint==1)
         fwd_su3_rk2(dn,eps);
      else
         fwd_su3_rk3(dn,eps);
   }

   Wact[in]=plaq_action_slices(data3.Wsl[in]);
   Yact[in]=ym_action_slices(data3.Ysl[in]);
   Qtop[in]=tcharge_slices(data3.Qsl[in]);
   if (bc_cstar()!=0)
   {
      for(x0=0;x0<N0;x0++)
      {
         data3.Wsl[in][x0]*=0.5;
         data3.Ysl[in][x0]*=0.5;
         data3.Qsl[in][x0]*=0.5;
      }
      Wact[in]*=0.5;
      Yact[in]*=0.5;
      Qtop[in]*=0.5;
   }

   for (in=0;in<nn;in++)
   {
      Uact[in]=u1_plaq_action_slices(data1.Usl[in]);
      Mact[in]=mxw_action_slices(data1.Msl[in]);
      Flux[0][in]=u1fluxes_slices(0,data1.Fsl[0][in]);
      Flux[1][in]=u1fluxes_slices(1,data1.Fsl[1][in]);
      Flux[2][in]=u1fluxes_slices(2,data1.Fsl[2][in]);
      Flux[3][in]=u1fluxes_slices(3,data1.Fsl[3][in]);
      Flux[4][in]=u1fluxes_slices(4,data1.Fsl[4][in]);
      Flux[5][in]=u1fluxes_slices(5,data1.Fsl[5][in]);
      if (bc_cstar()!=0)
      {
         for(x0=0;x0<N0;x0++)
         {
            data1.Usl[in][x0]*=0.5;
            data1.Msl[in][x0]*=0.5;
         }
         Uact[in]*=0.5;
         Mact[in]*=0.5;
      }

      if (flint==0)
         fwd_u1_euler(dn,eps);
      else if (flint==1)
         fwd_u1_rk2(dn,eps);
      else
         fwd_u1_rk3(dn,eps);
   }

   Uact[in]=u1_plaq_action_slices(data1.Usl[in]);
   Mact[in]=mxw_action_slices(data1.Msl[in]);
   Flux[0][in]=u1fluxes_slices(0,data1.Fsl[0][in]);
   Flux[1][in]=u1fluxes_slices(1,data1.Fsl[1][in]);
   Flux[2][in]=u1fluxes_slices(2,data1.Fsl[2][in]);
   Flux[3][in]=u1fluxes_slices(3,data1.Fsl[3][in]);
   Flux[4][in]=u1fluxes_slices(4,data1.Fsl[4][in]);
   Flux[5][in]=u1fluxes_slices(5,data1.Fsl[5][in]);
   if (bc_cstar()!=0)
   {
      for(x0=0;x0<N0;x0++)
      {
         data1.Usl[in][x0]*=0.5;
         data1.Msl[in][x0]*=0.5;
      }
      Uact[in]*=0.5;
      Mact[in]*=0.5;
   }
}


static void print_info(int icnfg)
{
   int isap,idfl,n;
   long ip;

   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0L)
         remove("STARTUP_ERROR");

      if (append)
         flog=freopen(log_file,"a",stdout);
      else
         flog=freopen(log_file,"w",stdout);

      error_root(flog==NULL,1,"print_info [iso1.c]","Unable to open log file");

      if (append)
         printf("Continuation run, start from configuration %s\n\n",cnfg);
      else
      {
         printf("\n");
         printf("Simulation of QCD+QED with Wilson quarks\n");
         printf("----------------------------------------\n\n");

         if (scnfg)
            printf("New run, start from configuration %s\n\n",cnfg);
         else
            printf("New run, start from random configuration\n\n");

         printf("Using the HMC algorithm\n");
         printf("Program version %s\n",openQCD_RELEASE);
      }

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noloc)
         printf("The local disks are not used\n");
      if (noexp)
         printf("The generated configurations are not exported\n");
      if (rmold)
         printf("Old configurations are deleted\n");
      printf("\n");

      if ((ipgrd[0]!=0)&&(ipgrd[1]!=0))
         printf("Process grid and process block size changed:\n");
      else if (ipgrd[0]!=0)
         printf("Process grid changed:\n");
      else if (ipgrd[1]!=0)
         printf("Process block size changed:\n");

      if ((append==0)||(ipgrd[0]!=0)||(ipgrd[1]!=0))
      {
         printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
         printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
         printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
         printf("%dx%dx%dx%d process block size\n\n",
                NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);
      }

      if (append==0)
         print_flds_bc_lat_parms();

      printf("Random number generator:\n");

      if (append)
      {
         if (cnfg[strlen(cnfg)-1]!='*')
         {
            if (norng)
               printf("level = %d, seed = %d, effective seed = %d\n\n",
                      level,seed,seed^(icnfg-1));
            else
            {
               printf("State of ranlxs and ranlxd reset to the\n");
               printf("last exported state\n\n");
            }
         }
         else
         {
            printf("State of ranlxs and ranlxd read from\n");
            printf("initial field-configuration file\n\n");
         }
      }
      else
         printf("level = %d, seed = %d\n\n",level,seed);

      if (append)
      {
         printf("Trajectories:\n");
         printf("ntr = %d\n\n",ntr);
      }
      else
      {
         print_hmc_parms();

         printf("Trajectories:\n");
         printf("nth = %d, ntr = %d\n",nth,ntr);

         if (noms)
         {
            printf("dtr_log = %d, dtr_cnfg = %d\n\n",
                   dtr_log,dtr_cnfg);
            printf("Wilson flow observables are not measured\n\n");
         }
         else
         {
            printf("dtr_log = %d, dtr_ms = %d, dtr_cnfg = %d\n\n",
                   dtr_log,dtr_ms,dtr_cnfg);
            printf("Online measurement of Wilson flow observables\n\n");
         }

         print_action_parms();
         print_rat_parms();
         print_mdint_parms();
         print_force_parms2();
         print_solver_parms(&isap,&idfl);

         if (isap)
            print_sap_parms(0);

         if (idfl)
            print_dfl_parms(1);

         if (noms==0)
         {
            printf("Wilson flow:\n");
            if (flint==0)
               printf("Euler integrator\n");
            else if (flint==1)
               printf("2nd order RK integrator\n");
            else
               printf("3rd order RK integrator\n");
            n=fdigits(file_head.eps);
            printf("eps = %.*f\n",IMAX(n,1),file_head.eps);
            printf("nstep = %d\n",file_head.dn*file_head.nn);
            printf("dnms = %d\n\n",file_head.dn);
         }
      }

      fflush(flog);
   }
}


static void print_log(dat_t *ndat)
{
   if (my_rank==0)
   {
      printf("Trajectory no %d\n",(*ndat).nt);
      printf("dH = %+.1e, ",(*ndat).dH);
      printf("iac = %d\n",(*ndat).iac);
      printf("Average SU(3) plaquette = %.6f\n",(*ndat).avpl3);
      printf("Average U(1) plaquette = %.6f\n",(*ndat).avpl1);
      print_all_avgstat();
   }
}


static void save_dat(int n,double siac,double wtcyc,double wtall,dat_t *ndat)
{
   int iw;
   FILE *fdat;

   if (my_rank==0)
   {
      fdat=fopen(dat_file,"ab");
      error_root(fdat==NULL,1,"save_dat [iso1.c]",
                 "Unable to open data file");

      iw=write_dat(fdat,1,ndat);
      error_root(iw!=1,1,"save_dat [iso1.c]",
                 "Incorrect write count");
      fclose(fdat);

      printf("Acceptance rate = %1.2f\n",siac/(double)(n+1));
      printf("Time per trajectory = %.2e sec (average = %.2e sec)\n\n",
             wtcyc/(double)(dtr_log),wtall/(double)(n+1));
      fflush(flog);
   }
}


static void save_msdat(int n,double wtms,double wtmsall)
{
   int nms,in,dn,nn,din;
   double eps;
   FILE *fdat;

   if (my_rank==0)
   {
      fdat=fopen(ms3dat_file,"ab");
      error_root(fdat==NULL,1,"save_msdat [iso1.c]",
                 "Unable to open data file (1)");
      write_data3(fdat);
      fclose(fdat);

      fdat=fopen(ms1dat_file,"ab");
      error_root(fdat==NULL,1,"save_msdat [iso1.c]",
                 "Unable to open data file (2)");
      write_data1(fdat);
      fclose(fdat);

      nms=(n+1-nth)/dtr_ms+(nth>0);
      dn=file_head.dn;
      nn=file_head.nn;
      eps=file_head.eps;

      din=nn/10;
      if (din<1)
         din=1;

      printf("SU(3) observables:\n\n");

      for (in=0;in<=nn;in+=din)
         printf("n = %3d, t = %.2e, Wact = %.6e, Yact = %.6e, Q = % .2e\n",
                in*dn,eps*(double)(in*dn),Wact[in],Yact[in],Qtop[in]);

      printf("\n");
      printf("U(1) observable:\n\n");

      for (in=0;in<=nn;in+=din)
         printf("n = %3d, t = %.2e, Uact = %.6e, Mact = %.6e, F(01,02,03,23,31,12) = ( %+.2e, %+.2e, %+.2e, %+.2e, %+.2e, %+.2e )\n",
                in*dn,eps*(double)(in*dn),Uact[in],Mact[in],Flux[0][in],Flux[1][in],Flux[2][in],Flux[3][in],Flux[4][in],Flux[5][in]);

      printf("\n");
      printf("Configuration fully processed in %.2e sec ",wtms);
      printf("(average = %.2e sec)\n",wtmsall/(double)(nms));
      printf("Measured data saved\n\n");
      fflush(flog);
   }
}


static void save_cnfg(int icnfg)
{
   int ie;

   ie=check_bc(0.0)^0x1;
   error_root(ie!=0,1,"save_cnfg [iso1.c]",
              "Unexpected boundary values");

   if (noloc==0)
   {
      sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,icnfg,my_rank);
      write_cnfg(cnfg_file);
   }

   if (noexp==0)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      export_cnfg(cnfg_file);
   }

   if (my_rank==0)
   {
      if ((noloc==0)&&(noexp==0))
         printf("Configuration no %d saved on the local disks "
                "and exported\n\n",icnfg);
      else if (noloc==0)
         printf("Configuration no %d saved on the local disks\n\n",icnfg);
      else if (noexp==0)
         printf("Configuration no %d exported\n\n",icnfg);
   }
}


static void check_endflag(int *iend)
{
   FILE *fend;
   
   if (my_rank==0)
   {
      fend=fopen(end_file,"r");

      if (fend!=NULL)
      {
         fclose(fend);
         remove(end_file);
         (*iend)=1;
         printf("End flag set, run stopped\n\n");
      }
      else
         (*iend)=0;
   }

   MPI_Bcast(iend,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void remove_cnfg(int icnfg)
{
   if ((rmold)&&(icnfg>=1))
   {
      if (noloc==0)
      {
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,icnfg,my_rank);
         remove(cnfg_file);
      }

      if ((noexp==0)&&(my_rank==0))
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
         remove(cnfg_file);
      }
   }
}


int main(int argc,char *argv[])
{
   int nl,icnfg;
   int nwud,nwad,nws,nwsd,nwv,nwvd;
   int n,iend,iac,i;
   double *act0,*act1,w0[3],w1[3],npl,siac;
   double wt1,wt2,wtcyc,wtall,wtms,wtmsall;
   su3_dble **usv;
   double **asv;
   dat_t ndat;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   if (noms==0)
      alloc_data();
   check_files(&nl,&icnfg);
   geometry();

   hmc_wsize(&nwud,&nwad,&nws,&nwsd,&nwv,&nwvd);
   alloc_wud(nwud);
   alloc_wad(nwad);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);
   if ((noms==0)&&(flint))
   {
      alloc_wf3d(1);
      alloc_wf1d(1);
   }

   act0=malloc(2*(hmc.nact+1)*sizeof(*act0));
   act1=act0+hmc.nact+1;
   error(act0==NULL,1,"main [iso1.c]","Unable to allocate action arrays");

   print_info(icnfg);
   hmc_sanity_check();
   set_mdsteps();
   setup_counters();
   setup_chrono();
   init_rng(icnfg);
   init_gflds();

   if (bc_type()==0)
      npl=(double)(6*(N0-1)*N1)*(double)(N2*N3);
   else
      npl=(double)(6*N0*N1)*(double)(N2*N3);

   iend=0;
   siac=0.0;
   wtcyc=0.0;
   wtall=0.0;
   wtms=0.0;
   wtmsall=0.0;

   for (n=0;(iend==0)&&(n<ntr);n++)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      iac=run_hmc(act0,act1);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();

      siac+=(double)(iac);
      wtcyc+=(wt2-wt1);

      if (((ntr-n-1)%dtr_log)==0)
      {
         w0[0]=0.0;

         for (i=0;i<=hmc.nact;i++)
            w0[0]+=(act1[i]-act0[i]);

         w0[1]=plaq_wsum_dble(0)/npl;
         w0[2]=u1_plaq_wsum_dble(0)/npl;

         MPI_Reduce(w0,w1,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         MPI_Bcast(w1,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

         ndat.nt=nl+n+1;
         ndat.iac=iac;
         ndat.dH=w1[0];
         ndat.avpl3=w1[1];
         ndat.avpl1=w1[2];

         print_log(&ndat);
         wtall+=wtcyc;
         save_dat(n,siac,wtcyc,wtall,&ndat);
         wtcyc=0.0;

         if ((noms==0)&&((n+1)>=nth)&&(((ntr-n-1)%dtr_ms)==0))
         {
            MPI_Barrier(MPI_COMM_WORLD);
            wt1=MPI_Wtime();

            usv=reserve_wud(1);
            asv=reserve_wad(1);
            store_gflds(usv[0],asv[0]);
            set_data(nl+n+1);
            recall_gflds(usv[0],asv[0]);
            release_wud();
            release_wad();

            MPI_Barrier(MPI_COMM_WORLD);
            wt2=MPI_Wtime();

            wtms=wt2-wt1;
            wtmsall+=wtms;
            save_msdat(n,wtms,wtmsall);
         }
      }

      if (((n+1)>=nth)&&(((ntr-n-1)%dtr_cnfg)==0))
      {
         save_cnfg(icnfg);
         export_ranlux(icnfg,rng_file);
         check_endflag(&iend);

         if (my_rank==0)
         {
            fflush(flog);
            copy_file(log_file,log_save);
            copy_file(dat_file,dat_save);
            if (noms==0)
            {
               copy_file(ms3dat_file,ms3dat_save);
               copy_file(ms1dat_file,ms1dat_save);
            }
            copy_file(rng_file,rng_save);
         }

         remove_cnfg(icnfg-1);
         icnfg+=1;
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
