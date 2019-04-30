
/*******************************************************************************
*
* File ym1.c
*
* Copyright (C) 2010-2013, 2016 Martin Luescher
*               2017            Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* HMC simulation program for the SU(3) gauge theory.
*
* Syntax: ym1 -i <filename> [-noloc] [-noexp] [-rmold] [-noms]
*                           [-c <filename> [-a [-norng]]]
*
* For usage instructions see the file README.ym1.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "archive.h"
#include "update.h"
#include "version.h"
#include "global.h"
#include "lib/main_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int my_rank,noloc,noexp,rmold,noms,norng,unit;
static int scnfg,append,endian;
static int level,seed;
static int nth,ntr,dtr_log,dtr_ms,dtr_cnfg;
static int ipgrd[2];

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char ms3dat_file[NAME_SIZE],ms3dat_save[NAME_SIZE];
static char rng_file[NAME_SIZE],rng_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],end_file[NAME_SIZE];
static char nbase[NAME_SIZE],cnfg[NAME_SIZE];
static FILE *flog=NULL;

static hmc_parms_t hmc;


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
      if ((noexp==0)||((scnfg)&&(strchr(cnfg,'*')==NULL)))
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
   check_dir_root(log_dir);
   check_dir_root(dat_dir);

   if (noloc==0)
      check_dir(loc_dir);

   if (noexp==0)
      check_dir_root(cnfg_dir);

   error_root(name_size("%s/%s.log~",log_dir,nbase)>=NAME_SIZE,1,
              "setup_files [ym1.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms3.dat~",dat_dir,nbase)>=NAME_SIZE,1,
              "setup_files [ym1.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.par",dat_dir,nbase);
   sprintf(dat_file,"%s/%s.dat",dat_dir,nbase);
   sprintf(ms3dat_file,"%s/%s.ms3.dat",dat_dir,nbase);
   sprintf(rng_file,"%s/%s.rng",dat_dir,nbase);
   sprintf(end_file,"%s/%s.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
   sprintf(ms3dat_save,"%s~",ms3dat_file);
   sprintf(rng_save,"%s~",rng_file);
}


static void read_flds_bc_lat_parms(FILE *fpar)
{
   set_flds_parms(1,0);
   read_bc_parms();
   read_glat_parms();

   if (append)
      check_flds_bc_lat_parms(fpar);
   else
      write_flds_bc_lat_parms(fpar);
}


static void read_hmc_parms(FILE *fpar)
{
   int iact[1];
   double tau;

   if (my_rank==0)
   {
      find_section("HMC parameters");
      read_line("tau","%lf",&tau);
   }

   MPI_Bcast(&tau,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   iact[0]=0;
   hmc=set_hmc_parms(1,iact,0,0,NULL,1,tau,0);

   if (append)
      check_hmc_parms(fpar);
   else
      write_hmc_parms(fpar);
}


static void read_integrator(FILE *fpar)
{
   int nstep,imd,ifr[1];
   double lambda;

   if (my_rank==0)
   {
      find_section("MD integrator");
      read_line("integrator","%s",line);
      lambda=0.0;

      if (strcmp(line,"LPFR")==0)
         imd=(int)(LPFR);
      else if (strcmp(line,"OMF2")==0)
      {
         imd=(int)(OMF2);
         read_line("lambda","%lf",&lambda);
      }
      else if (strcmp(line,"OMF4")==0)
         imd=(int)(OMF4);
      else
         error_root(1,1,"read_integrator [ym1.c]","Unknown integrator");

      read_line("nstep","%d",&nstep);
   }

   MPI_Bcast(&imd,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&lambda,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&nstep,1,MPI_INT,0,MPI_COMM_WORLD);

   ifr[0]=0;

   if (imd==(int)(LPFR))
      set_mdint_parms(0,LPFR,lambda,nstep,1,ifr);
   else if (imd==(int)(OMF2))
      set_mdint_parms(0,OMF2,lambda,nstep,1,ifr);
   else if (imd==(int)(OMF4))
      set_mdint_parms(0,OMF4,lambda,nstep,1,ifr);

   set_action_parms(0,ACG_SU3,0,0,NULL,NULL,NULL);
   set_force_parms(0,FRG_SU3,0,0,NULL,NULL,NULL,NULL);

   if (append)
   {
      check_mdint_parms(fpar);
      check_action_parms(fpar);
      check_force_parms(fpar);
   }
   else
   {
      write_mdint_parms(fpar);
      write_action_parms(fpar);
      write_force_parms(fpar);
   }
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

      error_root((append!=0)&&(nth!=0),1,"read_schedule [ym1.c]",
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

      error_root(ie!=0,1,"read_schedule [ym1.c]",
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
         error_root(ir!=3,1,"read_schedule [ym1.c]",
                    "Incorrect read count");

         if (endian==BIG_ENDIAN)
            bswap_int(3,istd);

         ie=0;
         ie|=(istd[0]!=(stdint_t)(dtr_log));
         ie|=(istd[1]!=(stdint_t)(dtr_ms));
         ie|=(istd[2]!=(stdint_t)(dtr_cnfg));

         error_root(ie!=0,1,"read_schedule [ym1.c]",
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
         error_root(iw!=3,1,"read_schedule [ym1.c]",
                    "Incorrect write count");
      }
   }
}


static void read_observables(FILE *fpar)
{
   int itmp;
   
   if (my_rank==0)
   {
      if (append)
      {
         read_little_int(1,fpar,1,&itmp);
         error_root(itmp!=(stdint_t)(noms==0),1,"read_observables [ym1.c]",
                    "Attempt to mix measurement with other runs");
      }
      else
         write_little_int(1,fpar,1,noms==0);
   }
   
   if (noms==0)
   {
      read_wflow_parms();
      if (append)
         check_wflow_parms(fpar);
      else
         write_wflow_parms(fpar);
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
      unit=find_opt(argc,argv,"-unit");
      scnfg=find_opt(argc,argv,"-c");
      append=find_opt(argc,argv,"-a");
      norng=find_opt(argc,argv,"-norng");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1))||(scnfg==(argc-1))||
                 ((append!=0)&&(scnfg==0)),1,"read_infile [ym1.c]",
                 "Syntax: ym1 -i <filename> [-noloc] [-noexp] "
                 "[-rmold] [-noms] [-unit] [-c <filename> [-a [-norng]]]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ym1.c]",
                 "Machine has unknown endianness");

      error_root((noexp)&&(noloc),1,"read_infile [ym1.c]",
		 "The concurrent use of -noloc and -noexp is not permitted");

      error_root((unit)&&(scnfg),1,"read_infile [iso1.c]",
		 "The concurrent use of -unit and -c is not permitted");

      if (scnfg)
      {
         strncpy(cnfg,argv[scnfg+1],NAME_SIZE-1);
         cnfg[NAME_SIZE-1]='\0';
      }
      else
         cnfg[0]='\0';

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ym1.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&noloc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&rmold,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noms,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&unit,1,MPI_INT,0,MPI_COMM_WORLD);
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

      error_root(fpar==NULL,1,"read_infile [ym1.c]",
                 "Unable to open parameter file");
   }

   read_flds_bc_lat_parms(fpar);
   read_hmc_parms(fpar);
   read_schedule(fpar);
   read_integrator(fpar);
   read_observables(fpar);

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
   FILE *fold=NULL;

   fold=fopen(log_file,"r");
   error_root(fold==NULL,1,"check_old_log [ym1.c]",
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

   error_root(ir!=1,1,"check_old_log [ym1.c]","Incorrect read count");

   error_root(ic!=(*icnfg),1,"check_old_log [ym1.c]",
              "Continuation run:\n"
              "Initial configuration is not the last one of the previous run");

   error_root(isv==0,1,"check_old_log [ym1.c]",
              "Continuation run:\n"
              "The log file extends beyond the last configuration save");
}


static void check_old_dat(int nl)
{
   int nt;
   FILE *fdat=NULL;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [ym1.c]",
              "Unable to open data file");
   nt=0;

   while (read_dat(fdat,&nt)==1)
      continue;

   fclose(fdat);

   error_root(nt!=nl,1,"check_old_dat [ym1.c]",
              "Continuation run: Incomplete or too many data records");
}


static void check_old_ms3dat(int nl)
{
   int ic,ir,nt,pnt,dnt;
   FILE *fdat=NULL;

   fdat=fopen(ms3dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_ms3dat [ym1.c]",
              "Unable to open data file");

   check_ms3dat_head(fdat);

   nt=0;
   dnt=0;
   pnt=0;

   for (ic=0;;ic++)
   {
      ir=read_ms3dat(fdat,&nt);

      if (ir==0)
      {
         error_root(ic==0,1,"check_old_ms3dat [ym1.c]",
                    "No data records found");
         break;
      }

      if (ic==1)
      {
         dnt=nt-pnt;
         error_root(dnt<1,1,"check_old_ms3dat [ym1.c]",
                    "Incorrect trajectory separation");
      }
      else if (ic>1)
         error_root(nt!=(pnt+dnt),1,"check_old_ms3dat [ym1.c]",
                    "Trajectory sequence is not equally spaced");

      pnt=nt;
   }

   fclose(fdat);

   error_root((nt!=nl)||((ic>1)&&(dnt!=dtr_ms)),1,
              "check_old_ms3dat [ym1.c]","Last trajectory numbers "
              "or the trajectory separations do not match");
}


static void check_files(int *nl,int *icnfg)
{
   int icmax,ic;
   FILE *fdat=NULL;

   ipgrd[0]=0;
   ipgrd[1]=0;

   if (my_rank==0)
   {
      if (noloc)
         error_root(cnfg[strlen(cnfg)-1]=='*',1,
                    "check_files [ym1.c]","Attempt to read an "
                    "imported configuration when -noloc is set");

      if (append)
      {
         error_root(strstr(cnfg,nbase)!=cnfg,1,"check_files [ym1.c]",
                    "Continuation run:\n"
                    "Run name does not match the previous one");
         error_root(sscanf(cnfg+strlen(nbase),"n%d",&ic)!=1,1,
                    "check_files [ym1.c]","Continuation run:\n"
                    "Unable to read configuration number from file name");

         check_old_log(ic,nl,icnfg);
         check_old_dat(*nl);
         if (noms==0)
            check_old_ms3dat(*nl);

         (*icnfg)+=1;
      }
      else
      {
         error_root(fopen(log_file,"r")!=NULL,1,
                    "check_files [ym1.c]",
                    "Attempt to overwrite old *.log file");

         error_root(fopen(dat_file,"rb")!=NULL,1,
                    "check_files [ym1.c]",
                    "Attempt to overwrite old *.dat file");

         if (noms==0)
         {
            error_root(fopen(ms3dat_file,"rb")!=NULL,1,
                       "check_files [ym1.c]",
                       "Attempt to overwrite old *.ms*.dat file");

            fdat=fopen(ms3dat_file,"wb");
            error_root(fdat==NULL,1,"check_files [ym1.c]",
                       "Unable to open measurement data file");
            write_ms3dat_head(fdat);
            fclose(fdat);
         }

         (*nl)=0;
         (*icnfg)=1;
      }

      icmax=(*icnfg)+(ntr-nth)/dtr_cnfg;

      if (noloc==0)
         error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,icmax,NPROC-1)>=
                    NAME_SIZE,1,"check_files [ym1.c]",
                    "loc_dir name is too long");

      if (noexp==0)
         error_root(name_size("%s/%sn%d",cnfg_dir,nbase,icmax)>=NAME_SIZE,1,
                    "check_files [ym1.c]","cnfg_dir name is too long");

      if (scnfg)
      {
         if (cnfg[strlen(cnfg)-1]=='*')
            error_root(name_size("%s/%s%d",loc_dir,cnfg,NPROC-1)>=NAME_SIZE,1,
                       "check_files [ym1.c]","loc_dir name is too long");
         else
            error_root(name_size("%s/%s",cnfg_dir,cnfg)>=NAME_SIZE,1,
                       "check_files [ym1.c]","cnfg_dir name is too long");
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
            error_root(ic!=(icnfg-1),1,"init_rng [ym1.c]",
                       "Configuration number mismatch (*.rng file)");
         }
      }
   }
   else
      start_ranlux(level,seed);
}


static void init_ud(void)
{
   int cnfg_type;
   char *p;

   if (scnfg)
   {
      if (cnfg[strlen(cnfg)-1]!='*')
      {
         sprintf(cnfg_file,"%s/%s",cnfg_dir,cnfg);
         cnfg_type=import_cnfg(cnfg_file);
         error(cnfg_type!=gauge(),1,"init_ud [ym1.c]",
               "Starting configuration is not SU(3)");
      }
      else
      {
         sprintf(line,"%s/%s",loc_dir,cnfg);
         p=line+strlen(line)-1;
         p[0]='\0';
         sprintf(cnfg_file,"%s_%d",line,my_rank);
         cnfg_type=read_cnfg(cnfg_file);
         error(cnfg_type!=gauge(),1,"init_ud [ym1.c]",
               "Starting configuration is not SU(3)");
      }
   }
   else if (!unit)
      random_ud();
}


static void print_info(int icnfg)
{
   int n;
   long ip;
   mdint_parms_t mdp;

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

      error_root(flog==NULL,1,"print_info [ym1.c]","Unable to open log file");

      if (append)
         printf("Continuation run, start from configuration %s\n\n",cnfg);
      else
      {
         printf("\n");
         printf("Simulation of the SU(3) gauge theory\n");
         printf("------------------------------------\n\n");

         if (scnfg)
            printf("New run, start from configuration %s\n\n",cnfg);
         else if (unit)
            printf("New run, start from unit configuration\n\n");
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
         printf("Trajectories:\n");
         n=fdigits(hmc.tau);
         printf("tau = %.*f\n",IMAX(n,1),hmc.tau);

         mdp=mdint_parms(0);

         if (mdp.integrator==LPFR)
            printf("Leapfrog integrator\n");
         else if (mdp.integrator==OMF2)
            printf("2nd order OMF integrator with lambda = %.4e\n",mdp.lambda);
         else if (mdp.integrator==OMF4)
            printf("4th order OMF integrator\n");

         printf("Number of steps = %d\n\n",mdp.nstep);

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

            print_wflow_parms();
         }
      }

      fflush(flog);
   }
}


static void save_dat(int n,double siac,double wtcyc,double wtall)
{
   FILE *fdat=NULL;

   if (my_rank==0)
   {
      fdat=fopen(dat_file,"ab");
      error_root(fdat==NULL,1,"save_dat [ym1.c]",
                 "Unable to open data file");
      write_dat(fdat);
      fclose(fdat);

      printf("Acceptance rate = %1.2f\n",siac/(double)(n+1));
      printf("Time per trajectory = %.2e sec (average = %.2e sec)\n\n",
             wtcyc/(double)(dtr_log),wtall/(double)(n+1));
      fflush(flog);
   }
}


static void save_msdat(int n,double wtms,double wtmsall)
{
   int nms;
   FILE *fdat=NULL;

   if (my_rank==0)
   {
      fdat=fopen(ms3dat_file,"ab");
      error_root(fdat==NULL,1,"save_msdat [ym1.c]",
                 "Unable to open data file (1)");
      write_ms3dat(fdat);
      fclose(fdat);

      nms=(n+1-nth)/dtr_ms+(nth>0);

      print_ms3dat();

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
   error_root(ie!=0,1,"save_cnfg [ym1.c]",
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
   FILE *fend=NULL;

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
   int n,iend,iac;
   double act0[2],act1[2],w0[2],w1[2],npl,siac;
   double wt1,wt2,wtcyc,wtall,wtms,wtmsall;
   wflow_parms_t wfp;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   check_files(&nl,&icnfg);
   geometry();
   
   if (noms==0)
      wfp=wflow_parms();

   alloc_wud(1);
   if ((noms==0)&&(wfp.flint))
         alloc_wf3d(1);

   print_info(icnfg);
   set_mdsteps();
   init_rng(icnfg);
   init_ud();

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
         w0[0]=(act1[0]-act0[0])+(act1[1]-act0[1]);
         w0[1]=plaq_wsum_dble(0)/npl;

         MPI_Reduce(w0,w1,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         MPI_Bcast(w1,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

         set_dat(nl+n+1,iac,w1[0],w1[1],0.0);

         if (my_rank==0)
         {
            print_dat();
            print_all_avgstat();
         }
         wtall+=wtcyc;
         save_dat(n,siac,wtcyc,wtall);
         wtcyc=0.0;

         if ((noms==0)&&((n+1)>=nth)&&(((ntr-n-1)%dtr_ms)==0))
         {
            MPI_Barrier(MPI_COMM_WORLD);
            wt1=MPI_Wtime();

            set_ms3dat(nl+n+1);

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
               copy_file(ms3dat_file,ms3dat_save);
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
