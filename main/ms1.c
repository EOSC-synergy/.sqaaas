
/*******************************************************************************
*
* File ms1.c
*
* Copyright (C) 2012-2014, 2016 Martin Luescher
*               2017, 2019      Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Stochastic estimation of reweighting factors.
*
* Syntax: ms1 -i <input file> [-noexp] [-a [-norng]]
*
* For usage instructions see the file README.ms1.
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
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "archive.h"
#include "dfl.h"
#include "update.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static struct
{
   int nrw;
   int *nfct,*nsrc;
} file_head;

static struct
{
   int nc;
   double ***sqn,***lnr;
} data;

static int my_rank,noexp,append,norng,endian;
static int first,last,step,level,seed;
static int ipgrd[2],**rwstat=NULL,*rlxs_state=NULL,*rlxd_state=NULL;

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char rng_file[NAME_SIZE],rng_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],nbase[NAME_SIZE];
static FILE *flog=NULL;


static void alloc_data(void)
{
   int nrw,*nfct,*nsrc;
   int i,irw,ifct,n1,n2,n3;
   double ***ppp,**pp,*p;

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;
   n1=nrw;
   n2=0;
   n3=0;

   for (irw=0;irw<nrw;irw++)
   {
      n2+=nfct[irw];
      n3+=(nfct[irw]*nsrc[irw]);
   }

   ppp=malloc(2*n1*sizeof(*ppp));
   pp=malloc(2*n2*sizeof(*pp));
   p=malloc(2*n3*sizeof(*p));
   error((ppp==NULL)||(pp==NULL)||(p==NULL),1,"alloc_data [ms1.c]",
         "Unable to allocate data arrays");

   data.sqn=ppp;
   data.lnr=ppp+nrw;

   for (i=0;i<2;i++)
   {
      for (irw=0;irw<nrw;irw++)
      {
         (*ppp)=pp;
         ppp+=1;

         for (ifct=0;ifct<nfct[irw];ifct++)
         {
            (*pp)=p;
            pp+=1;
            p+=nsrc[irw];
         }
      }
   }
}


static void write_file_head(FILE *fpar)
{
   int nrw,iw;

   nrw=file_head.nrw;
   iw=write_little_int(0,fpar,1,nrw);
   iw+=write_little_intarray(0,fpar,nrw,file_head.nfct);
   iw+=write_little_intarray(0,fpar,nrw,file_head.nsrc);

   error_root(iw!=(1+2*file_head.nrw),1,"write_file_head [ms1.c]",
              "Incorrect write count");
}


static void check_file_head(FILE *fpar)
{
   int nrw;

   nrw=file_head.nrw;
   check_little_int("check_file_head [ms1.c]",fpar,1,nrw);
   check_little_intarray("check_file_head [ms1.c]",fpar,nrw,file_head.nfct);
   check_little_intarray("check_file_head [ms1.c]",fpar,nrw,file_head.nsrc);
}


static void write_data(FILE *fdat)
{
   int iw,n;
   int nrw,*nfct,*nsrc,irw,ifct;

   iw=write_little_int(0,fdat,1,data.nc);

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;
   n=0;

   for (irw=0;irw<nrw;irw++)
   {
      for (ifct=0;ifct<nfct[irw];ifct++)
      {
         iw+=write_little_dblearray(0,fdat,nsrc[irw],data.sqn[irw][ifct]);
         iw+=write_little_dblearray(0,fdat,nsrc[irw],data.lnr[irw][ifct]);
         n+=nsrc[irw];
      }
   }

   error_root(iw!=(1+2*n),1,"write_data [ms1.c]",
              "Incorrect write count");
}


static int read_data(FILE *fdat,int *nt)
{
   int ir,n;
   int nrw,*nfct,*nsrc,irw,ifct;

   ir=read_little_int(0,fdat,1,&(data.nc));
   if (ir!=1) return 0;

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;
   n=0;

   for (irw=0;irw<nrw;irw++)
   {
      for (ifct=0;ifct<nfct[irw];ifct++)
      {
         ir+=read_little_dblearray(0,fdat,nsrc[irw],data.sqn[irw][ifct]);
         ir+=read_little_dblearray(0,fdat,nsrc[irw],data.lnr[irw][ifct]);
         n+=nsrc[irw];
      }
   }

   error_root(ir!=(1+2*n),1,"read_data [ms1.c]",
              "Read error or incomplete data record");

   (*nt)=data.nc;

   return 1;
}


static void read_dirs(void)
{
   int nrw,*nfct;

   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);
      read_line("dat_dir","%s",dat_dir);

      if (noexp)
      {
         read_line("loc_dir","%s",loc_dir);
         cnfg_dir[0]='\0';
      }
      else
      {
         read_line("cnfg_dir","%s",cnfg_dir);
         loc_dir[0]='\0';
      }

      find_section("Configurations");
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);
      read_line("nrw","%d",&nrw);

      error_root((last<first)||(step<1)||(((last-first)%step)!=0),1,
                 "read_dirs [ms1.c]","Improper configuration range");
      error_root(nrw<1,1,"read_dirs [ms1.c]",
                 "The number nrw or reweighting factors must be at least 1");
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nrw,1,MPI_INT,0,MPI_COMM_WORLD);

   nfct=malloc(2*nrw*sizeof(*nfct));
   error(nfct==NULL,1,"read_dirs [ms1.c]",
         "Unable to allocate data array");
   file_head.nrw=nrw;
   file_head.nfct=nfct;
   file_head.nsrc=nfct+nrw;
}


static void setup_files(void)
{
   if (noexp)
      error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,last,NPROC-1)>=NAME_SIZE,
                 1,"setup_files [ms1.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [ms1.c]","cnfg_dir name is too long");

   check_dir_root(log_dir);
   check_dir_root(dat_dir);
   error_root(name_size("%s/%s.ms1.log~",log_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms1.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms1.dat~",dat_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms1.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.ms1.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.ms1.par",dat_dir,nbase);
   sprintf(dat_file,"%s/%s.ms1.dat",dat_dir,nbase);
   sprintf(rng_file,"%s/%s.ms1.rng",dat_dir,nbase);
   sprintf(end_file,"%s/%s.ms1.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
   sprintf(rng_save,"%s~",rng_file);
}


static void read_flds_bc_lat_parms(FILE *fpar)
{
   int gg,nfl;

   if (my_rank==0)
   {
      find_section("Field parameters");
      read_line("gauge","%s",line);

      gg=0;
      if (strcmp(line,"SU(3)")==0)
         gg=1;
      else if (strcmp(line,"U(1)")==0)
         gg=2;
      else if (strcmp(line,"SU(3)xU(1)")==0)
         gg=3;
      else
         error_root(1,1,"read_flds_bc_lat_parms [ms1.c]",
                    "Unknown gauge group %s",line);
      read_line("nfl","%d",&nfl);
   }
   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nfl,1,MPI_INT,0,MPI_COMM_WORLD);

   set_flds_parms(gg,nfl);
   read_bc_parms();
   read_qlat_parms();

   if (append)
      check_flds_bc_lat_parms(fpar);
   else
      write_flds_bc_lat_parms(fpar);
}


static void read_rw_factors(FILE *fpar)
{
   int nrw,*nfct,*nsrc,irw,np,k;
   rw_parms_t rwp;
   rat_parms_t rp[2];

   nrw=file_head.nrw;
   nfct=file_head.nfct;
   nsrc=file_head.nsrc;

   for (irw=0;irw<nrw;irw++)
   {
      read_rw_parms(irw);
      rwp=rw_parms(irw);
      nsrc[irw]=rwp.nsrc;

      if (rwp.rwfact==RWRAT)
      {
         nfct[irw]=1;
         rp[0]=rat_parms(rwp.irp[0]);

         if (rp[0].degree==0)
         {
            read_rat_parms(rwp.irp[0]);
            rp[0]=rat_parms(rwp.irp[0]);
         }

         np=rp[0].degree;
         for (k=0;k<rwp.nfct;k++)
            np-=rwp.np[k];

         error_root(np!=0,1,"read_rw_factors [ms1.c]",
                    "Invalid pole decomposition in reweighting factor %d",irw);

      }
      else if (rwp.rwfact==RWRTM)
      {
         nfct[irw]=rwp.nfct;
         rp[0]=rat_parms(rwp.irp[0]);
         rp[1]=rat_parms(rwp.irp[1]);

         if (rp[0].degree==0)
         {
            read_rat_parms(rwp.irp[0]);
            rp[0]=rat_parms(rwp.irp[0]);
         }
         if (rp[1].degree==0)
         {
            read_rat_parms(rwp.irp[1]);
            rp[1]=rat_parms(rwp.irp[1]);
         }

         error_root(rp[0].degree!=rp[1].degree,1,"read_rw_factors [ms1.c]",
                    "Rational approximations in RWRTM reweighting factor %d "
                    "must have the same order",irw);

         error_root((rp[0].power[0]!=rp[1].power[0])||
                    (rp[0].power[1]!=rp[1].power[1]),1,
                    "read_rw_factors [ms1.c]",
                    "Rational approximations in RWRTM reweighting factor %d "
                    "must have the same power",irw);

         error_root(rp[1].mu0!=0.0,1,"read_rw_factors [ms1.c]",
                    "Rational approximations in RWRTM reweighting factor %d "
                    "must have zero twisted mass",irw);

         np=rp[0].degree;
         for (k=0;k<rwp.nfct;k++)
            np-=rwp.np[k];

         error_root(np!=0,1,"read_rw_factors [ms1.c]",
                    "Invalid pole decomposition in reweighting factor %d (%d,%d,%d)",irw,np,rp[0].degree,rwp.nfct);
      }
      else
         nfct[irw]=rwp.nfct;
   }

   if (append)
   {
      check_rw_parms(fpar);
      check_rat_parms(fpar);
   }
   else
   {
      write_rw_parms(fpar);
      write_rat_parms(fpar);
   }
}


static void read_solvers(FILE *fpar)
{
   int nrw,nfct,irw,ifct,isp;
   int isap,idfl;
   rw_parms_t rwp;
   solver_parms_t sp;

   nrw=file_head.nrw;
   isap=0;
   idfl=0;

   for (irw=0;irw<nrw;irw++)
   {
      rwp=rw_parms(irw);
      nfct=rwp.nfct;

      for (ifct=0;ifct<nfct;ifct++)
      {
         isp=rwp.isp[ifct];
         sp=solver_parms(isp);

         if (sp.solver==SOLVERS)
         {
            read_solver_parms(isp);
            sp=solver_parms(isp);

            if (sp.solver==SAP_GCR)
               isap=1;
            else if (sp.solver==DFL_SAP_GCR)
            {
               isap=1;
               idfl=1;
               
               if (dfl_gen_parms(sp.idfl).status!=DFL_DEF)
                  read_dfl_parms(sp.idfl);
            }
         }
      }
   }

   if (append)
      check_solver_parms(fpar);
   else
      write_solver_parms(fpar);

   if (isap)
   {
      read_sap_parms();
      if (append)
         check_sap_parms(fpar);
      else
         write_sap_parms(fpar);
   }

   if (idfl)
   {
      read_dfl_parms(-1);
      if (append)
         check_dfl_parms(fpar);
      else
         write_dfl_parms(fpar);
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
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [ms1.c]",
                 "Syntax: ms1 -i <input file> [-noexp] [-a [-norng]]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ms1.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");
      append=find_opt(argc,argv,"-a");
      norng=find_opt(argc,argv,"-norng");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ms1.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&norng,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();

   if (my_rank==0)
   {
      if (append)
         fpar=fopen(par_file,"rb");
      else
         fpar=fopen(par_file,"wb");

      error_root(fpar==NULL,1,"read_infile [ms1.c]",
                 "Unable to open parameter file");
   }

   if (my_rank==0)
   {
      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);
   }

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

   read_flds_bc_lat_parms(fpar);
   read_rw_factors(fpar);
   read_solvers(fpar);

   if (my_rank==0)
   {
      fclose(fin);
      fclose(fpar);

      if (append==0)
         copy_file(par_file,par_save);
   }
}


static void check_old_log(int *fst,int *lst,int *stp)
{
   int ie,ic,isv;
   int fc,lc,dc,pc;
   int np[4],bp[4];
   FILE *fold=NULL;

   fold=fopen(log_file,"r");
   error_root(fold==NULL,1,"check_old_log [ms1.c]",
              "Unable to open log file");

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;
   isv=0;

   while (fgets(line,NAME_SIZE,fold)!=NULL)
   {
      if (strstr(line,"process grid")!=NULL)
      {
         if (sscanf(line,"%dx%dx%dx%d process grid, %dx%dx%dx%d",
                    np,np+1,np+2,np+3,bp,bp+1,bp+2,bp+3)==8)
         {
            ipgrd[0]=((np[0]!=NPROC0)||(np[1]!=NPROC1)||
                      (np[2]!=NPROC2)||(np[3]!=NPROC3));
            ipgrd[1]=((bp[0]!=NPROC0_BLK)||(bp[1]!=NPROC1_BLK)||
                      (bp[2]!=NPROC2_BLK)||(bp[3]!=NPROC3_BLK));
         }
         else
            ie|=0x1;
      }
      else if (strstr(line,"fully processed")!=NULL)
      {
         pc=lc;

         if (sscanf(line,"Configuration no %d",&lc)==1)
         {
            ic+=1;
            isv=1;
         }
         else
            ie|=0x1;

         if (ic==1)
            fc=lc;
         else if (ic==2)
            dc=lc-fc;
         else if ((ic>2)&&(lc!=(pc+dc)))
            ie|=0x2;
      }
      else if (strstr(line,"Configuration no")!=NULL)
         isv=0;
   }

   fclose(fold);

   error_root((ie&0x1)!=0x0,1,"check_old_log [ms1.c]",
              "Incorrect read count");
   error_root((ie&0x2)!=0x0,1,"check_old_log [ms1.c]",
              "Configuration numbers are not equally spaced");
   error_root(isv==0,1,"check_old_log [ms1.c]",
              "Log file extends beyond the last configuration save");

   (*fst)=fc;
   (*lst)=lc;
   (*stp)=dc;
}


static void check_old_dat(int fst,int lst,int stp)
{
   int ie,ic;
   int fc,lc,dc,pc;
   FILE *fdat=NULL;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [ms1.c]",
              "Unable to open data file");

   check_file_head(fdat);

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;

   while (read_data(fdat,&lc)==1)
   {
      ic+=1;

      if (ic==1)
         fc=lc;
      else if (ic==2)
         dc=lc-fc;
      else if ((ic>2)&&(lc!=(pc+dc)))
         ie|=0x1;

      pc=lc;
   }

   fclose(fdat);

   error_root(ic==0,1,"check_old_dat [ms1.c]",
              "No data records found");
   error_root((ie&0x1)!=0x0,1,"check_old_dat [ms1.c]",
              "Configuration numbers are not equally spaced");
   error_root((fst!=fc)||(lst!=lc)||(stp!=dc),1,"check_old_dat [ms1.c]",
              "Configuration range is not as reported in the log file");
}


static void check_files(void)
{
   int fst,lst,stp;
   FILE *fdat=NULL;

   ipgrd[0]=0;
   ipgrd[1]=0;

   if (my_rank==0)
   {
      if (append)
      {
         check_old_log(&fst,&lst,&stp);
         check_old_dat(fst,lst,stp);

         error_root((fst!=lst)&&(stp!=step),1,"check_files [ms1.c]",
                    "Continuation run:\n"
                    "Previous run had a different configuration separation");
         error_root(first!=lst+step,1,"check_files [ms1.c]",
                    "Continuation run:\n"
                    "Configuration range does not continue the previous one");
      }
      else
      {
         error_root(fopen(log_file,"r")!=NULL,1,
                    "check_files [ms1.c]",
                    "Attempt to overwrite old *.log file");

         error_root(fopen(dat_file,"rb")!=NULL,1,
                    "check_files [ms1.c]",
                    "Attempt to overwrite old *.ms*.dat file");

         fdat=fopen(dat_file,"wb");
         error_root(fdat==NULL,1,"check_files [ms1.c]",
                    "Unable to open data file");
         write_file_head(fdat);
         fclose(fdat);
      }
   }
}


static void print_info(void)
{
   int isap,idfl;
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

      error_root(flog==NULL,1,"print_info [ms1.c]","Unable to open log file");
      printf("\n");

      if (append)
         printf("Continuation run\n\n");
      else
      {
         printf("Measurement of reweighting factors\n");
         printf("----------------------------------\n\n");
      }

      printf("Program version %s\n",openQCD_RELEASE);

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noexp)
         printf("Configurations are read in imported file format\n\n");
      else
         printf("Configurations are read in exported file format\n\n");

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

      if (append)
      {
         printf("Random number generator:\n");

         if (norng)
            printf("level = %d, seed = %d, effective seed = %d\n\n",
                   level,seed,seed^(first-step));
         else
         {
            printf("State of ranlxs and ranlxd reset to the\n");
            printf("last exported state\n\n");
         }
      }
      else
      {
         printf("Random number generator:\n");
         printf("level = %d, seed = %d\n\n",level,seed);

         print_bc_parms();
         print_flds_parms();
         print_lat_parms();
         print_rw_parms();
         print_rat_parms();
         print_solver_parms(&isap,&idfl);

         if (isap)
            print_sap_parms(0);

         if (idfl)
            print_dfl_parms(0);
      }

      printf("Configurations no %d -> %d in steps of %d\n\n",
             first,last,step);
      fflush(flog);
   }
}


static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpp;

   dp=dfl_parms();
   dpp=dfl_pro_parms();

   MAX(*nws,dp.Ns+2);
   MAX(*nwv,2*dpp.nkv+2);
   MAX(*nwvd,4);
}


static void solver_wsize(int isp,int nsd,int np,
                         int *nws,int *nwsd,int *nwv,int *nwvd)
{
   solver_parms_t sp;

   sp=solver_parms(isp);

   if (sp.solver==CGNE)
   {
      MAX(*nws,5);
      MAX(*nwsd,nsd+3);
   }
   else if (sp.solver==MSCG)
   {
      MAX(*nws,5);
      if (np>1)
      {
         MAX(*nwsd,nsd+np+3);
      }
      else
      {
         MAX(*nwsd,nsd+5);
      }
   }
   else if (sp.solver==SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+1);
      MAX(*nwsd,nsd+2);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      MAX(*nws,2*sp.nkv+2);
      MAX(*nwsd,nsd+4);
      dfl_wsize(nws,nwv,nwvd);
   }
}


static void reweight_wsize(int *nws,int *nwsd,int *nwv,int *nwvd)
{
   int nrw,nfct;
   int irw,ifct,nsd;
   int *np,*isp;
   rw_parms_t rwp;
   solver_parms_t sp;

   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;
   nrw=file_head.nrw;

   for (irw=0;irw<nrw;irw++)
   {
      rwp=rw_parms(irw);
      nfct=rwp.nfct;
      np=rwp.np;
      isp=rwp.isp;

      for (ifct=0;ifct<nfct;ifct++)
      {
         if (rwp.rwfact==RWRAT)
         {
            sp=solver_parms(isp[ifct]);

            if (sp.solver==MSCG)
               nsd=3+np[ifct];
            else
               nsd=5;

            solver_wsize(isp[ifct],nsd,np[ifct],nws,nwsd,nwv,nwvd);
         }
         else if (rwp.rwfact==RWRTM)
         {
            sp=solver_parms(isp[ifct]);

            if (sp.solver==MSCG)
               nsd=3+np[ifct];
            else
               nsd=4;

            solver_wsize(isp[ifct],nsd,np[ifct],nws,nwsd,nwv,nwvd);
         }
         else
         {
            nsd=2;
            solver_wsize(isp[ifct],nsd,0,nws,nwsd,nwv,nwvd);
         }
      }
   }
}


static void alloc_rwstat(void)
{
   int nrw,mfct;
   int irw,ifct;
   int **pp,*p;
   rw_parms_t rwp;

   nrw=file_head.nrw;
   mfct=0;

   for (irw=0;irw<nrw;irw++)
   {
      rwp=rw_parms(irw);

      if (mfct<rwp.nfct)
         mfct=rwp.nfct;
   }

   mfct*=2;
   pp=malloc(mfct*sizeof(*pp));
   p=malloc(3*mfct*sizeof(*p));
   error((pp==NULL)||(p==NULL),1,"alloc_rwstat [ms1.c]",
         "Unable to allocate status array");
   rwstat=pp;

   for (ifct=0;ifct<mfct;ifct++)
   {
      (*pp)=p;
      pp+=1;
      p+=3;
   }
}


static void print_rwstat(int irw)
{
   int nfct,nsrc,*isp;
   int nrs,ifct,j;
   rw_parms_t rwp;
   solver_parms_t sp;

   rwp=rw_parms(irw);
   nfct=rwp.nfct;
   nsrc=rwp.nsrc;
   isp=rwp.isp;
   nrs=0;

   for (ifct=0;ifct<nfct;ifct++)
   {
      if (ifct==0)
         printf("RWF %d: status = ",irw);
      else
         printf(";");

      sp=solver_parms(isp[ifct]);

      if (sp.solver==DFL_SAP_GCR)
      {
         for (j=0;j<2;j++)
            rwstat[ifct][j]=(rwstat[ifct][j]+(nsrc/2))/nsrc;

         printf("%d,%d",rwstat[ifct][0],rwstat[ifct][1]);
         nrs+=rwstat[ifct][2];
      }
      else
      {
         rwstat[ifct][0]=(rwstat[ifct][0]+(nsrc/2))/nsrc;
         printf("%d",rwstat[ifct][0]);
      }
   }

   if (nrs)
      printf(" (no of subspace regenerations = %d)\n",nrs);
   else
      printf("\n");
}


static void set_data(int nc)
{
   int nrw,nfct,nsrc;
   int irp1,irp2,np1,np2;
   int irw,ifct,isrc,isp,j;
   double mu1,mu2,*sqn,*lnr;
   rw_parms_t rwp;
   dirac_parms_t dp;

   nrw=file_head.nrw;
   data.nc=nc;

   for (irw=0;irw<nrw;irw++)
   {
      rwp=rw_parms(irw);
      dp=qlat_parms(rwp.ifl);
      set_dirac_parms1(&dp);
      nsrc=rwp.nsrc;
      nfct=rwp.nfct;

      for (ifct=0;ifct<(2*nfct);ifct++)
      {
         for (j=0;j<3;j++)
            rwstat[ifct][j]=0;
      }

      if (rwp.rwfact==RWRAT)
      {
         sqn=data.sqn[irw][0];
         lnr=data.lnr[irw][0];

         for (isrc=0;isrc<nsrc;isrc++)
         {
            lnr[isrc]=rwrat(rwp.irp[0],nfct,rwp.np,rwp.isp,sqn+isrc,rwstat+nfct);

            for (ifct=0;ifct<nfct;ifct++)
            {
               for (j=0;j<3;j++)
                  rwstat[ifct][j]+=rwstat[nfct+ifct][j];
            }
         }
      }
      else if (rwp.rwfact==RWRTM)
      {
         np1=0;
         for (ifct=0;ifct<nfct;ifct++)
         {
            sqn=data.sqn[irw][ifct];
            lnr=data.lnr[irw][ifct];
      
            irp1=rwp.irp[0];
            irp2=rwp.irp[1];
            isp=rwp.isp[ifct];
      
            np2=np1+rwp.np[ifct]-1;
            for (isrc=0;isrc<nsrc;isrc++)
            {
               lnr[isrc]=rwrtm(irp1,irp2,np1,np2,isp,sqn+isrc,rwstat[nfct]);

               for (j=0;j<3;j++)
                  rwstat[ifct][j]+=rwstat[nfct][j];
            }
            np1=np2+1;
         }
      }
      else
      {
         for (ifct=0;ifct<nfct;ifct++)
         {
            sqn=data.sqn[irw][ifct];
            lnr=data.lnr[irw][ifct];

            if (ifct>0)
               mu1=rwp.mu[ifct-1];
            else
               mu1=0.0;

            mu2=rwp.mu[ifct];
            isp=rwp.isp[ifct];

            for (isrc=0;isrc<nsrc;isrc++)
            {
               if (rwp.rwfact==RWTM1)
                  lnr[isrc]=rwtm1(mu1,mu2,isp,sqn+isrc,rwstat[nfct]);
               else if (rwp.rwfact==RWTM1_EO)
                  lnr[isrc]=rwtm1eo(mu1,mu2,isp,sqn+isrc,rwstat[nfct]);
               else if (rwp.rwfact==RWTM2)
                  lnr[isrc]=rwtm2(mu1,mu2,isp,sqn+isrc,rwstat[nfct]);
               else if (rwp.rwfact==RWTM2_EO)
                  lnr[isrc]=rwtm2eo(mu1,mu2,isp,sqn+isrc,rwstat[nfct]);

               for (j=0;j<3;j++)
                  rwstat[ifct][j]+=rwstat[nfct][j];
            }
         }
      }

      if (my_rank==0)
      {
         print_rwstat(irw);

         if (rwp.rwfact==RWRAT)
            nfct=1;
         else
            nfct=rwp.nfct;

         for (ifct=0;ifct<nfct;ifct++)
         {
            lnr=data.lnr[irw][ifct];

            if (nfct==1)
               printf("RWF %d: -ln(r) = %.4e",irw,lnr[0]);
            else
               printf("RWF %d, factor %d: -ln(r) = %.4e",irw,ifct,lnr[0]);

            if (nsrc<=4)
            {
               for (isrc=1;isrc<nsrc;isrc++)
                  printf(",%.4e",lnr[isrc]);
            }
            else
            {
               printf(",%.4e,...",lnr[1]);

               for (isrc=(nsrc-2);isrc<nsrc;isrc++)
                  printf(",%.4e",lnr[isrc]);
            }

            printf("\n");
         }
      }
   }
}


static void init_rng(void)
{
   int ic;

   if (append)
   {
      if (norng)
         start_ranlux(level,seed^(first-step));
      else
      {
         ic=import_ranlux(rng_file);
         error_root(ic!=(first-step),1,"init_rng [ms1.c]",
                    "Configuration number mismatch (*.rng file)");
      }
   }
   else
      start_ranlux(level,seed);
}


static void save_ranlux(void)
{
   int nlxs,nlxd;

   if (rlxs_state==NULL)
   {
      nlxs=rlxs_size();
      nlxd=rlxd_size();

      rlxs_state=malloc((nlxs+nlxd)*sizeof(int));
      rlxd_state=rlxs_state+nlxs;

      error(rlxs_state==NULL,1,"save_ranlux [ms1.c]",
            "Unable to allocate state arrays");
   }

   rlxs_get(rlxs_state);
   rlxd_get(rlxd_state);
}


static void restore_ranlux(void)
{
   rlxs_reset(rlxs_state);
   rlxd_reset(rlxd_state);
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


int main(int argc,char *argv[])
{
   int nc,iend,status,idfl;
   int nws,nwsd,nwv,nwvd;
   int cnfg_type;
   double wt1,wt2,wtavg;
   dfl_parms_t dfl;
   dflst_t dfl_status;
   FILE *fdat=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   alloc_data();
   check_files();
   print_info();
   dfl=dfl_parms();

   geometry();
   init_rng();

   reweight_wsize(&nws,&nwsd,&nwv,&nwvd);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);
   alloc_rwstat();

   iend=0;
   wtavg=0.0;

   for (nc=first;(iend==0)&&(nc<=last);nc+=step)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      if (my_rank==0)
         printf("Configuration no %d\n",nc);

      if (noexp)
      {
         save_ranlux();
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,nc,my_rank);
         read_cnfg(cnfg_file);
         restore_ranlux();
      }
      else
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,nc);
         cnfg_type=import_cnfg(cnfg_file);
         error(cnfg_type!=gauge(),1,"main [ms1.c]",
               "Imported configuration is not of the correct type");
      }

      if (dfl.Ns)
      {
         idfl=0;
         while(1)
         {
            dfl_status=dfl_gen_parms(idfl).status;
            if(dfl_status==DFL_OUTOFRANGE) break;
            if(dfl_status==DFL_DEF)
            {
               dfl_modes(idfl,&status);
               error_root(status<0,1,"main [ms1.c]",
                          "Generation of deflation subspace %d failed (status = %d)",
                          idfl,status);
            }
            idfl++;
         }
      }

      set_data(nc);

      if (my_rank==0)
      {
         fdat=fopen(dat_file,"ab");
         error_root(fdat==NULL,1,"main [ms1.c]",
                    "Unable to open dat file");
         write_data(fdat);
         fclose(fdat);
      }

      export_ranlux(nc,rng_file);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtavg+=(wt2-wt1);

      if (my_rank==0)
      {
         printf("Configuration no %d fully processed in %.2e sec ",
                nc,wt2-wt1);
         printf("(average = %.2e sec)\n\n",
                wtavg/(double)((nc-first)/step+1));
      }

      check_endflag(&iend);

      if (my_rank==0)
      {
         fflush(flog);
         copy_file(log_file,log_save);
         copy_file(dat_file,dat_save);
         copy_file(rng_file,rng_save);
      }
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
