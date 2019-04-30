
/*******************************************************************************
*
* File ms3.c
*
* Copyright (C) 2012, 2013 Martin Luescher
*               2017       Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of SU(3) Wilson-flow observables.
*
* Syntax: ms3 -i <input file> [-noexp] [-a]
*
* For usage instructions see the file README.ms3.
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
#include "version.h"
#include "global.h"
#include "lib/main_utils.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int my_rank,noexp,append,endian;
static int first,last,step;
static int ipgrd[2];

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char ms3dat_file[NAME_SIZE],ms3dat_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],nbase[NAME_SIZE];
static FILE *flog=NULL;


static void read_dirs(void)
{
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

      error_root((last<first)||(step<1)||(((last-first)%step)!=0),1,
                 "read_dirs [ms3.c]","Improper configuration range");
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   if (noexp)
      error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,last,NPROC-1)>=NAME_SIZE,
                 1,"setup_files [ms3.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [ms3.c]","cnfg_dir name is too long");

   check_dir_root(log_dir);
   check_dir_root(dat_dir);
   error_root(name_size("%s/%s.ms3.log~",log_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms3.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms3.dat~",dat_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms3.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.ms3.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.ms3.par",dat_dir,nbase);
   sprintf(ms3dat_file,"%s/%s.ms3.dat",dat_dir,nbase);
   sprintf(end_file,"%s/%s.ms3.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(ms3dat_save,"%s~",ms3dat_file);
}


static void read_flds_bc_lat_parms(FILE *fpar)
{
   set_flds_parms(1,0);
   read_bc_parms();

   if (append)
      check_flds_bc_lat_parms(fpar);
   else
      write_flds_bc_lat_parms(fpar);
}


static void read_observables(FILE *fpar)
{
   read_wflow_parms();
   if (append)
      check_wflow_parms(fpar);
   else
      write_wflow_parms(fpar);
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

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [ms3.c]",
                 "Syntax: ms3 -i <input file> [-noexp] [-a]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ms3.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");
      append=find_opt(argc,argv,"-a");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ms3.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();

   if (my_rank==0)
   {
      if (append)
         fpar=fopen(par_file,"rb");
      else
         fpar=fopen(par_file,"wb");

      error_root(fpar==NULL,1,"read_infile [ms3.c]",
                 "Unable to open parameter file");
   }

   read_flds_bc_lat_parms(fpar);
   read_observables(fpar);

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
   error_root(fold==NULL,1,"check_old_log [ms3.c]",
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

   error_root((ie&0x1)!=0x0,1,"check_old_log [ms3.c]",
              "Incorrect read count");
   error_root((ie&0x2)!=0x0,1,"check_old_log [ms3.c]",
              "Configuration numbers are not equally spaced");
   error_root(isv==0,1,"check_old_log [ms3.c]",
              "Log file extends beyond the last configuration save");

   (*fst)=fc;
   (*lst)=lc;
   (*stp)=dc;
}


static void check_old_ms3dat(int fst,int lst,int stp)
{
   int ie,ic;
   int fc,lc,dc,pc;
   FILE *fdat=NULL;

   fdat=fopen(ms3dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [ms3.c]",
              "Unable to open data file");

   check_ms3dat_head(fdat);

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;

   while (read_ms3dat(fdat,&lc)==1)
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

   error_root(ic==0,1,"check_old_dat [ms3.c]",
              "No data records found");
   error_root((ie&0x1)!=0x0,1,"check_old_dat [ms3.c]",
              "Configuration numbers are not equally spaced");
   error_root((fst!=fc)||(lst!=lc)||(stp!=dc),1,"check_old_dat [ms3.c]",
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
         check_old_ms3dat(fst,lst,stp);

         error_root((fst!=lst)&&(stp!=step),1,"check_files [ms3.c]",
                    "Continuation run:\n"
                    "Previous run had a different configuration separation");
         error_root(first!=lst+step,1,"check_files [ms3.c]",
                    "Continuation run:\n"
                    "Configuration range does not continue the previous one");
      }
      else
      {
         error_root(fopen(log_file,"r")!=NULL,1,
                    "check_files [ms3.c]",
                    "Attempt to overwrite old *.log file");

         error_root(fopen(ms3dat_file,"rb")!=NULL,1,
                    "check_files [ms3.c]",
                    "Attempt to overwrite old *.ms*.dat file");

         fdat=fopen(ms3dat_file,"wb");
         error_root(fdat==NULL,1,"check_files [ms3.c]",
                    "Unable to open measurement data file");
         write_ms3dat_head(fdat);
         fclose(fdat);
      }
   }
}


static void print_info(void)
{
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

      error_root(flog==NULL,1,"print_info [ms3.c]","Unable to open log file");
      printf("\n");

      if (append)
         printf("Continuation run\n\n");
      else
      {
         printf("Computation of SU(3) Wilson flow observables\n");
         printf("--------------------------------------------\n\n");
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

      if (append==0)
      {
         print_bc_parms();
         print_flds_parms();
         print_wflow_parms();
      }

      printf("Configurations no %d -> %d in steps of %d\n\n",
             first,last,step);
      fflush(flog);
   }
}


static void save_msdat(void)
{
   FILE *fdat=NULL;

   if (my_rank==0)
   {
      fdat=fopen(ms3dat_file,"ab");
      error_root(fdat==NULL,1,"save_msdat [ms3.c]",
                 "Unable to open data file");
      write_ms3dat(fdat);
      fclose(fdat);
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


int main(int argc,char *argv[])
{
   int nc,iend,cnfg_type;
   double wt1,wt2,wtavg;
   wflow_parms_t wfp;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   check_files();
   print_info();

   geometry();
   
   wfp=wflow_parms();
   alloc_wud(1);
   if (wfp.flint)
      alloc_wf3d(1);

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
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,nc,my_rank);
         read_cnfg(cnfg_file);
      }
      else
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,nc);
         cnfg_type=import_cnfg(cnfg_file);
         error((cnfg_type&1)==0,1,"main [ms3.c]",
               "Imported configuration does not contain an SU(3) gauge field");
      }

      set_ms3dat(nc);
      save_msdat();
      print_ms3dat();

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtavg+=(wt2-wt1);

      if (my_rank==0)
      {
         printf("Configuration no %d fully processed in %.2e sec ",
                nc,wt2-wt1);
         printf("(average = %.2e sec)\n\n",
                wtavg/(double)((nc-first)/step+1));
         fflush(flog);

         copy_file(log_file,log_save);
         copy_file(ms3dat_file,ms3dat_save);
      }

      check_endflag(&iend);
   }

   if (my_rank==0)
   {
      fflush(flog);
      copy_file(log_file,log_save);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
