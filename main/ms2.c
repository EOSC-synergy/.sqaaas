
/*******************************************************************************
*
* File ms2.c
*
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*               2017, 2019       Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the spectral range of the hermitian Dirac operator.
*
* Syntax: ms2 -i <input file> [-noexp] [-a]
*
* For usage instructions see the file README.ms2.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "sap.h"
#include "dfl.h"
#include "ratfcts.h"
#include "forces.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static int my_rank,noexp,endian,append;
static int first,last,step;
static int n_np_ra,n_np_rb;
static int *np_ra,*np_rb;
static int *rlxs_state=NULL,*rlxd_state=NULL;

static char log_dir[NAME_SIZE],loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char cnfg_file[NAME_SIZE],nbase[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL;


static void read_dirs(void)
{
   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);

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
                 "read_dirs [ms2.c]","Improper configuration range");
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
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
                 1,"setup_files [ms2.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [ms2.c]","cnfg_dir name is too long");

   check_dir_root(log_dir);
   error_root(name_size("%s/%s.ms2.log~",log_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms2.c]","log_dir name is too long");

   sprintf(log_file,"%s/%s.ms2.log",log_dir,nbase);
   sprintf(end_file,"%s/%s.ms2.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
}


static void read_flds_bc_lat_parms(void)
{
   int gg,nfl;
   char line[NAME_SIZE];

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
}


static void read_solvers(void)
{
   int nfl,ifl,isp;
   int isap,idfl;
   solver_parms_t sp;

   nfl=flds_parms().nfl;
   isap=0;
   idfl=0;

   for (ifl=0;ifl<nfl;ifl++)
   {
      isp=ifl;
      sp=solver_parms(isp);

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

   if (isap)
      read_sap_parms();

   if (idfl)
      read_dfl_parms(-1);
}


int intcmp (const void *a, const void *b) {
   return (*(int*)a - *(int*)b);
}


static void read_infile(int argc,char *argv[])
{
   int ifile,k;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [ms2.c]",
                 "Syntax: ms2 -i <input file> [-noexp] [-a]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ms2.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");
      append=find_opt(argc,argv,"-a");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ms2.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();
   read_flds_bc_lat_parms();

   if (my_rank==0)
   {
      find_section("Power method");
      n_np_ra=count_tokens("np_ra");
      n_np_rb=count_tokens("np_rb");
      error_root((n_np_ra<1),1,"read_infile [ms2.c]",
                 "np_ra token empty or absent");
      error_root((n_np_rb<1),1,"read_infile [ms2.c]",
                 "np_rb token empty or absent");
   }

   MPI_Bcast(&n_np_ra,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&n_np_rb,1,MPI_INT,0,MPI_COMM_WORLD);

   np_ra=malloc(sizeof(int)*n_np_ra);
   np_rb=malloc(sizeof(int)*n_np_rb);
   error((np_ra==NULL)||(np_rb==NULL),1,"read_infile [ms2.c]",
         "Unable to allocate np_ra or np_rb arrays");

   if (my_rank==0)
   {
      read_iprms("np_ra",n_np_ra,np_ra);
      read_iprms("np_rb",n_np_rb,np_rb);

      for(k=0;k<n_np_ra;k++)
      {
         error_root((np_ra[k]<1),1,"read_infile [ms2.c]",
                    "Power method iteration numbers must be positive integers");
      }
      for(k=0;k<n_np_rb;k++)
      {
         error_root((np_rb[k]<1),1,"read_infile [ms2.c]",
                    "Power method iteration numbers must be positive integers");
      }
   }

   MPI_Bcast(np_ra,n_np_ra,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(np_rb,n_np_rb,MPI_INT,0,MPI_COMM_WORLD);

   qsort(np_ra,n_np_ra,sizeof(int),&intcmp);
   qsort(np_rb,n_np_rb,sizeof(int),&intcmp);

   read_solvers();

   if (my_rank==0)
      fclose(fin);
}


static void check_files(void)
{
   if (my_rank==0)
   {
      if(!append)
      {
         fin=fopen(log_file,"r");
         error_root(fin!=NULL,1,"check_files [ms2.c]",
                    "Attempt to overwrite old *.log file");
      }
   }
}


static void print_info(void)
{
   int isap,idfl,k;
   long ip;

   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0L)
         remove("STARTUP_ERROR");

      flog=freopen(log_file,"a",stdout);
      error_root(flog==NULL,1,"print_info [ms2.c]","Unable to open log file");
      printf("\n");

      printf("Spectral range of the hermitian Dirac operator\n");
      printf("----------------------------------------------\n\n");

      printf("Program version %s\n",openQCD_RELEASE);

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noexp)
         printf("Configurations are read in imported file format\n\n");
      else
         printf("Configurations are read in exported file format\n\n");

      printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
      printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d process block size\n",
             NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);
      printf("SF boundary conditions on the quark fields\n\n");

      print_bc_parms();
      print_flds_parms();
      print_lat_parms();

      printf("Power method:\n");
      printf("np_ra =");
      for(k=0;k<n_np_ra;k++)
         printf(" %d",np_ra[k]);
      printf("\nnp_rb =");
      for(k=0;k<n_np_rb;k++)
         printf(" %d",np_rb[k]);
      printf("\n\n");

      print_solver_parms(&isap,&idfl);

      if (isap)
         print_sap_parms(0);

      if (idfl)
         print_dfl_parms(0);

      printf("Configurations no %d -> %d in steps of %d\n\n",
             first,last,step);
      fflush(flog);
   }
}


static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dfp;
   dfl_pro_parms_t dpp;

   dfp=dfl_parms();
   dpp=dfl_pro_parms();

   MAX(*nws,dfp.Ns+2);
   MAX(*nwv,2*dpp.nkv+2);
   MAX(*nwvd,4);
}


static void wsize(int *nws,int *nwsd,int *nwv,int *nwvd)
{
   int nsd,ifl,nfl,isp;
   solver_parms_t sp;

   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;

   nfl=flds_parms().nfl;
   for (ifl=0;ifl<nfl;ifl++)
   {
      isp=ifl;
      sp=solver_parms(isp);

      if (sp.solver==CGNE)
      {
         nsd=1;
         MAX(*nws,5);
         MAX(*nwsd,nsd+3);
      }
      else if (sp.solver==SAP_GCR)
      {
         nsd=2;
         MAX(*nws,2*sp.nkv+1);
         MAX(*nwsd,nsd+2);
      }
      else if (sp.solver==DFL_SAP_GCR)
      {
         nsd=2;
         MAX(*nws,2*sp.nkv+2);
         MAX(*nwsd,nsd+4);
         dfl_wsize(nws,nwv,nwvd);
      }
      else
         error_root(1,1,"wsize [ms2.c]",
                    "Unknown or unsupported solver");
   }
}


static void power1(double *ra,int ifl,int *status)
{
   int nsd,j,k,l,stat[6];
   double r;
   spinor_dble **wsd;
   dirac_parms_t dp;
   solver_parms_t sp;
   sap_parms_t sap;

   dp=qlat_parms(ifl);
   set_dirac_parms1(&dp);
   sp=solver_parms(ifl);

   if (sp.solver==CGNE)
   {
      nsd=1;
      status[0]=0;
   }
   else if (sp.solver==SAP_GCR)
   {
      nsd=2;
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);
      status[0]=0;
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      nsd=2;
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      for (l=0;l<3;l++)
         status[l]=0;
   }
   else
   {
      nsd=1;
      error_root(1,1,"power1 [ms2.c]",
                 "Unknown or unsupported solver");
   }

   wsd=reserve_wsd(nsd);
   random_sd(VOLUME/2,wsd[0],1.0);
   bnd_sd2zero(EVEN_PTS,wsd[0]);
   r=normalize_dble(VOLUME/2,1,wsd[0]);

   j=0;
   for (k=0;k<np_ra[n_np_ra-1];k++)
   {
      if (sp.solver==CGNE)
      {
         tmcgeo(sp.nmx,sp.res,0.0,wsd[0],wsd[0],stat);

         error_root(stat[0]<0,1,"power1 [ms2.c]",
                    "CGNE solver failed (status = %d)",stat[0]);

         if (status[0]<stat[0])
            status[0]=stat[0];
      }
      else if (sp.solver==SAP_GCR)
      {
         mulg5_dble(VOLUME/2,wsd[0]);
         set_sd2zero(VOLUME/2,wsd[0]+(VOLUME/2));
         sap_gcr(sp.nkv,sp.nmx,sp.res,0.0,wsd[0],wsd[1],stat);
         mulg5_dble(VOLUME/2,wsd[1]);
         set_sd2zero(VOLUME/2,wsd[1]+(VOLUME/2));
         sap_gcr(sp.nkv,sp.nmx,sp.res,0.0,wsd[1],wsd[0],stat+1);

         error_root((stat[0]<0)||(stat[1]<0),1,"power1 [ms2.c]",
                    "SAP_GCR solver failed (status = %d;%d)",
                    stat[0],stat[1]);

         for (l=0;l<2;l++)
         {
            if (status[0]<stat[l])
               status[0]=stat[l];
         }
      }
      else if (sp.solver==DFL_SAP_GCR)
      {
         mulg5_dble(VOLUME/2,wsd[0]);
         set_sd2zero(VOLUME/2,wsd[0]+(VOLUME/2));
         dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,0.0,wsd[0],wsd[1],stat);
         mulg5_dble(VOLUME/2,wsd[1]);
         set_sd2zero(VOLUME/2,wsd[1]+(VOLUME/2));
         dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,0.0,wsd[1],wsd[0],stat+3);

         error_root((stat[0]<0)||(stat[1]<0)||(stat[3]<0)||(stat[4]<0),1,
                    "power1 [ms2.c]","DFL_SAP_GCR solver failed "
                    "(status = %d,%d,%d;%d,%d,%d)",
                    stat[0],stat[1],stat[2],stat[3],stat[4],stat[5]);

         for (l=0;l<2;l++)
         {
            if (status[l]<stat[l])
               status[l]=stat[l];

            if (status[l]<stat[l+3])
               status[l]=stat[l+3];
         }

         status[2]+=(stat[2]!=0);
         status[2]+=(stat[5]!=0);
      }

      r=normalize_dble(VOLUME/2,1,wsd[0]);

      if(k+1==np_ra[j])
      {
         ra[j]=1.0/sqrt(r);
         j++;
      }
   }

   release_wsd();
}


static void power2(double *rb,int ifl)
{
   int j,k;
   double r;
   spinor_dble **wsd;
   dirac_parms_t dp;

   dp=qlat_parms(ifl);
   set_dirac_parms1(&dp);
   sw_term(ODD_PTS);

   wsd=reserve_wsd(2);
   random_sd(VOLUME/2,wsd[0],1.0);
   bnd_sd2zero(EVEN_PTS,wsd[0]);
   r=normalize_dble(VOLUME/2,1,wsd[0]);

   j=0;
   for (k=0;k<np_rb[n_np_rb-1];k++)
   {
      Dwhat_dble(0.0,wsd[0],wsd[1]);
      mulg5_dble(VOLUME/2,wsd[1]);
      Dwhat_dble(0.0,wsd[1],wsd[0]);
      mulg5_dble(VOLUME/2,wsd[0]);

      r=normalize_dble(VOLUME/2,1,wsd[0]);

      if(k+1==np_rb[j])
      {
         rb[j]=sqrt(r);
         j++;
      }
   }

   release_wsd();
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

      error(rlxs_state==NULL,1,"save_ranlux [ms2.c]",
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
   int k,nc,iend,status[3];
   int nws,nwsd,nwv,nwvd,idfl,nfl,ifl;
   int cnfg_type;
   double *ra,*rb;
   double wt1,wt2,wtavg;
   dfl_parms_t dfl;
   dflst_t dfl_status;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   check_files();
   print_info();
   dfl=dfl_parms();

   geometry();
   start_ranlux(0,1234);

   wsize(&nws,&nwsd,&nwv,&nwvd);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);

   ra=malloc(sizeof(double)*n_np_ra);
   rb=malloc(sizeof(double)*n_np_rb);
   error((ra==NULL)||(rb==NULL),1,"main [ms2.c]",
         "Unable to allocate ra or rb arrays");

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
         error(cnfg_type!=gauge(),1,"main [ms2.c]",
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
               dfl_modes(idfl,status);
               error_root(status[0]<0,1,"main [ms2.c]",
                          "Generation of deflation subspace %d failed (status = %d)",
                          idfl,status[0]);
            }
            idfl++;
         }
      }

      nfl=flds_parms().nfl;
      for (ifl=0;ifl<nfl;ifl++)
      {
         power1(ra,ifl,status);
         power2(rb,ifl);

         if (my_rank==0)
         {
            printf("nc= %6d   ifl= %3d   ra=",nc,ifl);
            for(k=0;k<n_np_ra;k++)
               printf(" %.6e",ra[k]);
            printf("\n");

            printf("nc= %6d   ifl= %3d   rb=",nc,ifl);
            for(k=0;k<n_np_rb;k++)
               printf(" %.6e",rb[k]);
            printf("\n");

            printf("nc= %6d   ifl= %3d   status= ",nc,ifl);
            if (dfl.Ns)
               printf("%d,%d,%d\n",
                      status[0],status[1],status[2]);
            else
               printf("%d\n",status[0]);
         }
      }

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
