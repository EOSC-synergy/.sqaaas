
/*******************************************************************************
*
* File ms6.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* Based on openQCD-1.6/main/ms4.c
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of CA_mu(y0,x0) and CP(y0,x0) correlators.
*
* Syntax: ms6 -i <input file> [-noexp] [-a]
*
* For the definition of the correlators and for usage instructions see the file 
* README.ms6.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "uflds.h"
#include "u1flds.h"
#include "sflds.h"
#include "linalg.h"
#include "dirac.h"
#include "sap.h"
#include "dfl.h"
#include "forces.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static struct
{
   int tmax,x0,nsrc;
   int nfl,qhat[2];
   double kappa[2],mu[2];
   double su3csw[2],u1csw[2];
   double cF[2],cF_prime[2];
   double th1[2],th2[2],th3[2];
} file_head;

static struct
{
   int nc;
   double *P,*A0,*A1,*A2,*A3;
} data;

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static int my_rank,noexp,append,endian;
static int first,last,step;
static int ipgrd[2],level,seed;
static int *rlxs_state=NULL,*rlxd_state=NULL;
static double mus;

static int isx_re[L0],isx_im[L0],init=0;
static double corr_re[2*N0],corr_im[2*N0];

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],nbase[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL,*fdat=NULL,*fend=NULL;


static void alloc_data(void)
{
   int tmax;
   double *p;

   tmax=file_head.tmax;

   p=amalloc(5*tmax*sizeof(*p),4);
   error(p==NULL,1,"alloc_data [ms6.c]","Unable to allocate data arrays");

   data.P=p;
   data.A0=p+tmax;
   data.A1=p+2*tmax;
   data.A2=p+3*tmax;
   data.A3=p+4*tmax;
}


static void write_file_head(void)
{
   int ifl;

   write_little_int(fdat,4,
		    file_head.tmax,
		    file_head.x0,
		    file_head.nsrc,
		    file_head.nfl);

   for (ifl=0;ifl<file_head.nfl;++ifl)
   {
      write_little_int(fdat,1,file_head.qhat[ifl]);

      write_little_dble(fdat,9,
			file_head.kappa[ifl],
			file_head.mu[ifl],
			file_head.su3csw[ifl],
			file_head.u1csw[ifl],
			file_head.cF[ifl],
			file_head.cF_prime[ifl],
			file_head.th1[ifl],
			file_head.th2[ifl],
			file_head.th3[ifl]);
   }
}


static void check_file_head(void)
{
   int ifl;

   check_fpar_int("ms6",fdat,4,
		  file_head.tmax,
		  file_head.x0,
		  file_head.nsrc,
		  file_head.nfl);

   for (ifl=0;ifl<file_head.nfl;++ifl)
   {
      check_fpar_int("ms6",fdat,1,file_head.qhat[ifl]);

      check_fpar_dble("ms6",fdat,9,
		      file_head.kappa[ifl],
		      file_head.mu[ifl],
		      file_head.su3csw[ifl],
		      file_head.u1csw[ifl],
		      file_head.cF[ifl],
		      file_head.cF_prime[ifl],
		      file_head.th1[ifl],
		      file_head.th2[ifl],
		      file_head.th3[ifl]);
   }
}


static void write_data(void)
{
   int tmax;

   tmax=file_head.tmax;

   write_little_int(fdat,1,data.nc);
   
   write_little_dblearray(fdat,tmax,data.P);
   write_little_dblearray(fdat,tmax,data.A0);
   write_little_dblearray(fdat,tmax,data.A1);
   write_little_dblearray(fdat,tmax,data.A2);
   write_little_dblearray(fdat,tmax,data.A3);
}


static int read_data(void)
{
   int ir,t,tmax;
   stdint_t istd[1];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);

   if (ir!=1)
      return 0;

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   data.nc=(int)(istd[0]);

   tmax=file_head.tmax;

   for (t=0;t<tmax;t++)
   {
      ir+=fread(dstd,sizeof(double),1,fdat);
      
      if (endian==BIG_ENDIAN)
	 bswap_double(1,dstd);

      data.P[t]=dstd[0];
   }

   for (t=0;t<tmax;t++)
   {
      ir+=fread(dstd,sizeof(double),1,fdat);
      
      if (endian==BIG_ENDIAN)
	 bswap_double(1,dstd);

      data.A0[t]=dstd[0];
   }
   for (t=0;t<tmax;t++)
   {
      ir+=fread(dstd,sizeof(double),1,fdat);
      
      if (endian==BIG_ENDIAN)
	 bswap_double(1,dstd);

      data.A1[t]=dstd[0];
   }
   for (t=0;t<tmax;t++)
   {
      ir+=fread(dstd,sizeof(double),1,fdat);
      
      if (endian==BIG_ENDIAN)
	 bswap_double(1,dstd);

      data.A2[t]=dstd[0];
   }
   for (t=0;t<tmax;t++)
   {
      ir+=fread(dstd,sizeof(double),1,fdat);
      
      if (endian==BIG_ENDIAN)
	 bswap_double(1,dstd);

      data.A3[t]=dstd[0];
   }

   error_root(ir!=(1+5*tmax),1,"read_data [ms6.c]",
              "Read error or incomplete data record");

   return 1;
}


static void save_data(void)
{
   if (my_rank==0)
   {
      fdat=fopen(dat_file,"ab");
      error_root(fdat==NULL,1,"save_data [ms6.c]",
                 "Unable to open data file");
      write_data();
      fclose(fdat);
   }
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
                 "read_dirs [ms6.c]","Improper configuration range");

      find_section("Random number generator");
      read_line("level","%d",&level);
      read_line("seed","%d",&seed);
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&level,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   if (noexp)
      error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,last,NPROC-1)>=NAME_SIZE,
                 1,"setup_files [ms6.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [ms6.c]","cnfg_dir name is too long");

   check_dir_root(log_dir);
   check_dir_root(dat_dir);
   error_root(name_size("%s/%s.ms6.log~",log_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms6.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms6.dat~",dat_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms6.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.ms6.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.ms6.par",dat_dir,nbase);
   sprintf(dat_file,"%s/%s.ms6.dat",dat_dir,nbase);
   sprintf(end_file,"%s/%s.ms6.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
}


static void read_flds_bc_lat_parms(void)
{
   int gg,ifl,bc,cs;
   double phi[2],phi_prime[2];
   char line[NAME_SIZE];

   if (my_rank==0)
   {
      find_section("Gauge group");
      read_line("gauge","%s",line);

      gg=0;
      if (strcmp(line,"SU(3)")==0)
         gg=1;
      else if (strcmp(line,"U(1)")==0)
         gg=2;
      else if (strcmp(line,"SU(3)xU(1)")==0)
         gg=3;
      else
         error_root(1,1,"read_flds_bc_lat_parms [ms6.c]",
                    "Unknown gauge group %s",line);

      find_section("Correlators");
      read_line("nfl","%d",&(file_head.nfl));
      error_root( (file_head.nfl!=1) && (file_head.nfl!=2),1,"read_flds_bc_lat_parms [ms6.c]",
		 "Number of flavours must be 1 or 2 %s",line);

   }
   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(file_head.nfl),1,MPI_INT,0,MPI_COMM_WORLD);

   set_flds_parms(gg,file_head.nfl);

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
         error_root(1,1,"read_flds_bc_lat_parms [ms4.c]",
                    "Unknown time boundary condition type %s",line);
      
      read_line("cstar","%d",&cs);

      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;
      if ((cs==0)&&(bc==1)&&((gg&1)!=0))
         read_dprms("phi",2,phi);
      if ((bc==1)||(bc==2))
      {
         if ((cs==0)&&((gg&1)!=0))
            read_dprms("phi'",2,phi_prime);
      }
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_bc_parms(bc,0,cs,phi,phi_prime);
   
   for (ifl=0;ifl<file_head.nfl;ifl++)
   {
      file_head.qhat[ifl]=0;
      file_head.su3csw[ifl]=0.0;
      file_head.u1csw[ifl]=0.0;
      file_head.cF[ifl]=0.0;
      file_head.cF_prime[ifl]=0.0;
      file_head.th1[ifl]=0.0;
      file_head.th2[ifl]=0.0;
      file_head.th3[ifl]=0.0;
      
      sprintf(line,"Flavour %d",ifl);
      if (my_rank==0)
      {
         find_section(line);
         read_line("kappa","%lf",file_head.kappa+ifl);
         read_line("mu","%lf",file_head.mu+ifl);
         if (gg==1)
         {
            read_line("su3csw","%lf",file_head.su3csw+ifl);
         }
         else if (gg==2)
         {
            read_line("qhat","%d",file_head.qhat+ifl);
            read_line("u1csw","%lf",file_head.u1csw+ifl);
         }
         else if (gg==3)
         {
            read_line("qhat","%d",file_head.qhat+ifl);         
            read_line("su3csw","%lf",file_head.su3csw+ifl);
            read_line("u1csw","%lf",file_head.u1csw+ifl);
         }
         if (bc!=3) read_line("cF","%lf",file_head.cF+ifl);
         if (bc==2) read_line("cF'","%lf",file_head.cF_prime+ifl);
         if (cs==0) read_line("theta","%lf %lf %lf",
			      file_head.th1+ifl,
			      file_head.th2+ifl,
			      file_head.th3+ifl);
      }
      
      MPI_Bcast(file_head.qhat+ifl,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.kappa+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.mu+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.su3csw+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.u1csw+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.cF+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.cF_prime+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.th1+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.th2+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(file_head.th3+ifl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      set_qlat_parms(ifl,
		     file_head.kappa[ifl],
		     file_head.qhat[ifl],
		     file_head.su3csw[ifl],
		     file_head.u1csw[ifl],
		     file_head.cF[ifl],
		     file_head.cF_prime[ifl],
		     file_head.th1[ifl],
		     file_head.th2[ifl],
		     file_head.th3[ifl]);
   }
   
   error( (file_head.nfl==2)&&(file_head.qhat[0]!=file_head.qhat[1]),1,
	  "read_flds_bc_lat_parms [ms6.c]",
	  "The two flavours must have the same electric charge");
   
   error( (file_head.nfl==2) &&
	  (file_head.kappa[0]==file_head.kappa[1]) &&
	  (file_head.mu[0]==file_head.mu[1]) &&
	  (file_head.qhat[0]==file_head.qhat[1]) &&
	  (file_head.su3csw[0]==file_head.su3csw[1]) &&
	  (file_head.u1csw[0]==file_head.u1csw[1]) &&
	  (file_head.cF[0]==file_head.cF[1]) &&
	  (file_head.cF_prime[0]==file_head.cF_prime[1]) &&
	  (file_head.th1[0]==file_head.th1[1]) &&
	  (file_head.th2[0]==file_head.th2[1]) &&
	  (file_head.th3[0]==file_head.th3[1]),1,
	  "read_flds_bc_lat_parms [ms6.c]",
	  "The two flavours are identical, in this case you have to set nfl=1");

   if (my_rank==0)
   {
      find_section("Source fields");
      read_line("x0","%d",&(file_head.x0));
      read_line("nsrc","%d",&(file_head.nsrc));

      error_root( (file_head.x0<0)||(file_head.x0>=N0),1,"read_fld_bc_lat_parms [ms6.c]",
                 "Specified time x0 is out of range");

      error_root(((file_head.x0==0)&&(bc!=3))||((file_head.x0==(N0-1))&&(bc==0)),1,
                 "read_bc_parms [ms6.c]","Incompatible choice of boundary "
                 "conditions and source time");
   }

   MPI_Bcast(&(file_head.x0),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(file_head.nsrc),1,MPI_INT,0,MPI_COMM_WORLD);

   file_head.tmax=N0;
   
   if (append)
      check_flds_bc_lat_parms(fdat);
   else
      write_flds_bc_lat_parms(fdat);
}


static void read_sap_parms(void)
{
   int bs[4];

   if (my_rank==0)
   {
      find_section("SAP");
      read_line("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3);
   }

   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   set_sap_parms(bs,1,4,5);
}


static void read_dfl_parms(void)
{
   int bs[4],Ns;
   int ninv,nmr,ncy,nkv,nmx,qhat;
   double kappa,mu,su3csw,u1csw,cF,cF_prime,th1,th2,th3,res;

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
      cF=cF_prime=0.0;
      th1=th2=th3=0;
      su3csw=u1csw=0.0;

      find_section("Deflation subspace generation");
      read_line("kappa","%lf",&kappa);
      read_line("mu","%lf",&mu);
      if(gauge()==1)
	 read_line("csw","%lf",&su3csw);
      else if(gauge()==2)
      {
	 read_line("qhat","%d",&qhat);
	 read_line("csw","%lf",&su3csw);	 
      }
      else if(gauge()==3)
      {
	 read_line("qhat","%d",&qhat);
	 read_line("su3csw","%lf",&su3csw);	 
	 read_line("u1csw","%lf",&u1csw);	 
      }
      if (bc_type()!=3) read_line("cF","%lf",&cF);
      if (bc_type()==2) read_line("cF'","%lf",&cF_prime);
      if (bc_cstar()==0) read_line("theta","%lf %lf %lf",&th1,&th2,&th3);

      read_line("ninv","%d",&ninv);
      read_line("nmr","%d",&nmr);
      read_line("ncy","%d",&ncy);
   }

   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&qhat,1,MPI_INT,0,MPI_COMM_WORLD);
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
}


static void read_solver(void)
{
   solver_parms_t sp;

   read_solver_parms(0);
   sp=solver_parms(0);

   if ((sp.solver==SAP_GCR)||(sp.solver==DFL_SAP_GCR))
      read_sap_parms();

   if (sp.solver==DFL_SAP_GCR)
      read_dfl_parms();
   
   if (append)
      check_solver_parms(fdat);
   else
      write_solver_parms(fdat);
}


static void read_infile(int argc,char *argv[])
{
   int ifile;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [ms6.c]",
                 "Syntax: ms6 -i <input file> [-noexp]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ms6.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");
      append=find_opt(argc,argv,"-a");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ms6.c]",
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
         fdat=fopen(par_file,"rb");
      else
         fdat=fopen(par_file,"wb");

      error_root(fdat==NULL,1,"read_infile [ms3.c]",
                 "Unable to open parameter file");
   }

   read_flds_bc_lat_parms();
   read_solver();

   if (my_rank==0)
   {
      fclose(fin);
      fclose(fdat);

      if (append==0)
         copy_file(par_file,par_save);
   }
}


static void check_old_log(int *fst,int *lst,int *stp)
{
   int ie,ic,isv;
   int fc,lc,dc,pc;
   int np[4],bp[4];

   fend=fopen(log_file,"r");
   error_root(fend==NULL,1,"check_old_log [ms6.c]",
              "Unable to open log file");

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;
   isv=0;

   while (fgets(line,NAME_SIZE,fend)!=NULL)
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

   fclose(fend);

   error_root((ie&0x1)!=0x0,1,"check_old_log [ms6.c]",
              "Incorrect read count");
   error_root((ie&0x2)!=0x0,1,"check_old_log [ms6.c]",
              "Configuration numbers are not equally spaced");
   error_root(isv==0,1,"check_old_log [ms6.c]",
              "Log file extends beyond the last configuration save");

   (*fst)=fc;
   (*lst)=lc;
   (*stp)=dc;
}


static void check_old_dat(int fst,int lst,int stp)
{
   int ie,ic;
   int fc,lc,dc,pc;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [ms6.c]",
              "Unable to open data file");

   check_file_head();

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;

   while (read_data()==1)
   {
      pc=lc;
      lc=data.nc;
      ic+=1;

      if (ic==1)
         fc=lc;
      else if (ic==2)
         dc=lc-fc;
      else if ((ic>2)&&(lc!=(pc+dc)))
         ie|=0x1;
   }

   fclose(fdat);

   error_root(ic==0,1,"check_old_dat [ms6.c]",
              "No data records found");
   error_root((ie&0x1)!=0x0,1,"check_old_dat [ms6.c]",
              "Configuration numbers are not equally spaced");
   error_root((fst!=fc)||(lst!=lc)||(stp!=dc),1,"check_old_dat [ms6.c]",
              "Configuration range is not as reported in the log file");
}


static void check_files(void)
{
   int fst,lst,stp;

   ipgrd[0]=0;
   ipgrd[1]=0;

   if (my_rank==0)
   {
      if (append)
      {
         check_old_log(&fst,&lst,&stp);
         check_old_dat(fst,lst,stp);

         error_root((fst!=lst)&&(stp!=step),1,"check_files [ms6.c]",
                    "Continuation run:\n"
                    "Previous run had a different configuration separation");
         error_root(first!=lst+step,1,"check_files [ms6.c]",
                    "Continuation run:\n"
                    "Configuration range does not continue the previous one");
      }
      else
      {
         fin=fopen(log_file,"r");
         fdat=fopen(dat_file,"rb");

         error_root((fin!=NULL)||(fdat!=NULL),1,"check_files [ms6.c]",
                    "Attempt to overwrite old *.log or *.dat file");

         fdat=fopen(dat_file,"wb");
         error_root(fdat==NULL,1,"check_files [ms6.c]",
                    "Unable to open data file");
         write_file_head();
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

      error_root(flog==NULL,1,"print_info [ms6.c]","Unable to open log file");
      printf("\n");

      if (append)
         printf("Continuation run\n\n");
      else
      {
	 printf("Computation of CA_mu(y0,x0) and CP(y0,x0) correlators\n");
	 printf("-----------------------------------------------------\n\n");
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
         print_lat_parms();
	 
	 printf("Random number generator:\n");
	 printf("level = %d, seed = %d\n\n",level,seed);
	 
	 printf("Source fields:\n");
	 printf("x0 = %d\n",file_head.x0);
	 printf("nsrc = %d\n\n",file_head.nsrc);
	 
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


static void wsize(int *nws,int *nwsd,int *nwv,int *nwvd)
{
   int nsd;
   solver_parms_t sp;

   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;

   sp=solver_parms(0);
   nsd=2;

   if (sp.solver==CGNE)
   {
      MAX(*nws,5);
      MAX(*nwsd,nsd+3);
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
   else
      error_root(1,1,"wsize [ms6.c]",
                 "Unknown or unsupported solver");
}

static void random_source(spinor_dble *eta)
{
   int i,y0,iy,ix;
   double twopi,r[12];
   complex_dble *c;

   twopi=8.0*atan(1.0);

   set_sd2zero(VOLUME,eta);
   y0=file_head.x0-cpr[0]*L0;

   if ((y0>=0)&&(y0<L0))
   {
      for (iy=0;iy<(L1*L2*L3);iy++)
      {
         ix=ipt[iy+y0*L1*L2*L3];
	 c=(complex_dble*)(eta+ix);
	 ranlxd(r,12);
	 for(i=0;i<12;++i)
	 {
	    c[i].re=cos(twopi*r[i]);
	    c[i].im=sin(twopi*r[i]);
	 }
      }
   }
}


static void solve_dirac(spinor_dble *eta,spinor_dble *psi,int *status)
{
   solver_parms_t sp;
   sap_parms_t sap;

   sp=solver_parms(0);

   if (sp.solver==CGNE)
   {
      mulg5_dble(VOLUME,eta);

      tmcg(sp.nmx,sp.res,mus,eta,eta,status);

      error_root(status[0]<0,1,"solve_dirac [ms6.c]",
                 "CGNE solver failed (status = %d)",status[0]);

      Dw_dble(-mus,eta,psi);
      mulg5_dble(VOLUME,psi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      sap_gcr(sp.nkv,sp.nmx,sp.res,mus,eta,psi,status);

      error_root(status[0]<0,1,"solve_dirac [ms6.c]",
                 "SAP_GCR solver failed (status = %d)",status[0]);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.res,mus,eta,psi,status);

      error_root((status[0]<0)||(status[1]<0),1,
                 "solve_dirac [ms6.c]","DFL_SAP_GCR solver failed "
                 "(status = %d,%d,%d)",status[0],status[1],status[2]);
   }
   else
      error_root(1,1,"solve_dirac [ms6.c]",
                 "Unknown or unsupported solver");
}


static void scalar_product(int idirac,spinor_dble *r,spinor_dble *s, complex_dble *spt)
{  
   spt->re=0.0;
   spt->im=0.0;

   if (idirac==0)
   {
      spt->re+=(_vector_prod_re((*r).c1,(*s).c1));
      spt->re+=(_vector_prod_re((*r).c2,(*s).c2));
      spt->re+=(_vector_prod_re((*r).c3,(*s).c3));
      spt->re+=(_vector_prod_re((*r).c4,(*s).c4));
      
      spt->im+=(_vector_prod_im((*r).c1,(*s).c1));
      spt->im+=(_vector_prod_im((*r).c2,(*s).c2));
      spt->im+=(_vector_prod_im((*r).c3,(*s).c3));
      spt->im+=(_vector_prod_im((*r).c4,(*s).c4));
   }
   else if (idirac==1)
   {
      spt->re-=(_vector_prod_re((*r).c1,(*s).c3));
      spt->re-=(_vector_prod_re((*r).c2,(*s).c4));
      spt->re-=(_vector_prod_re((*r).c3,(*s).c1));
      spt->re-=(_vector_prod_re((*r).c4,(*s).c2));
      
      spt->im-=(_vector_prod_im((*r).c1,(*s).c3));
      spt->im-=(_vector_prod_im((*r).c2,(*s).c4));
      spt->im-=(_vector_prod_im((*r).c3,(*s).c1));
      spt->im-=(_vector_prod_im((*r).c4,(*s).c2));
   }
   else if (idirac==2)
   {
      spt->re+=(_vector_prod_im((*r).c1,(*s).c4));
      spt->re+=(_vector_prod_im((*r).c2,(*s).c3));
      spt->re-=(_vector_prod_im((*r).c3,(*s).c2));
      spt->re-=(_vector_prod_im((*r).c4,(*s).c1));
      
      spt->im-=(_vector_prod_re((*r).c1,(*s).c4));
      spt->im-=(_vector_prod_re((*r).c2,(*s).c3));
      spt->im+=(_vector_prod_re((*r).c3,(*s).c2));
      spt->im+=(_vector_prod_re((*r).c4,(*s).c1));
   }
   else if (idirac==3)
   {
      spt->re-=(_vector_prod_re((*r).c1,(*s).c4));
      spt->re+=(_vector_prod_re((*r).c2,(*s).c3));
      spt->re+=(_vector_prod_re((*r).c3,(*s).c2));
      spt->re-=(_vector_prod_re((*r).c4,(*s).c1));
      
      spt->im-=(_vector_prod_im((*r).c1,(*s).c4));
      spt->im+=(_vector_prod_im((*r).c2,(*s).c3));
      spt->im+=(_vector_prod_im((*r).c3,(*s).c2));
      spt->im-=(_vector_prod_im((*r).c4,(*s).c1));
   }
   else
   {
      spt->re+=(_vector_prod_im((*r).c1,(*s).c3));
      spt->re-=(_vector_prod_im((*r).c2,(*s).c4));
      spt->re-=(_vector_prod_im((*r).c3,(*s).c1));
      spt->re+=(_vector_prod_im((*r).c4,(*s).c2));
      
      spt->im-=(_vector_prod_re((*r).c1,(*s).c3));
      spt->im+=(_vector_prod_re((*r).c2,(*s).c4));
      spt->im+=(_vector_prod_re((*r).c3,(*s).c1));
      spt->im-=(_vector_prod_re((*r).c4,(*s).c2));
   }
}


void slices(int idirac,spinor_dble *psi1,spinor_dble *psi2)
{
   int bc,ix,t,t0,tmx;
   complex_dble spt;

   if (init==0)
   {
      for (t=0;t<L0;t++)
      {
         isx_re[t]=init_hsum(1);
         isx_im[t]=init_hsum(1);
      }

      init=1;
   }

   bc=bc_type();
   if (bc==0)
      tmx=N0-1;
   else
      tmx=N0;
   t0=cpr[0]*L0;

   for (t=0;t<L0;t++)
   {
      reset_hsum(isx_re[t]);
      reset_hsum(isx_im[t]);
   }

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (((t>0)&&(t<tmx))||(bc==3))
      {
         t-=t0;
         scalar_product(idirac,psi1+ix,psi2+ix,&spt);
         add_to_hsum(isx_re[t],&(spt.re));
         add_to_hsum(isx_im[t],&(spt.im));
      }
   }

   for (t=0;t<2*N0;t++)
   {
      corr_re[t]=0.0;
      corr_im[t]=0.0;
   }

   for (t=0;t<L0;t++)
   {
      local_hsum(isx_re[t],&(spt.re));
      local_hsum(isx_im[t],&(spt.im));
      corr_re[t+t0+N0]=spt.re;
      corr_im[t+t0+N0]=spt.im;
   }

   if (NPROC>1)
   {
      MPI_Reduce(corr_re+N0,corr_re,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(corr_re,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);

      MPI_Reduce(corr_im+N0,corr_im,N0,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(corr_im,N0,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   else
   {
      for (t=0;t<N0;t++)
      {
         corr_re[t]=corr_re[t+N0];
         corr_im[t]=corr_im[t+N0];
      }
   }
}



static void correlators(int nc)
{
   int t;
   int isrc,ifl,ifl2,stat[6];
   double wt1,wt2,wt[2],wtsum,norm;
   spinor_dble **wsd;
   dirac_parms_t dp;
   dfl_parms_t dfl;

   if(file_head.nfl==1)
   {
      wsd=reserve_wsd(2);
      ifl2=1;
   }
   else
   {
      wsd=reserve_wsd(3);
      ifl2=2;
   }

   norm=(double)(file_head.nsrc*N1*N2*N3);

   data.nc=nc;
   for (t=0;t<file_head.tmax;++t)
   {
      data.P[t]=0.0;
      data.A0[t]=0.0;
      data.A1[t]=0.0;
      data.A2[t]=0.0;
      data.A3[t]=0.0;
   }

   dfl=dfl_parms();
   if (dfl.Ns)
   {
      dfl_modes(stat);
      error_root(stat[0]<0,1,"main [ms6.c]",
		 "Deflation subspace generation failed (status = %d)",
		 stat[0]);
      
      if (my_rank==0)
	 printf("Deflation subspace generation: status = %d\n",stat[0]);
   }

   wtsum=0.0;
   for (isrc=0;isrc<file_head.nsrc;isrc++)
   {
      random_source(wsd[0]);

      for (ifl=0;ifl<file_head.nfl;++ifl)
      {
	 MPI_Barrier(MPI_COMM_WORLD);
	 wt1=MPI_Wtime();

	 dp=qlat_parms(ifl);
	 set_dirac_parms1(&dp);
	 mus=file_head.mu[ifl];
	 solve_dirac(wsd[0],wsd[1+ifl],stat+3*ifl);

	 MPI_Barrier(MPI_COMM_WORLD);
	 wt2=MPI_Wtime();
	 wt[ifl]=(wt2-wt1);
	 wtsum+=wt[ifl];
      }

      if (my_rank==0)
      {
         printf("Source n. %d, computation of propagators completed\n",isrc);

	 for (ifl=0;ifl<file_head.nfl;++ifl)
	 {
	    printf("propagator_%d: ",ifl+1);
	    if (dfl.Ns)
	    {
	       printf("status = %d,%d",stat[0+3*ifl],stat[1+3*ifl]);

	       if (stat[2+3*ifl])
		  printf(" (no of subspace regenerations = %d) ",stat[2+3*ifl]);
	       else
		  printf(" ");
	    }
	    else
	       printf("status = %d ",stat[0+3*ifl]);
	    printf("time %.2e sec\n",wt[ifl]);
	 }
	 printf("\n");
	 fflush(flog);
      }

      slices(0,wsd[1],wsd[ifl2]);
      for (t=0;t<file_head.tmax;++t)
	 data.P[t]+=corr_re[t]/norm;

      slices(1,wsd[1],wsd[ifl2]);
      for (t=0;t<file_head.tmax;++t)
	 data.A0[t]+=corr_re[t]/norm;

      slices(2,wsd[1],wsd[ifl2]);
      for (t=0;t<file_head.tmax;++t)
	 data.A1[t]+=corr_im[t]/norm;

      slices(3,wsd[1],wsd[ifl2]);
      for (t=0;t<file_head.tmax;++t)
	 data.A2[t]+=corr_im[t]/norm;

      slices(4,wsd[1],wsd[ifl2]);
      for (t=0;t<file_head.tmax;++t)
	 data.A3[t]+=corr_im[t]/norm;
   }
   
   wtsum/=(double)(file_head.nsrc);
   if(file_head.nfl==2)
      wtsum*=0.5;
   
   if (my_rank==0)
   {
      printf("Average time per inversion %.2e sec\n\n",wtsum);
      fflush(flog);
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

      error(rlxs_state==NULL,1,"save_ranlux [ms6.c]",
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


static void print_log(void)
{
   int t;

   if (my_rank==0)
      for (t=0;t<file_head.tmax;++t)
         printf("t = %3d  CP = %14.6e   CA_0 = %14.6e   CA_1 = %14.6e   CA_2 = %14.6e   CA_3 = %14.6e\n",
                t,data.P[t],data.A0[t],data.A1[t],data.A2[t],data.A3[t]);
}


static void check_endflag(int *iend)
{
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
   int nc,iend;
   int nws,nwsd,nwv,nwvd,n;
   double wt1,wt2,wtavg;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   alloc_data();
   check_files();
   print_info();

   if (append)
      start_ranlux(level,seed^(first-1));
   else
      start_ranlux(level,seed);

   geometry();

   wsize(&nws,&nwsd,&nwv,&nwvd);
   alloc_ws(nws);
   alloc_wsd(nwsd);
   alloc_wv(nwv);
   alloc_wvd(nwvd);

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
	 
         import_cnfg(cnfg_file);
	 
	 udfld(); 
	 set_flags(UPDATED_UD); 
	 adfld(); 
	 set_flags(UPDATED_AD);
      }
      
      correlators(nc);
      save_data();
      print_log();

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtavg+=(wt2-wt1);

      if (my_rank==0)
      {
         n=(nc-first)/step+1;

         printf("Configuration no %d fully processed in %.2e sec ",
                nc,wt2-wt1);
         printf("(average = %.2e sec)\n\n",wtavg/(double)(n));

         fflush(flog);
         copy_file(log_file,log_save);
         copy_file(dat_file,dat_save);
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
