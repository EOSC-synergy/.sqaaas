
/*******************************************************************************
*
* File check6.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* Based on openQCD-1.6/main/ms4.c
* Copyright (C) 2012, 2013, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge invariance of the CA_mu(y0,x0) and CP(y0,x0) correlators
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
#include "global.h"
#include "cstates.h"
#include "gflds_utils.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static struct
{
   int tmax,x0;
   int nfl,qhat[2];
   double kappa[2],mu[2];
   double su3csw[2],u1csw[2];
   double cF[2],cF_prime[2];
   double th1[2],th2[2],th3[2];
} file_head;

static struct
{
   double *P,*A0,*A1,*A2,*A3;
} data;

#define MAX(n,m) if ((n)<(m)) (n)=(m)

static int my_rank,rflag;
static double mus;

static int isx_re[L0],isx_im[L0],init=0;
static double corr_re[2*N0],corr_im[2*N0];

static char line[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL;

static spinor_dble **wsd;

static void alloc_data(void)
{
   int t,tmax;
   double *p;

   tmax=file_head.tmax;

   p=amalloc(10*tmax*sizeof(*p),4);
   error(p==NULL,1,"alloc_data [check6.c]","Unable to allocate data arrays");

   data.P=p;
   data.A0=p+tmax;
   data.A1=p+2*tmax;
   data.A2=p+3*tmax;
   data.A3=p+4*tmax;

   for (t=0;t<5*file_head.tmax;++t)
      data.P[t]=0.0;
}


static void read_flds_bc_lat_parms(void)
{
   int ifl;

   if (my_rank==0)
   {
      find_section("Gauge configuration");
      read_line("type","%s",&line);

      rflag=2;
      if ((strcmp(line,"cold")==0)||(strcmp(line,"0")==0))
         rflag=0;
      else if ((strcmp(line,"hot")==0)||(strcmp(line,"1")==0))
         rflag=1;
      else
      error_root(1,1,"read_flds_bc_lat_parms [ms4.c]",
                  "Unknown time of gauge configuration %s",line);

      find_section("Correlators");
      read_line("nfl","%d",&(file_head.nfl));
      error_root( (file_head.nfl!=1) && (file_head.nfl!=2),1,"read_flds_bc_lat_parms [check6.c]",
                  "Number of flavours must be 1 or 2 %s",line);
   }
   MPI_Bcast(&rflag,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(file_head.nfl),1,MPI_INT,0,MPI_COMM_WORLD);

   set_flds_parms(2,file_head.nfl);

   read_bc_parms();

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
         read_line("qhat","%d",file_head.qhat+ifl);
         read_line("u1csw","%lf",file_head.u1csw+ifl);

         if (bc_type()!=3) read_line("cF","%lf",file_head.cF+ifl);
         if (bc_type()==2) read_line("cF'","%lf",file_head.cF_prime+ifl);
         if (bc_cstar()==0) read_line("theta","%lf %lf %lf",
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

   error((file_head.nfl==2) &&
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
         "read_flds_bc_lat_parms [check6.c]",
         "The two flavours are identical, in this case you have to set nfl=1");

   if (my_rank==0)
   {
      find_section("Source fields");
      read_line("x0","%d",&(file_head.x0));

      error_root( (file_head.x0<0)||(file_head.x0>=N0),1,"read_fld_bc_lat_parms [check6.c]",
                  "Specified time x0 is out of range");

      error_root(((file_head.x0==0)&&(bc_type()!=3))||((file_head.x0==(N0-1))&&(bc_type()==0)),1,
                  "read_bc_parms [check6.c]","Incompatible choice of boundary "
                  "conditions and source time");
   }

   MPI_Bcast(&(file_head.x0),1,MPI_INT,0,MPI_COMM_WORLD);

   file_head.tmax=N0;
}


static void read_infile(void)
{
   solver_parms_t sp;
   int isap,idfl;

   if (my_rank==0)
   {
      flog=freopen("check6.log","w",stdout);
      error_root(flog==NULL,1,"print_info [check6.c]","Unable to open log file");
      printf("\n");

      fin=freopen("check6.in","r",stdin);
      error_root(fin==NULL,1,"read_infile [check6.c]",
                  "Unable to open input file");
   }

   read_flds_bc_lat_parms();
   read_solver_parms(0);
   sp=solver_parms(0);

   isap=0;
   idfl=0;
   if (sp.solver==SAP_GCR)
      isap=1;
   else if (sp.solver==DFL_SAP_GCR)
   {
      isap=1;
      idfl=1;
      
      if (dfl_gen_parms(sp.idfl).status!=DFL_DEF)
         read_dfl_parms(sp.idfl);
   }
   
   if (isap)
      read_sap_parms();
   
   if (idfl)
      read_dfl_parms(-1);

   if (my_rank==0)
      fclose(fin);
}


static void print_info(void)
{
   int isap,idfl;

   if (my_rank==0)
   {
      printf("Gauge invariance of the CA_mu(y0,x0) and CP(y0,x0) correlators\n");
      printf("--------------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
      printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d process block size\n\n",
      NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);

      print_bc_parms();
      print_flds_parms();
      print_lat_parms();

      if (rflag)
         printf("Random (hot) gauge configuration\n");
      else
         printf("Free (cold) gauge configuration\n");

      printf("Source fields:\n");
      printf("x0 = %d\n\n",file_head.x0);

      print_solver_parms(&isap,&idfl);

      if (isap)
         print_sap_parms(0);

      if (idfl)
         print_dfl_parms(0);

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
      error_root(1,1,"wsize [check6.c]",
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

      error_root(status[0]<0,1,"solve_dirac [check6.c]",
                  "CGNE solver failed (status = %d)",status[0]);

      Dw_dble(-mus,eta,psi);
      mulg5_dble(VOLUME,psi);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      sap_gcr(sp.nkv,sp.nmx,sp.res,mus,eta,psi,status);

      error_root(status[0]<0,1,"solve_dirac [check6.c]",
                  "SAP_GCR solver failed (status = %d)",status[0]);
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      dfl_sap_gcr2(sp.idfl,sp.nkv,sp.nmx,sp.res,mus,eta,psi,status);

      error_root((status[0]<0)||(status[1]<0),1,
                  "solve_dirac [check6.c]","DFL_SAP_GCR solver failed "
                  "(status = %d,%d,%d)",status[0],status[1],status[2]);
   }
   else
   error_root(1,1,"solve_dirac [check6.c]",
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



static void correlators(double sign)
{
   int t,coulomb;
   int ifl,ifl2,stat[6],idfl;
   dirac_parms_t dp;
   dfl_parms_t dfl;
   dflst_t dfl_status;

   dfl=dfl_parms();
   if (dfl.Ns)
   {
      idfl=0;
      while(1)
      {
         dfl_status=dfl_gen_parms(idfl).status;
         if(dfl_status==DFL_OUTOFRANGE) break;
         if(dfl_status==DFL_DEF)
         {
            dfl_modes2(idfl,stat);
            error_root(stat[0]<0,1,"main [ms6.c]",
                        "Generation of deflation subspace %d failed (status = %d)",
            idfl,stat[0]);

            if (my_rank==0)
               printf("Generation of deflation subspace %d: status = %d\n\n",idfl,stat[0]);
         }
         idfl++;
      }
   }

   coulomb=rflag;

   for (ifl=0;ifl<file_head.nfl;++ifl)
   {
      dp=qlat_parms(ifl);
      set_dirac_parms1(&dp);
      mus=file_head.mu[ifl];

      mul_cfactor_muaverage(1,coulomb,wsd[0],wsd[1]);
      solve_dirac(wsd[1],wsd[2+ifl],stat+3*ifl);
      mul_cfactor_muaverage(0,coulomb,wsd[2+ifl],wsd[2+ifl]);
   }

   if (my_rank==0)
   {
      printf("\n");
      printf("Computation of propagators completed\n");

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
         printf("\n");
      }
      printf("\n");
      fflush(flog);
   }


   if(file_head.nfl==1)
      ifl2=2;
   else
      ifl2=3;

   slices(0,wsd[2],wsd[ifl2]);
   for (t=0;t<file_head.tmax;++t)
      data.P[t]+=sign*corr_re[t];

   slices(1,wsd[2],wsd[ifl2]);
   for (t=0;t<file_head.tmax;++t)
      data.A0[t]+=sign*corr_re[t];

   slices(2,wsd[2],wsd[ifl2]);
   for (t=0;t<file_head.tmax;++t)
      data.A1[t]+=sign*corr_im[t];

   slices(3,wsd[2],wsd[ifl2]);
   for (t=0;t<file_head.tmax;++t)
      data.A2[t]+=sign*corr_im[t];

   slices(4,wsd[2],wsd[ifl2]);
   for (t=0;t<file_head.tmax;++t)
      data.A3[t]+=sign*corr_im[t];
}


int main(int argc,char *argv[])
{
   int t,nws,nwsd,nwv,nwvd;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile();
   alloc_data();
   print_info();

   start_ranlux(0,12345);

   geometry();

   wsize(&nws,&nwsd,&nwv,&nwvd);
   alloc_ws(nws);
   alloc_wsd(nwsd+4);
   alloc_wv(nwv);
   alloc_wvd(nwvd);

   wsd=reserve_wsd(4);
   random_source(wsd[0]);
   
   adfld();

   if(rflag)
   {
      random_gflds();
      orbi_cpy_ad();
   }

   correlators(1.0);

   if (my_rank==0)
   {
      for (t=0;t<file_head.tmax;++t)
      printf("t = %3d  CP = %14.6e   CA_0 = %14.6e   CA_1 = %14.6e   CA_2 = %14.6e   CA_3 = %14.6e\n",
               t,data.P[t],data.A0[t],data.A1[t],data.A2[t],data.A3[t]);
      fflush(flog);
   }

   for (t=0;t<5*file_head.tmax;++t)
      data.P[t+5*file_head.tmax]=data.P[t];

   random_g();
   orbi_cpy_g();
   transform_gflds();

   correlators(-1.0);

   if (my_rank==0)
   {
      for (t=0;t<file_head.tmax;++t)
         printf("t = %3d  dev_CP = %14.6e   deve_CA_0 = %14.6e   dev_CA_1 = %14.6e   dev_CA_2 = %14.6e   dev_CA_3 = %14.6e\n",
                  t,data.P[t],data.A0[t],data.A1[t],data.A2[t],data.A3[t]);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
