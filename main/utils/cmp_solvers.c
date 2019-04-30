
/*******************************************************************************
*
* File cmp_solvers.c
*
* Copyright (C) 2011-2013, 2016 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Compare performance of DFL+SAP+GCR and CG solvers.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "forces.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "uflds.h"
#include "u1flds.h"
#include "sflds.h"
#include "linalg.h"
#include "dirac.h"
#include "dfl.h"
#include "sap.h"
#include "global.h"

int my_rank;

int first,last,step;
int nsp,nsrc,*ifl=NULL,*eoflg=NULL;
double *mu=NULL;
char cnfg_dir[NAME_SIZE],cnfg_file[NAME_SIZE],nbase[NAME_SIZE];


void print_solver_parms2(void)
{
   int i;
   solver_parms_t sp;
   
   if (my_rank==0)
   {
      for (i=0;i<nsp;i++)
      {
         sp=solver_parms(i);
         printf("Solver %d:\n",i);
         
         if (sp.solver==CGNE)
         {
            printf("CGNE solver\n");
            printf("ifl = %d\n",ifl[i]);
            printf("mu = %.3e\n",mu[i]);
            printf("eoflg = %d\n",eoflg[i]);
            printf("nmx = %d\n",sp.nmx);
            printf("res = %.1e\n\n",sp.res);
         }
         else if (sp.solver==SAP_GCR)
         {
            printf("SAP_GCR solver\n");
            printf("ifl = %d\n",ifl[i]);
            printf("mu = %.3e\n",mu[i]);
            printf("eoflg = %d\n",eoflg[i]);
            printf("nkv = %d\n",sp.nkv);
            printf("isolv = %d\n",sp.isolv);
            printf("nmr = %d\n",sp.nmr);              
            printf("ncy = %d\n",sp.ncy);
            printf("nmx = %d\n",sp.nmx);
            printf("res = %.1e\n\n",sp.res);
         }
         else if (sp.solver==DFL_SAP_GCR)
         {
            printf("DFL_SAP_GCR solver\n");
            printf("ifl = %d\n",ifl[i]);
            printf("mu = %.3e\n",mu[i]);
            printf("eoflg = %d\n",eoflg[i]);
            printf("nkv = %d\n",sp.nkv);
            printf("isolv = %d\n",sp.isolv);
            printf("nmr = %d\n",sp.nmr);              
            printf("ncy = %d\n",sp.ncy);
            printf("nmx = %d\n",sp.nmx);
            printf("idfl = %d\n",sp.idfl);
            printf("res = %.1e\n\n",sp.res);
         }
         else
            printf("UNKNOWN solver\n\n");
      }
   }
}



void read_input(void)
{
   int gg,nfl;
   int isp;
   int isap,idfl;
   double invqel;
   char line[NAME_SIZE];
   solver_parms_t sp;

   if (my_rank==0)
   {
      find_section("Configurations");
      read_line("name","%s",nbase);
      read_line("cnfg_dir","%s",cnfg_dir);
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);

      find_section("Parameters");
      read_line("gauge","%s",line);
      read_line("n_flavours","%d",&nfl);
      read_line("n_solvers","%d",&nsp);
      read_line("n_sources","%d",&nsrc);

      gg=0;
      if (strcmp(line,"SU(3)")==0)
         gg=1;
      else if (strcmp(line,"U(1)")==0)
         gg=2;
      else if (strcmp(line,"SU(3)xU(1)")==0)
         gg=3;
      else
         error_root(1,1,"read_input [cmp_solvers.c]",
                    "Unknown gauge group %s",line);

      if((gg&2)!=0)
         read_line("invqel","%lf",&invqel);
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(&gg,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nfl,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nsp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&invqel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


   set_flds_parms(gg,nfl);
   read_bc_parms();
   set_u1lat_parms(0,DBL_MAX,invqel,0.0,1.0,1.0,1.0,0);
   read_qlat_parms();
   
   ifl=malloc(nsp*sizeof(*ifl));
   error(ifl==NULL,1,"read_input [cmp_solvers.c]",
         "Unable to allocate ifl array");
   
   mu=malloc(nsp*sizeof(*mu));
   error(mu==NULL,1,"read_input [cmp_solvers.c]",
         "Unable to allocate mu array");
   
   eoflg=malloc(nsp*sizeof(*eoflg));
   error(eoflg==NULL,1,"read_input [cmp_solvers.c]",
         "Unable to allocate eoflg array");

   isap=0;
   idfl=0;
   for (isp=0;isp<nsp;isp++)
   {
      read_solver_parms(isp);

      sp=solver_parms(isp);

      read_line("ifl","%d",ifl+isp);

      eoflg[isp]=1;
      if (sp.solver==CGNE)
         read_line("eoflg","%d",eoflg+isp);

      read_line("mu","%lf",mu+isp);

      error_root((ifl[isp]<0)||(ifl[isp]>=nfl),1,"read_input [cmp_solvers.c]",
                 "Solver %d: ifl is out of range",isp);

      error_root((eoflg[isp]<0)||(eoflg[isp]>1),1,"read_input [cmp_solvers.c]",
                 "Solver %d: eoflg is out of range",isp);

      error_root(sp.solver==MSCG,1,"read_input [cmp_solvers.c]",
                 "MSCG solver not supported");

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

   MPI_Bcast(ifl,nsp,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(eoflg,nsp,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(mu,nsp,MPI_DOUBLE,0,MPI_COMM_WORLD);

   if (isap)
      read_sap_parms();

   if (idfl)
      read_dfl_parms(-1);

   
      
   print_flds_bc_lat_parms();
   print_solver_parms2();

   if (isap)
      print_sap_parms(0);

   if (idfl)
      print_dfl_parms(0);
   
}



void alloc_workspace(void)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpr;
   solver_parms_t sp;
   int isp,mnkv,idfl,ns;

   dp=dfl_parms();
   dpr=dfl_pro_parms();

   idfl=0;
   mnkv=0;
   for (isp=0;isp<nsp;isp++)
   {
      sp=solver_parms(isp);
      if (sp.nkv>mnkv)
         mnkv=sp.nkv;
      if (sp.solver==DFL_SAP_GCR)
         idfl=1;
   }
   ns=5;
   if((2*mnkv+2)>ns) ns=2*mnkv+2;
   if((idfl==1)&&(dp.Ns+2>ns)) ns=dp.Ns+2;
   alloc_ws(ns);
   alloc_wsd(6);
   if (idfl)
   {
      alloc_wv(2*dpr.nkv+2);
      alloc_wvd(4);
   }
}



void set_dfl_modes(void)
{
   int idfl;
   int status[3];
   dflst_t dfl_status;
   double wt1,wt2,wdt;

   if (my_rank==0)
      printf("\n\n");

   idfl=0;
   while(1)
   {
      dfl_status=dfl_gen_parms(idfl).status;
      if(dfl_status==DFL_OUTOFRANGE) break;
      if(dfl_status==DFL_DEF)
      {
         MPI_Barrier(MPI_COMM_WORLD);
         wt1=MPI_Wtime();
         
         dfl_modes2(idfl,status);

         wt2=MPI_Wtime();
         wdt=wt2-wt1;
         
         error_root((status[1]<0)||((status[1]==0)&&(status[0]<0)),1,
                    "set_dfl_modes [cmp_solvers.c]","Generation of deflation subspace %d "
                    "failed (status = %d;%d)",idfl,status[0],status[1]);

         if (my_rank==0)
            printf("Deflation subspace %d generated in %.2e sec\n",idfl,wdt);
      }
      idfl++;
   }
}



void apply_solver(int isp,spinor_dble *in,double *rho,int *stat,double *wdt)
{
   spinor_dble *out;
   solver_parms_t sp;
   dirac_parms_t dp;
   double wt1,wt2;

   out=(reserve_wsd(1))[0];
   
   sp=solver_parms(isp);

   dp=qlat_parms(ifl[isp]);
   set_dirac_parms1(&dp);
   
   set_tm_parms(eoflg[isp]);

   if(sp.solver==CGNE)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      set_sd2zero(VOLUME,out);
      if(eoflg[isp]==0)
         (*rho)=tmcg(sp.nmx,sp.res,mu[isp],in,out,stat);
      else
         (*rho)=tmcgeo(sp.nmx,sp.res,mu[isp],in,out,stat);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      (*wdt)=wt2-wt1;
   }
   else if(sp.solver==SAP_GCR)
   {
      set_sd2zero(VOLUME,out);
      sap_gcr(sp.nkv,sp.nmx,sp.res,mu[isp],in,out,stat);

      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      set_sd2zero(VOLUME,out);
      (*rho)=sap_gcr(sp.nkv,sp.nmx,sp.res,mu[isp],in,out,stat);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      (*wdt)=wt2-wt1;
   }
   else if(sp.solver==DFL_SAP_GCR)
   {
      set_sd2zero(VOLUME,out);
      dfl_sap_gcr(sp.idfl,sp.nkv,sp.nmx,sp.res,mu[isp],in,out,stat);

      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      set_sd2zero(VOLUME,out);
      (*rho)=dfl_sap_gcr(sp.idfl,sp.nkv,sp.nmx,sp.res,mu[isp],in,out,stat);

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      (*wdt)=wt2-wt1;
   }
   else
   error(ifl==NULL,1,"apply_solver [cmp_solvers.c]",
         "UNKNOWN solver");   

   release_wsd();
}



int main(int argc,char *argv[])
{
   int nsize,opt;
   int icnfg,ncnfg,isp,isrc;
   int status[2],succ;
   double nrm,rho,wdt;
   int *avgstat[2],*nsucc;
   double *avgres,*avgwdt;
   spinor_dble *in,*sv,**sws;
   solver_parms_t sp;
   FILE *flog=NULL,*fin=NULL;
   int cnfg_type;
   char input_file[NAME_SIZE],log_file[NAME_SIZE];

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      opt=find_opt(argc,argv,"-i");
      if (opt!=0)
         strcpy(input_file,argv[opt+1]);
      else
         strcpy(input_file,"cmp_solvers.in");

      opt=find_opt(argc,argv,"-l");
      if (opt!=0)
         strcpy(log_file,argv[opt+1]);
      else
         strcpy(log_file,"cmp_solvers.log");

      flog=freopen(log_file,"w",stdout);
      fin=freopen(input_file,"r",stdin);

      printf("\n");
      printf("Compare performance of various solvers\n");
      printf("--------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
   }

   read_input();
   
   if (my_rank==0)
      fclose(fin);

   if (my_rank==0)
   {
      printf("Configurations %sn%d -> %sn%d in steps of %d\n",
             nbase,first,nbase,last,step);
      printf("Number of stocastic sources %d\n\n",nsrc);
      fflush(flog);
   }



   start_ranlux(0,1234);
   geometry();
   alloc_workspace();



   sws=reserve_wsd(2);
   in=sws[0];
   sv=sws[1];

   error_root(((last-first)%step)!=0,1,"main [cmp_solvers.c]",
              "last-first is not a multiple of step");

   nsize=name_size("%s/%sn%d",cnfg_dir,nbase,last);
   error_root(nsize>=NAME_SIZE,1,"main [cmp_solvers.c]",
              "cnfg_dir name is too long");



   avgstat[0]=malloc(nsp*sizeof(int));
   avgstat[1]=malloc(nsp*sizeof(int));
   nsucc=malloc(nsp*sizeof(int));
   avgres=malloc(nsp*sizeof(double));
   avgwdt=malloc(nsp*sizeof(double));

   for (isp=0;isp<nsp;isp++)
   {
      avgstat[0][isp]=0;
      avgstat[1][isp]=0;
      nsucc[isp]=0;
      avgres[isp]=0.0;
      avgwdt[isp]=0.0;
   }


   ncnfg=0;
   for (icnfg=first;icnfg<=last;icnfg+=step)
   {
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      cnfg_type=import_cnfg(cnfg_file);
      error_root(cnfg_type!=gauge(),1,"main [cmp_solvers.c]",
                 "Wrong gauge group for imported configuration");
      
      ncnfg+=1;

      set_dfl_modes();
      
      if (my_rank==0)
      {
         printf("\n\nicnfg     isp    isrc  |  st[0]  st[1]         res      time\n");
         fflush(flog);
      }

      for (isrc=0;isrc<nsrc;isrc++)
      {
         random_sd(VOLUME,in,1.0);
         assign_sd2sd(VOLUME,in,sv);

         for (isp=0;isp<nsp;isp++)
         {
            assign_sd2sd(VOLUME,sv,in);
            bnd_sd2zero(ALL_PTS,in);
            nrm=sqrt(norm_square_dble(VOLUME,1,in));

            apply_solver(isp,in,&rho,status,&wdt);
            
            sp=solver_parms(isp);
            if(sp.solver!=DFL_SAP_GCR) status[1]=0;

            if (my_rank==0)
            {
               printf("%5d   %5d   %5d  |  %5d  %5d    %.2e  %.2e\n",
                       icnfg,isp,isrc,status[0],status[1],rho/nrm,wdt);
               fflush(flog);
            }
               
            succ=0;
            if((status[0]>=0)&&(status[1]>=0))
               succ=1;
            
            if(succ)
            {
               avgstat[0][isp]+=status[0];
               avgstat[1][isp]+=status[1];
               nsucc[isp]++;
            }
            avgres[isp]+=rho/nrm;
            avgwdt[isp]+=wdt;
         }
      }
   }

   for (isp=0;isp<nsp;isp++)
   {
      if (nsucc[isp]!=0)
      {
         avgstat[0][isp]/=nsucc[isp];
         avgstat[1][isp]/=nsucc[isp];
      }
      avgres[isp]/=(ncnfg*nsrc);
      avgwdt[isp]/=(ncnfg*nsrc);
   }


   if (my_rank==0)
   {
      printf("\n");
      printf("Summary of results\n");
      printf("------------------\n\n");

      printf("Processed %d configurations\n",ncnfg);
      printf("Used %d random spinors per configuration\n\n",nsrc);
   
      for (isp=0;isp<nsp;isp++)
      {
         printf("Solver %d:\n",isp);
         printf("Success rate = %.2f\n",((double)nsucc[isp])/(ncnfg*nsrc));
         if (nsucc[isp]!=0)
         {
            printf("Average status (over successful attempts) = %d %d\n",avgstat[0][isp],avgstat[1][isp]);
         }
         printf("Average relative residue (over all attempts) = %.2e\n",avgres[isp]);
         printf("Average time (over all attempts, w/o preparatory steps) = %.2e sec\n",avgwdt[isp]);
         printf("\n");
      }

      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
