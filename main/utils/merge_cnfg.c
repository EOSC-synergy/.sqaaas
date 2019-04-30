
/*******************************************************************************
*
* File merge_cnfg.c
*
* Copyright (C) 2018 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Merge QCD and QED configurations.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "global.h"



int main(int argc,char *argv[])
{
   int my_rank,bc,cs,opt;
   int first,last;
   char cnfg_dir[NAME_SIZE],nbase[NAME_SIZE];
   int qcd_first,qcd_last,qcd_step;
   char qcd_cnfg_dir[NAME_SIZE],qcd_nbase[NAME_SIZE];
   int qed_first,qed_last,qed_step;
   char qed_cnfg_dir[NAME_SIZE],qed_nbase[NAME_SIZE];
   double su3phi[2],su3phi_prime[2],u1phi,u1phi_prime;
   char cnfg_file[NAME_SIZE],input_file[NAME_SIZE],log_file[NAME_SIZE];
   int nsize,icnfg;
   FILE *flog=NULL,*fin=NULL;
   int cnfg_type;
   char line[NAME_SIZE];

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      opt=find_opt(argc,argv,"-i");
      if (opt!=0)
         strcpy(input_file,argv[opt+1]);
      else
         strcpy(input_file,"merge_cnfg.in");

      opt=find_opt(argc,argv,"-l");
      if (opt!=0)
         strcpy(log_file,argv[opt+1]);
      else
         strcpy(log_file,"merge_cnfg.log");

      flog=freopen(log_file,"w",stdout);
      fin=freopen(input_file,"r",stdin);

      printf("\n");
      printf("Merge QCD and QED configurations\n");
      printf("--------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

      find_section("QCD+QED configurations");
      read_line("name","%s",nbase);
      read_line("cnfg_dir","%s",cnfg_dir);

      find_section("QCD configurations");
      read_line("name","%s",qcd_nbase);
      read_line("cnfg_dir","%s",qcd_cnfg_dir);
      read_line("first","%d",&qcd_first);
      read_line("last","%d",&qcd_last);
      read_line("step","%d",&qcd_step);

      find_section("QED configurations");
      read_line("name","%s",qed_nbase);
      read_line("cnfg_dir","%s",qed_cnfg_dir);
      read_line("first","%d",&qed_first);
      read_line("step","%d",&qed_step);

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

      su3phi[0]=0.0;
      su3phi[1]=0.0;
      su3phi_prime[0]=0.0;
      su3phi_prime[1]=0.0;
      u1phi=0.0;
      u1phi_prime=0.0;

      if (bc==1)
      {
         read_dprms("su3phi",2,su3phi);
         read_dprms("u1phi",1,&u1phi);
      }
      if ((bc==1)||(bc==2))
      {
         read_dprms("su3phi'",2,su3phi_prime);
         read_dprms("u1phi'",1,&u1phi_prime);
      }

      fclose(fin);
   }

   MPI_Bcast(qcd_nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(qcd_cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&qcd_first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&qcd_last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&qcd_step,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(qed_nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(qed_cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&qed_first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&qed_step,1,MPI_INT,0,MPI_COMM_WORLD);
   qed_last=qed_first+qed_step*((qcd_last-qcd_first)/qcd_step);

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(su3phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(su3phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&u1phi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&u1phi_prime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   set_flds_parms(3,0);
   print_flds_parms();

   set_bc_parms(bc,cs,su3phi,su3phi_prime,u1phi,u1phi_prime);
   print_bc_parms();

   start_ranlux(0,1234);
   geometry();

   error_root(((qcd_last-qcd_first)%qcd_step)!=0,1,"main [merge_cnfg.c]",
              "qcd_last-qcd_first is not a multiple of step");

   nsize=name_size("%s/%sn%d",qcd_cnfg_dir,qcd_nbase,qcd_last);
   error_root(nsize>=NAME_SIZE,1,"main [merge_cnfg.c]",
              "qcd_cnfg_dir name is too long");

   nsize=name_size("%s/%sn%d",qed_cnfg_dir,qed_nbase,qed_last);
   error_root(nsize>=NAME_SIZE,1,"main [merge_cnfg.c]",
              "qed_cnfg_dir name is too long");

   first=1;
   last=1+(qcd_last-qcd_first)/qcd_step;
   
   for (icnfg=first;icnfg<=last;icnfg++)
   {
      sprintf(cnfg_file,"%s/%sn%d",qcd_cnfg_dir,qcd_nbase,
              qcd_first+(icnfg-1)*qcd_step);
              
      if (my_rank==0)
      {
         printf("Read  %s\n",cnfg_file);
         fflush(flog);
      }
      cnfg_type=import_cnfg(cnfg_file);
      error_root(cnfg_type!=1,1,"main [merge_cnfg.c]",
                 "Wrong gauge group for imported configuration -- expected SU(3)");

      sprintf(cnfg_file,"%s/%sn%d",qed_cnfg_dir,qed_nbase,
              qed_first+(icnfg-1)*qed_step);

      if (my_rank==0)
      {
         printf("Read  %s\n",cnfg_file);
         fflush(flog);
      }
      cnfg_type=import_cnfg(cnfg_file);
      error_root(cnfg_type!=2,1,"main [merge_cnfg.c]",
                 "Wrong gauge group for imported configuration -- expected U(1)");
      
      sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,icnfg);
      if (my_rank==0)
      {
         printf("Write %s\n\n",cnfg_file);
         fflush(flog);
      }
      export_cnfg(cnfg_file);
   }

   if (my_rank==0)
      fclose(flog);

   MPI_Finalize();
   exit(0);
}
