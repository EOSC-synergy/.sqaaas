
/*******************************************************************************
*
* File force_parms.c
*
* Copyright (C) 2011, 2012 Martin Luescher
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Force parameter data base
*
* The externally accessible functions are
*
*   force_parms_t set_force_parms(int ifr,force_t force,int ipf,int ifl,
*                                 int *irat,int *imu,int *isp,int *ncr)
*     Sets the parameters in the force parameter set number ifr and returns
*     a structure containing them (see the notes).
*
*   force_parms_t force_parms(int ifr)
*     Returns a structure containing the force parameter set number ifr
*     (see the notes).
*
*   void read_force_parms(int ifr)
*     On process 0, this program scans stdin for a line starting with the
*     string "[Force <int>]" (after any number of blanks), where <int> is
*     the integer value passed by the argument. An error occurs if no such
*     line or more than one is found. The lines
*
*       force   <force_t>
*       ipf     <int>
*       ifl     <int>
*       irat    <int> <int> <int>
*       imu     <int> [<int>]
*       isp     <int> [<int>]
*       ncr     <int> [<int>]
*
*     are then read using read_line() [utils/mutils.c]. Depending on the
*     value of "force", some lines are not read and can be omitted in the
*     input file. The number of integer items on the lines with tag "imu"
*     and "isp" and "ncr" depends on the force too. The data are then added
*     to the data base by calling set_force_parms(ifr,...).
*
*   void read_force_parms2(int ifr)
*     Same as read_force_parms() except that only the lines
*
*       force   <force_t>
*       isp     <int> [<int>]
*       ncr     <int> [<int>]
*
*     are read from stdin. All other force parameters are inferred from
*     the parameters of the action no ifr so that the force is the one
*     deriving from that action. An error occurs if the parameters of the
*     action no ifr have not previously been added to the data base or
*     if the force and action types do not match.
*
*   void print_force_parms(void)
*     Prints the parameters of the defined forces to stdout on MPI
*     process 0.
*
*   void print_force_parms2(void)
*     Prints the parameters of the defined forces to stdout on MPI
*     process 0 in a short format corresponding to read_force_parms2().
*
*   void write_force_parms(FILE *fdat)
*     Writes the parameters of the defined forces to the file fdat on
*     MPI process 0.
*
*   void check_force_parms(FILE *fdat)
*     Compares the parameters of the defined forces with those stored
*     on the file fdat on MPI process 0, assuming the latter were written
*     to the file by the program write_force_parms().
*
* Notes:
*
* For a description of the supported forces and their parameters see
* forces/README.forces.
*
* The elements of a structure of type force_parms_t are
*
*   force   Force program used. This parameter is an enum type with
*           one of the following values:
*
*            FRG_SU3         (program force0() [forces/force0.c]),
*
*            FRF_TM1         (program force1() [forces/force1.c]),
*
*            FRF_TM1_EO      (program force4() [forces/force4.c]),
*
*            FRF_TM1_EO_SDET (program force4() [forces/force4.c]),
*
*            FRF_TM2         (program force2() [forces/force2.c]),
*
*            FRF_TM2_EO      (program force5() [forces/force5.c]),
*
*            FRF_RAT         (program force3() [forces/force3.c]),
*
*            FRF_RAT_SDET    (program force3() [forces/force3.c]),
*
*            FRG_U1          (program force6() [forces/force6.c]),
*
*   ipf     Pseudo-fermion field index (see mdflds/mdflds.c),
*
*   ifl     Index of quark flavour in parameter data base
*           (see flags/lat_parms.c),
*
*   irat    Indices specifying a rational function (see ratfcts/ratfcts.c),
*
*   imu     Twisted mass indices (see flags/hmc_parms.c),
*
*   isp     Solver parameter set indices (see flags/solver_parms.c),
*
*   ncr     Chronological solver stack sizes (see update/chrono.c),
*
*   icr     Chronological solver stack indices (set internally).
*
* Depending on the force, some parameters are not used and are set to zero
* by set_force_parms() independently of the values of the arguments. In
* particular, for a given force, only the required number of integers are
* read from the arrays imu, isp and ncr passed to the program.
*
* The number of twisted mass indices is 1 and 2 in the case of the forces
* FRF_TM1* and FRF_TM2*, respectively. These forces require a chronological
* solver stack size to be specified and 1 solver parameter set to be used
* for the solution of the Dirac equation with twisted mass index imu[0].
*
* Up to 32 force parameter sets, labeled by an index ifr=0,1,..,31, can
* be specified. Once a set is specified, it cannot be changed by calling
* set_force_parms() again. Force parameters must be globally the same.
*
* Except for force_parms(), the programs in this module perform global
* operations and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define FORCE_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"

#define IFRMAX 128

static int init=0,icr=0;
static force_parms_t fp[IFRMAX+1]={{FORCES,0,0,{0,0,0},{0,0,0,0},{0,0,0,0},
                                    {0,0,0,0},{0,0,0,0}}};


static void init_fp(void)
{
   int i;

   for (i=1;i<=IFRMAX;i++)
      fp[i]=fp[0];

   init=1;
}


force_parms_t set_force_parms(int ifr,force_t force,int ipf,int ifl,
                              int *irat,int *imu,int *isp,int *ncr)
{
   int i,ie;
   int rat[3],mu[4],sp[4],nc[4],ic[4];

   if (init==0)
      init_fp();

   for (i=0;i<3;i++)
      rat[i]=0;

   for (i=0;i<4;i++)
   {
      mu[i]=0;
      sp[i]=0;
      nc[i]=0;
      ic[i]=0;
   }

   if ((force==FRG_SU3)||(force==FRG_U1)||(force==FORCES))
   {
      ipf=0;
      ifl=0;
   }
   else if ((force==FRF_TM1)||(force==FRF_TM1_EO)||(force==FRF_TM1_EO_SDET))
   {
      mu[0]=imu[0];
      sp[0]=isp[0];

      if (ncr[0]>0)
      {
         icr+=1;
         nc[0]=ncr[0];
         ic[0]=icr;
      }
   }
   else if ((force==FRF_TM2)||(force==FRF_TM2_EO))
   {
      mu[0]=imu[0];
      mu[1]=imu[1];
      sp[0]=isp[0];

      if (ncr[0]>0)
      {
         icr+=1;
         nc[0]=ncr[0];
         ic[0]=icr;
      }
   }
   else if ((force==FRF_RAT)||(force==FRF_RAT_SDET))
   {
      rat[0]=irat[0];
      rat[1]=irat[1];
      rat[2]=irat[2];
      sp[0]=isp[0];
   }

   check_global_int("set_force_parms",23,
                    ifr,(int)(force),ipf,ifl,
                    rat[0],rat[1],rat[2],
                    mu[0],mu[1],mu[2],mu[3],
                    sp[0],sp[1],sp[2],sp[3],
                    nc[0],nc[1],nc[2],nc[3],
                    ic[0],ic[1],ic[2],ic[3]);

   ie=0;
   ie|=((ifr<0)||(ifr>=IFRMAX));
   ie|=(force==FORCES);
   ie|=((ipf<0)||(ifl<0));

   for (i=0;i<3;i++)
      ie|=(rat[i]<0);

   for (i=0;i<4;i++)
   {
      ie|=(mu[i]<0);
      ie|=(sp[i]<0);
      ie|=(nc[i]<0);
   }

   error_root(ie!=0,1,"set_force_parms [force_parms.c]",
              "Parameters are out of range");

   error_root(fp[ifr].force!=FORCES,1,"set_force_parms [force_parms.c]",
              "Attempt to reset already specified force parameters");

   fp[ifr].force=force;
   fp[ifr].ipf=ipf;
   fp[ifr].ifl=ifl;

   for (i=0;i<3;i++)
      fp[ifr].irat[i]=rat[i];

   for (i=0;i<4;i++)
   {
      fp[ifr].imu[i]=mu[i];
      fp[ifr].isp[i]=sp[i];
      fp[ifr].ncr[i]=nc[i];
      fp[ifr].icr[i]=ic[i];
   }

   return fp[ifr];
}


force_parms_t force_parms(int ifr)
{
   if (init==0)
      init_fp();

   if ((ifr>=0)&&(ifr<IFRMAX))
      return fp[ifr];
   else
   {
      error_loc(1,1,"force_parms [force_parms.c]",
                "Force index is out of range");
      return fp[IFRMAX];
   }
}


void read_force_parms(int ifr)
{
   int my_rank,i,force;
   int ipf,ifl,irat[3],imu[4],isp[4],ncr[4];
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   force=-1;
   ipf=0;
   ifl=0;

   for (i=0;i<3;i++)
      irat[i]=0;

   for (i=0;i<4;i++)
   {
      imu[i]=0;
      isp[i]=0;
      ncr[i]=0;
   }

   if (my_rank==0)
   {
      sprintf(line,"Force %d",ifr);
      find_section(line);
      read_line("force","%s",line);

      if (strcmp(line,"FRG_SU3")==0)
         force=FRG_SU3;
      else if (strcmp(line,"FRG_U1")==0)
         force=FRG_U1;
      else if (strcmp(line,"FRF_TM1")==0)
      {
         force=FRF_TM1;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("imu","%d",imu);
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (strcmp(line,"FRF_TM1_EO")==0)
      {
         force=FRF_TM1_EO;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("imu","%d",imu);
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (strcmp(line,"FRF_TM1_EO_SDET")==0)
      {
         force=FRF_TM1_EO_SDET;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("imu","%d",imu);
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (strcmp(line,"FRF_TM2")==0)
      {
         force=FRF_TM2;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("imu","%d %d",imu,imu+1);
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (strcmp(line,"FRF_TM2_EO")==0)
      {
         force=FRF_TM2_EO;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("imu","%d %d",imu,imu+1);
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (strcmp(line,"FRF_RAT")==0)
      {
         force=FRF_RAT;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("irat","%d %d %d",irat,irat+1,irat+2);
         read_line("isp","%d",isp);
      }
      else if (strcmp(line,"FRF_RAT_SDET")==0)
      {
         force=FRF_RAT_SDET;
         read_line("ipf","%d",&ipf);
         read_line("ifl","%d",&ifl);
         read_line("irat","%d %d %d",irat,irat+1,irat+2);
         read_line("isp","%d",isp);
      }
      else
         error_root(1,1,"read_force_parms [force_parms.c]",
                    "Unknown force %s",line);
   }

   if (NPROC>1)
   {
      MPI_Bcast(&force,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ipf,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ifl,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(irat,3,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(imu,4,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(isp,4,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(ncr,4,MPI_INT,0,MPI_COMM_WORLD);
   }

   set_force_parms(ifr,force,ipf,ifl,irat,imu,isp,ncr);
}


void read_force_parms2(int ifr)
{
   int my_rank,i,ie,force;
   int ipf,ifl,irat[3],imu[4],isp[4],ncr[4];
   char line[NAME_SIZE];
   action_parms_t ap;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   ie=0;
   force=-1;
   ipf=0;
   ifl=0;

   for (i=0;i<3;i++)
      irat[i]=0;

   for (i=0;i<4;i++)
   {
      imu[i]=0;
      isp[i]=0;
      ncr[i]=0;
   }

   if (my_rank==0)
   {
      ap=action_parms(ifr);
      error_root(ap.action==ACTIONS,1,"read_force_parms2 [force_parms.c]",
                 "Undefined action");

      sprintf(line,"Force %d",ifr);
      find_section(line);
      read_line("force","%s",line);

      if (ap.action==ACG_SU3)
      {
         ie=strcmp(line,"FRG_SU3");
         force=FRG_SU3;
      }
      else if (ap.action==ACG_U1)
      {
         ie=strcmp(line,"FRG_U1");
         force=FRG_U1;
      }
      else if (ap.action==ACF_TM1)
      {
         ie=strcmp(line,"FRF_TM1");
         force=FRF_TM1;
         ipf=ap.ipf;
         ifl=ap.ifl;
         imu[0]=ap.imu[0];
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (ap.action==ACF_TM1_EO)
      {
         ie=strcmp(line,"FRF_TM1_EO");
         force=FRF_TM1_EO;
         ipf=ap.ipf;
         ifl=ap.ifl;
         imu[0]=ap.imu[0];
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (ap.action==ACF_TM1_EO_SDET)
      {
         ie=strcmp(line,"FRF_TM1_EO_SDET");
         force=FRF_TM1_EO_SDET;
         ipf=ap.ipf;
         ifl=ap.ifl;
         imu[0]=ap.imu[0];
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (ap.action==ACF_TM2)
      {
         ie=strcmp(line,"FRF_TM2");
         force=FRF_TM2;
         ipf=ap.ipf;
         ifl=ap.ifl;
         imu[0]=ap.imu[0];
         imu[1]=ap.imu[1];
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (ap.action==ACF_TM2_EO)
      {
         ie=strcmp(line,"FRF_TM2_EO");
         force=FRF_TM2_EO;
         ipf=ap.ipf;
         ifl=ap.ifl;
         imu[0]=ap.imu[0];
         imu[1]=ap.imu[1];
         read_line("isp","%d",isp);
         read_line("ncr","%d",ncr);
      }
      else if (ap.action==ACF_RAT)
      {
         ie=strcmp(line,"FRF_RAT");
         force=FRF_RAT;
         ipf=ap.ipf;
         ifl=ap.ifl;
         irat[0]=ap.irat[0];
         irat[1]=ap.irat[1];
         irat[2]=ap.irat[2];
         read_line("isp","%d",isp);
      }
      else if (ap.action==ACF_RAT_SDET)
      {
         ie=strcmp(line,"FRF_RAT_SDET");
         force=FRF_RAT_SDET;
         ipf=ap.ipf;
         ifl=ap.ifl;
         irat[0]=ap.irat[0];
         irat[1]=ap.irat[1];
         irat[2]=ap.irat[2];
         read_line("isp","%d",isp);
      }
      else
         error_root(1,1,"read_force_parms2 [force_parms.c]",
                    "Unknown action");

      error_root(ie!=0,1,"read_force_parms2 [force_parms.c]",
                    "Force and action types do not match");
   }

   if (NPROC>1)
   {
      MPI_Bcast(&force,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ipf,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ifl,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(irat,3,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(imu,4,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(isp,4,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(ncr,4,MPI_INT,0,MPI_COMM_WORLD);
   }

   set_force_parms(ifr,force,ipf,ifl,irat,imu,isp,ncr);
}


void print_force_parms(void)
{
   int my_rank,i;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if ((my_rank==0)&&(init==1))
   {
      for (i=0;i<IFRMAX;i++)
      {
         if (fp[i].force!=FORCES)
         {
            printf("Force %d:\n",i);

            if (fp[i].force==FRG_SU3)
               printf("FRG_SU3 force\n\n");
            else if (fp[i].force==FRG_U1)
               printf("FRG_U1 force\n\n");
            else if (fp[i].force==FRF_TM1)
            {
               printf("FRF_TM1 force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("imu = %d\n",fp[i].imu[0]);
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM1_EO)
            {
               printf("FRF_TM1_EO force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("imu = %d\n",fp[i].imu[0]);
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM1_EO_SDET)
            {
               printf("FRF_TM1_EO_SDET force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("imu = %d\n",fp[i].imu[0]);
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM2)
            {
               printf("FRF_TM2 force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("imu = %d %d\n",fp[i].imu[0],fp[i].imu[1]);
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM2_EO)
            {
               printf("FRF_TM2_EO force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("imu = %d %d\n",fp[i].imu[0],fp[i].imu[1]);
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_RAT)
            {
               printf("FRF_RAT force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("irat = %d %d %d\n",
                      fp[i].irat[0],fp[i].irat[1],fp[i].irat[2]);
               printf("isp = %d\n\n",fp[i].isp[0]);
            }
            else if (fp[i].force==FRF_RAT_SDET)
            {
               printf("FRF_RAT_SDET force\n");
               printf("ipf = %d\n",fp[i].ipf);
               printf("ifl = %d\n",fp[i].ifl);
               printf("irat = %d %d %d\n",
                      fp[i].irat[0],fp[i].irat[1],fp[i].irat[2]);
               printf("isp = %d\n\n",fp[i].isp[0]);
            }
            else
               printf("UNKNOWN force\n\n");
         }
      }
   }
}


void print_force_parms2(void)
{
   int my_rank,i;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if ((my_rank==0)&&(init==1))
   {
      for (i=0;i<IFRMAX;i++)
      {
         if (fp[i].force!=FORCES)
         {
            printf("Force %d:\n",i);

            if (fp[i].force==FRG_SU3)
               printf("FRG_SU3 force\n\n");
            else if (fp[i].force==FRG_U1)
               printf("FRG_U1 force\n\n");
            else if (fp[i].force==FRF_TM1)
            {
               printf("FRF_TM1 force\n");
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM1_EO)
            {
               printf("FRF_TM1_EO force\n");
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM1_EO_SDET)
            {
               printf("FRF_TM1_EO_SDET force\n");
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM2)
            {
               printf("FRF_TM2 force\n");
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_TM2_EO)
            {
               printf("FRF_TM2_EO force\n");
               printf("isp = %d\n",fp[i].isp[0]);
               printf("ncr = %d\n\n",fp[i].ncr[0]);
            }
            else if (fp[i].force==FRF_RAT)
            {
               printf("FRF_RAT force\n");
               printf("isp = %d\n\n",fp[i].isp[0]);
            }
            else if (fp[i].force==FRF_RAT_SDET)
            {
               printf("FRF_RAT_SDET force\n");
               printf("isp = %d\n\n",fp[i].isp[0]);
            }
            else
               printf("UNKNOWN force\n\n");
         }
      }
   }
}


void write_force_parms(FILE *fdat)
{
   int i;

   if (init==1)
   {
      for (i=0;i<IFRMAX;i++)
      {
         if (fp[i].force!=FORCES)
         {
            write_little_int(1,fdat,23,
                        i,fp[i].force,fp[i].ipf,fp[i].ifl,
                        fp[i].irat[0],fp[i].irat[1],fp[i].irat[2],
                        fp[i].imu[0],fp[i].imu[1],fp[i].imu[2],fp[i].imu[3],
                        fp[i].isp[0],fp[i].isp[1],fp[i].isp[2],fp[i].isp[3],
                        fp[i].ncr[0],fp[i].ncr[1],fp[i].ncr[2],fp[i].ncr[3],
                        fp[i].icr[0],fp[i].icr[1],fp[i].icr[2],fp[i].icr[3]);
         }
      }
   }
}


void check_force_parms(FILE *fdat)
{
   int i;

   if (init==1)
   {
      for (i=0;i<IFRMAX;i++)
      {
         if (fp[i].force!=FORCES)
         {
            check_little_int("check_force_parms",fdat,23,
                        i,fp[i].force,fp[i].ipf,fp[i].ifl,
                        fp[i].irat[0],fp[i].irat[1],fp[i].irat[2],
                        fp[i].imu[0],fp[i].imu[1],fp[i].imu[2],fp[i].imu[3],
                        fp[i].isp[0],fp[i].isp[1],fp[i].isp[2],fp[i].isp[3],
                        fp[i].ncr[0],fp[i].ncr[1],fp[i].ncr[2],fp[i].ncr[3],
                        fp[i].icr[0],fp[i].icr[1],fp[i].icr[2],fp[i].icr[3]);
         }
      }
   }
}
