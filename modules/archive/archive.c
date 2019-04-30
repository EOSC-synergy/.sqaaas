
/*******************************************************************************
*
* File archive.c
*
* Copyright (C) 2017 Marina Marinkovic, Agostino Patella
*
* Based on openQCD-1.6/modules/archive/archive.c
* Copyright (C) 2005, 2007, 2009-2014 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Programs to read and write SU(3), U(1) and SU(3)xU(1) gauge-field
* configurations.
*
* The externally accessible functions are
*
*   void write_cnfg(char *out)
*     Writes the lattice sizes, the process grid sizes, the coordinates of the
*     calling process, the state of the random number generators, the local
*     plaquette sums and the local double-precision (SU(3)xU(1) or U(1)/SU(3)
*     only) gauge field to the file "out".
*
*   int read_cnfg(char *in)
*     Detects the type of gauge configuration written in the file "in"
*     (i.e. (SU(3)xU(1) or U(1)/SU(3) only). It reads only the gauge fields
*     that are active, as returned by the gauge() program defined in
*     'modules/flags/lat_parms.c'. The program then resets the random number
*     generator and checks that the restored field is compatible with the chosen
*     boundary conditions (see the notes). It returns 1 if "in" contains only
*     the SU(3) field, 2 if "in" contains only the U(1) field, and 3 if "in"
*     contains both. The program assumes that the file "in" was written by the
*     program write_cnfg().
*
*   void export_cnfg(char *out)
*     Writes the lattice sizes and the global double-precision U(1) and/or SU(3)
*     gauge field to the file "out" from process 0 in the universal format
*     specified below (see the notes).
*
*   int import_cnfg(char *in)
*     Detects the type of gauge configuration written in the file "in"
*     (i.e. (SU(3)xU(1) or U(1)/SU(3) only). It reads only the gauge fields
*     that are active, as returned by the gauge() program defined in
*     'modules/flags/lat_parms.c'. The program then resets the random number
*     generator and checks that the restored field is compatible with the chosen
*     boundary conditions. The U(1) and SU(3) fields are periodically extended
*     if needed (see the notes). It returns 1 if "in" contains only the SU(3)
*     field, 2 if "in" contains only the U(1) field, and 3 if "in" contains
*     both. The program assumes that the file "in" was written by the program
*     export_cnfg() in the universal format, and file "in" is accessed by
*     process 0 only.
*
* Notes:
*
* The program export_cnfg() first writes the lattice sizes, the average
* SU(3) plaquette and the SU(3) link variables if the SU(3) gauge field is
* active, the average U(1) plaquette and non-compact U(1) link variables if the
* U(1) gauge field is active. The order of writing the fields is the following:
* link variables in the directions +0,-0,...,+3,-3 at the first odd point, the
* second odd point, and so on. The order of the point (x0,x1,x2,x3) with
* Cartesian coordinates in the range 0<=x0<N0,...,0<=x3<N3
* is determined by the index
*
*   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0,
*
* where N0,N1,N2,N3 are the global lattice sizes (N0=NPROC0*L0, etc.). The
* average U(1) field and SU(3) plaquettes are calculated by summing the
* plaquette values over all points in the lattice, including the space-like
* plaquettes at time N0 if SF or open-SF boundary conditions are chosen, and
* dividing the sum by 6*N0*N1*N2*N3.
*
* Independently of the machine, the export function writes the data to the
* output file in little-endian byte order. Integers and double-precision
* numbers on the output file occupy 4 and 8 bytes, respectively, the latter
* being formatted according to the IEEE-754 standard. The import function
* assumes the data on the input file to be little endian and converts them
* to big-endian order if the machine is big endian. Exported configurations
* can thus be safely exchanged between different machines.
*
* Let N0,N1,N2,N3 be the current lattice sizes, and let n0,n1,n2,n3 be the
* lattice sizes read from the configuration file. N0 must be an integer multiple
* of n0 if the temporal direction has periodic boundary conditions, and N0 must
* be equal to n0 otherwise. N1 must be an integer multiple of n1 if the spatial
* directions have periodic boundary conditions, and N1 must be an even multiple
* of n1 in case of C* boundary conditions. Nk must be an integer multiple of nk
* for k=2,3. The imported configuration is extended by periodicity on the whole
* lattice. In case of C* bundary conditions, the imported configuration is
* extended by periodicity on the fundamental domain and extended by means of the
* orbifold constraint to the mirror domain.
*
* Compatibility of a configuration with the chosen boundary conditions is
* established by calling check_bc() [lattice/bcnds.c] for SU(3) fields and
* check_ad_bc() [u1flds/ad_bcnds.c], with a tolerance on the boundary link
* variables of 64.0*DBL_EPSILON, and by checking that the average SU(3)
* plaquette and average U(1) field coincide with the value read from the
* configuration file. On exit both read_cnfg() and import_cnfg() set the
* boundary values of the field (if any) to the ones stored in the parameter
* data base so as to guarantee that they are bit-identical to the latter.
*
* While both active and inactive (i.e. as given by the gauge() program) gauge
* fields are read, checks are performed only for active gauge fields.
*
* All programs in this module may involve global communications and must be
* called simultaneously on all processes.
*
*******************************************************************************/

#define ARCHIVE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "u1flds.h"
#include "archive.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int endian,ns,nd,*state=NULL,*statebuf=NULL;
static su3_dble *ubuf=NULL,*vbuf,*udb;
static double *abuf=NULL,*vabuf,*adb;

static void alloc_state(void)
{
   int n;
   
   ns=rlxs_size();
   nd=rlxd_size();
   
   if (ns<nd)
      n=nd;
   else
      n=ns;
   
   state=malloc(n*sizeof(int));
   error(state==NULL,1,"alloc_state [archive.c]",
         "Unable to allocate auxiliary array");
}

static void alloc_statebuf(void)
{
   int n;

   ns=rlxs_size();
   nd=rlxd_size();
   n=ns+nd;

   statebuf=malloc(n*sizeof(int));
   error(statebuf==NULL,1,"alloc_state [archive.c]",
         "Unable to allocate auxiliary array for the rlx state");
}

static void check_machine(void)
{
   error_root(sizeof(stdint_t)!=4,1,"check_machine [archive.c]",
              "Size of a stdint_t integer is not 4");
   error_root(sizeof(double)!=8,1,"check_machine [archive.c]",
              "Size of a double is not 8");

   endian=endianness();
   error_root(endian==UNKNOWN_ENDIAN,1,"check_machine [archive.c]",
              "Unkown endianness");
}

static void alloc_ubuf(int my_rank)
{
   if (my_rank==0)
   {
      ubuf=amalloc(4*(L3+N3)*sizeof(su3_dble),ALIGN);
      vbuf=ubuf+4*L3;
   }
   else
   {
      ubuf=amalloc(4*L3*sizeof(su3_dble),ALIGN);
      vbuf=NULL;
   }

   error(ubuf==NULL,1,"alloc_ubuf [archive.c]",
         "Unable to allocate auxiliary array for SU(3) gauge fields");
}

static void alloc_abuf(int my_rank)
{
   if (my_rank==0)
   {
      abuf=amalloc(4*(L3+N3)*sizeof(double),ALIGN);
      vabuf=abuf+4*L3;
   }
   else
   {
      abuf=amalloc(4*L3*sizeof(double),ALIGN);
      vabuf=NULL;
   }

   error(abuf==NULL,1,"alloc_abuf [archive.c]",
         "Unable to allocate auxiliary array for U(1) gauge fields");
}

static int which_cnfg_type(char *in,int loc)
{
   int ldat[16],ie,cnfg_type=0;
   long long int ir=0,size1=0,size2=0,size3=0,globvol,np[4],d;
   int my_rank;
   stdint_t lsize[4];
   FILE *fin;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   if (ubuf==NULL)
      alloc_ubuf(my_rank);

   if (loc==1)
   {
      if (statebuf==NULL)
         alloc_statebuf();
      
      fin=fopen(in,"rb");
      error_loc(fin==NULL,1,"which_cnfg_type [archive.c]",
                "Unable to open configuration file from local disks");
      
      ir=fread(ldat,sizeof(int),16,fin);

      ie=0;
      ie|=((ldat[0]!=NPROC0)||(ldat[1]!=NPROC1)||
           (ldat[2]!=NPROC2)||(ldat[3]!=NPROC3));
      ie|=((ldat[4]!=L0)||(ldat[5]!=L1)||
           (ldat[6]!=L2)||(ldat[7]!=L3));
      ie|=((ldat[8]!=NPROC0_BLK)||(ldat[9]!=NPROC1_BLK)||
           (ldat[10]!=NPROC2_BLK)||(ldat[11]!=NPROC3_BLK));
      ie|=((ldat[12]!=cpr[0])||(ldat[13]!=cpr[1])||
           (ldat[14]!=cpr[2])||(ldat[15]!=cpr[3]));
      error(ie!=0,1,"which_cnfg_type [archive.c]","Unexpected lattice data");
      
      ir+=fread(statebuf,sizeof(int),ns+nd,fin);
      error_root(ir!=(ns+nd+16),1,"which_cnfg_type [archive.c]","Incorrect read count (1)");

      size1=(18*4*VOLUME+1);
      size2=(4*VOLUME+1);
      size3=(19*4*VOLUME+2);

      ir=0;
      while((d=fread(ubuf,sizeof(double),4*18*L3,fin))==4*18*L3)
         ir+=d;
      ir+=d;

      if (ir==size3)
         cnfg_type=3;
      else if (ir==size2)
         cnfg_type=2;
      else if (ir==size1)
         cnfg_type=1;
      else
         error_root(1,1,"which_cnfg_type [archive.c]","Incorrect read count (2)");

      fclose(fin);
      
      check_global_int("which_cnfg_type [archive.c]",1,cnfg_type);
   }
   else if (loc==0)
   {
      if (my_rank==0)
      {
         fin=fopen(in,"rb");
         error_loc(fin==NULL,1,"which_cnfg_type [archive.c]",
                   "Unable to open exported configuration file");
      
         ir=fread(lsize,sizeof(stdint_t),4,fin);
         error_root(ir!=4,1,"which_cnfg_type [archive.c]","Incorrect read count (3)");

         check_machine();
         if (endian==BIG_ENDIAN)
            bswap_int(4,lsize);

         np[0]=(int)(lsize[0]);
         np[1]=(int)(lsize[1]);
         np[2]=(int)(lsize[2]);
         np[3]=(int)(lsize[3]);

         error_root((np[0]<1)||((N0%np[0])!=0)||
                    (np[1]<1)||((N1%np[1])!=0)||
                    (np[2]<1)||((N2%np[2])!=0)||
                    (np[3]<1)||((N3%np[3])!=0),1,"which_cnfg_type [archive.c]",
                    "Unexpected or incompatible lattice sizes");

         error_root((np[0]!=N0)&&(bc_type()!=3),1,"which_cnfg_type [archive.c]",
                    "Periodic extension in time is only possible when\n"
                    "periodic boundary conditions are chosen");

         globvol=np[0]*np[1]*np[2]*np[3];
         size1=(18*4*globvol+1);
         size2=(4*globvol+1);
         size3=(19*4*globvol+2);

         ir=0;
         while((d=fread(ubuf,sizeof(double),4*18*(N3+L3),fin))==4*18*(N3+L3))
            ir+=d;
         ir+=d;

         if (ir==size3)
            cnfg_type=3;
         else if (ir==size2)
            cnfg_type=2;
         else if (ir==size1)
            cnfg_type=1;
         else
            error_root(1,1,"which_cnfg_type [archive.c]","Incorrect read count (4)");

         fclose(fin);
      }
      
      MPI_Bcast(&cnfg_type,1,MPI_INT,0,MPI_COMM_WORLD);
   }
   else
      error(1,1,"which_cnfg_type [archive.c]","Second argument out of range");

   return cnfg_type;
}


void write_cnfg(char *out)
{
   int ldat[16],iw,cnfg_type;
   double su3plaq=0.0,u1plaq=0.0;
   FILE *fout;

   fout=fopen(out,"wb");
   error_loc(fout==NULL,1,"write_cnfg [archive.c]",
             "Unable to open output file");

   if (state==NULL)
      alloc_state();

   cnfg_type=gauge();

   if ((cnfg_type&1)!=0)
   {
      udb=udfld();
      su3plaq=plaq_sum_dble(0);
   }
   if ((cnfg_type&2)!=0)
   {
      adb=adfld();
      u1plaq=u1_plaq_sum_dble(0);
   }

   ldat[0]=NPROC0;
   ldat[1]=NPROC1;
   ldat[2]=NPROC2;
   ldat[3]=NPROC3;

   ldat[4]=L0;
   ldat[5]=L1;
   ldat[6]=L2;
   ldat[7]=L3;

   ldat[8]=NPROC0_BLK;
   ldat[9]=NPROC1_BLK;
   ldat[10]=NPROC2_BLK;
   ldat[11]=NPROC3_BLK;

   ldat[12]=cpr[0];
   ldat[13]=cpr[1];
   ldat[14]=cpr[2];
   ldat[15]=cpr[3];

   iw=fwrite(ldat,sizeof(int),16,fout);
   rlxs_get(state);
   iw+=fwrite(state,sizeof(int),ns,fout);
   rlxd_get(state);
   iw+=fwrite(state,sizeof(int),nd,fout);
   if ((cnfg_type&1)!=0)
   {
      iw+=fwrite(&su3plaq,sizeof(double),1,fout);
      iw+=fwrite(udb,sizeof(su3_dble),4*VOLUME,fout);
   }
   if ((cnfg_type&2)!=0)
   {
      iw+=fwrite(&u1plaq,sizeof(double),1,fout);
      iw+=fwrite(adb,sizeof(double),4*VOLUME,fout);
   }

   if (cnfg_type!=3)
      error_loc(iw!=(17+ns+nd+4*VOLUME),1,"write_cnfg [archive.c]",
                "Incorrect write count");
   else
      error_loc(iw!=(18+ns+nd+8*VOLUME),1,"write_cnfg [archive.c]",
                "Incorrect write count");

   fclose(fout);
}


int read_cnfg(char *in)
{
   int ldat[16],ir,ie,cnfg_type;
   double nplaq,su3plaq0,su3plaq1,u1plaq0,u1plaq1,eps;
   FILE *fin;

   if (state==NULL)
      alloc_state();

   cnfg_type=which_cnfg_type(in,1);

   fin=fopen(in,"rb");
   ir=fread(ldat,sizeof(int),16,fin);
   ir+=fread(state,sizeof(int),ns,fin);
   rlxs_reset(state);
   ir+=fread(state,sizeof(int),nd,fin);
   rlxd_reset(state);

   nplaq=(double)(6*L0*L1)*(double)(L2*L3);
   eps=sqrt(nplaq)*DBL_EPSILON;

   udb=NULL;
   if ((cnfg_type&1)!=0)
   {
      udb=udfld();
   }

   adb=NULL;
   if ((cnfg_type&2)!=0)
   {
      adb=adfld();
   }

   ie=1;
   if ((cnfg_type&1)!=0)
   {
      ir+=fread(&su3plaq0,sizeof(double),1,fin);
      ir+=fread(udb,sizeof(su3_dble),4*VOLUME,fin);
      set_flags(UPDATED_UD);
      
      if ((gauge()&1)!=0)
      {
         ie|=check_bc(64.0*DBL_EPSILON);
         error_root(ie!=1,1,"read_cnfg [archive.c]",
                    "Incompatible boundary conditions of the SU(3) field");
         ie=0;
         su3plaq0/=nplaq;
         su3plaq1=plaq_sum_dble(0)/nplaq;
         ie|=(fabs(su3plaq1-su3plaq0)>eps);
         set_bc();
         su3plaq1=plaq_sum_dble(0)/nplaq;
         ie|=(fabs(su3plaq1-su3plaq0)>eps);
         error_loc(ie!=0,1,"read_cnfg [archive.c]",
                   "Incorrect average plaquette of the SU(3) field");
      }
   }

   ie=1;
   if ((cnfg_type&2)!=0)
   {
      ir+=fread(&u1plaq0,sizeof(double),1,fin);
      ir+=fread(adb,sizeof(double),4*VOLUME,fin);
      set_flags(UPDATED_AD);
      
      if ((gauge()&2)!=0)
      {
         ie|=check_ad_bc(64.0*DBL_EPSILON);
         error_root(ie!=1,1,"read_cnfg [archive.c]",
                    "Incompatible boundary conditions of the U(1) field");
         ie=0;
         u1plaq0/=nplaq;
         u1plaq1=u1_plaq_sum_dble(0)/nplaq;
         ie|=(fabs(u1plaq1-u1plaq0)>eps);
         set_ad_bc();
         u1plaq1=u1_plaq_sum_dble(0)/nplaq;
         ie|=(fabs(u1plaq1-u1plaq0)>eps);
         error_loc(ie!=0,1,"read_cnfg [archive.c]",
                   "Incorrect average plaquette of the U(1) field");
      }
   }
   fclose(fin);

   if (cnfg_type!=3)
      error_loc(ir!=(17+ns+nd+4*VOLUME),1,"read_cnfg [archive.c]",
                "Incorrect read count");
   else
      error_loc(ir!=(18+ns+nd+8*VOLUME),1,"read_cnfg [archive.c]",
                "Incorrect read count");

   return cnfg_type;
}



static void get_su3_links(int iy)
{
   int y3,ifc;
   su3_dble *u,*v;

   v=ubuf;
   iy*=L3;

   if (ipt[iy]<(VOLUME/2))
   iy+=1;

   for (y3=0;y3<L3;y3+=2)
   {
      u=udb+8*(ipt[iy+y3]-(VOLUME/2));

      for (ifc=0;ifc<8;ifc++)
      {
         v[0]=u[0];
         v+=1;
         u+=1;
      }
   }
}

static void get_u1_links(int iy)
{
   int y3,ifc;
   double *au,*av;

   av=abuf;
   iy*=L3;

   if (ipt[iy]<(VOLUME/2))
   iy+=1;

   for (y3=0;y3<L3;y3+=2)
   {
      au=adb+8*(ipt[iy+y3]-(VOLUME/2));

      for (ifc=0;ifc<8;ifc++)
      {
         av[0]=au[0];
         av+=1;
         au+=1;
      }
   }
}


static void set_su3_links(int iy)
{
   int y3,ifc;
   su3_dble *u,*v;

   v=ubuf;
   iy*=L3;

   if (ipt[iy]<(VOLUME/2))
   iy+=1;

   for (y3=0;y3<L3;y3+=2)
   {
      u=udb+8*(ipt[iy+y3]-(VOLUME/2));

      for (ifc=0;ifc<8;ifc++)
      {
         u[0]=v[0];
         v+=1;
         u+=1;
      }
   }
}

static void set_u1_links(int iy)
{
   int y3,ifc;
   double *au,*av;

   av=abuf;
   iy*=L3;

   if (ipt[iy]<(VOLUME/2))
   iy+=1;

   for (y3=0;y3<L3;y3+=2)
   {
      au=adb+8*(ipt[iy+y3]-(VOLUME/2));

      for (ifc=0;ifc<8;ifc++)
      {
         au[0]=av[0];
         av+=1;
         au+=1;
      }
   }
}


static void cstar_su3_dble(int vol,su3_dble* u)
{
  su3_dble *um; 

  um=u+vol;
  for (;u<um;u++)
  {    
     (*u).c11.im *= -1.0;
     (*u).c12.im *= -1.0;
     (*u).c13.im *= -1.0;
     (*u).c21.im *= -1.0;
     (*u).c22.im *= -1.0;
     (*u).c23.im *= -1.0;
     (*u).c31.im *= -1.0;
     (*u).c32.im *= -1.0;
     (*u).c33.im *= -1.0;
  }
}

static void cstar_double(int vol,double* a)
{
  double *am; 

  am=a+vol;
  for (;a<am;a++)  (*a) *= -1.0;
}


void export_cnfg(char *out)
{
   int my_rank,np[4],n,iw;
   int iwa,dmy,tag0,tag1,cnfg_type=0;
   int x0,x1,x2,x3,y0,y1,y2,ix,iy;
   stdint_t lsize[4];
   double nplaq,su3plaq,u1plaq;
   MPI_Status stat;
   FILE *fout=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   check_machine();

   cnfg_type=gauge();
   nplaq=(double)(6*N0*N1)*(double)(N2*N3);

   su3plaq=0.0;
   u1plaq=0.0;
   if ((cnfg_type&1)!=0)
   {
      udb=udfld();

      if ((bc_cstar()!=0))
      {
         if (cpr[1]>=NPROC1/2)
            cstar_su3_dble(4*VOLUME,udb);
         set_flags(UPDATED_UD);
      }

      su3plaq=plaq_sum_dble(1)/nplaq;

      if ((bc_cstar()!=0))
      {
         if (cpr[1]>=NPROC1/2)
            cstar_su3_dble(4*VOLUME,udb);
         set_flags(UPDATED_UD);
      }

      if (ubuf==NULL)
         alloc_ubuf(my_rank);
   }
   if ((cnfg_type&2)!=0)
   {
      adb=adfld();

      if ((bc_cstar()!=0))
      {
         if (cpr[1]>=NPROC1/2)
            cstar_double(4*VOLUME,adb);
         set_flags(UPDATED_AD);
      }

      u1plaq=u1_plaq_sum_dble(1)/nplaq;

      if ((bc_cstar()!=0))
      {
         if (cpr[1]>=NPROC1/2)
            cstar_double(4*VOLUME,adb);
         set_flags(UPDATED_AD);
      }

      if (abuf==NULL)
         alloc_abuf(my_rank);
   }

   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();

   if (my_rank==0)
   {
      fout=fopen(out,"wb");
      error_root(fout==NULL,1,"export_cnfg [archive.c]",
                 "Unable to open output file");

      lsize[0]=(stdint_t)(N0);
      lsize[1]=(stdint_t)(N1);
      if (bc_cstar()!=0) lsize[1]/=2;
      lsize[2]=(stdint_t)(N2);
      lsize[3]=(stdint_t)(N3);


      if (endian==BIG_ENDIAN)
      {
         bswap_int(4,lsize);

         if ((cnfg_type&1)!=0)
         {
            bswap_double(1,&su3plaq);
         }
         if ((cnfg_type&2)!=0)
            bswap_double(1,&u1plaq);
      }

      iw=fwrite(lsize,sizeof(stdint_t),4,fout);

      error_root(iw!=4,1,"export_cnfg [archive.c]","Incorrect write count");
   }

   iwa=0;

   if ((cnfg_type&1)!=0)
   {
      if (my_rank==0)
      {
         iw=fwrite(&su3plaq,sizeof(double),1,fout);
         error_root(iw!=1,1,"export_cnfg [archive.c]","Incorrect write count");
      }
      for (ix=0;ix<(N0*N1*N2);ix++)
      {
         x0=ix/(N1*N2);
         x1=(ix/N2)%N1;
         x2=ix%N2;

         y0=x0%L0;
         y1=x1%L1;
         y2=x2%L2;
         iy=y2+L2*y1+L1*L2*y0;

         np[0]=x0/L0;
         np[1]=x1/L1;
         np[2]=x2/L2;

         if ((bc_cstar()!=0)&&(np[1]>=NPROC1/2))
            continue;

         for (x3=0;x3<N3;x3+=L3)
         {
            np[3]=x3/L3;
            n=ipr_global(np);
            if (my_rank==n)
               get_su3_links(iy);

            if (n>0)
            {
               if (my_rank==0)
               {
                  MPI_Send(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD);
                  MPI_Recv(ubuf,4*L3*18,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD,&stat);
               }
               else if (my_rank==n)
               {
                  MPI_Recv(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD,&stat);
                  MPI_Send(ubuf,4*L3*18,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD);
               }
            }

            if (my_rank==0)
            {
               if (endian==BIG_ENDIAN)
               bswap_double(4*L3*18,ubuf);
               iw=fwrite(ubuf,sizeof(su3_dble),4*L3,fout);
               iwa|=(iw!=(4*L3));
            }
         }
      }
   }

   if ((cnfg_type&2)!=0)
   {
      if (my_rank==0)
      {
         iw=fwrite(&u1plaq,sizeof(double),1,fout);
         error_root(iw!=1,1,"export_cnfg [archive.c]","Incorrect write count");
      }

      for (ix=0;ix<(N0*N1*N2);ix++)
      {
         x0=ix/(N1*N2);
         x1=(ix/N2)%N1;
         x2=ix%N2;

         y0=x0%L0;
         y1=x1%L1;
         y2=x2%L2;
         iy=y2+L2*y1+L1*L2*y0;

         np[0]=x0/L0;
         np[1]=x1/L1;
         np[2]=x2/L2;

         if ((bc_cstar()!=0)&&(np[1]>=NPROC1/2))
            continue;

         for (x3=0;x3<N3;x3+=L3)
         {
            np[3]=x3/L3;
            n=ipr_global(np);
            if (my_rank==n)
            get_u1_links(iy);

            if (n>0)
            {
               if (my_rank==0)
               {
                  MPI_Send(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD);
                  MPI_Recv(abuf,4*L3,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD,&stat);
               }
               else if (my_rank==n)
               {
                  MPI_Recv(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD,&stat);
                  MPI_Send(abuf,4*L3,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD);
               }
            }

            if (my_rank==0)
            {
               if (endian==BIG_ENDIAN)
               bswap_double(4*L3,abuf);
               iw=fwrite(abuf,sizeof(double),4*L3,fout);
               iwa|=(iw!=(4*L3));
            }
         }
      }
   }


   if (my_rank==0)
   {
      error_root(iwa!=0,1,"export_cnfg [archive.c]",
                 "Incorrect write count");
      fclose(fout);
   }
}


int import_cnfg(char *in)
{
   int my_rank,np[4],ir,ie,cnfg_type;
   int ira,dmy,tag0,tag1;
   int n0,n1,n2,n3,nc0,nc1,nc2,nc3;
   int x0,x1,x2,y0,y1,y2,y3,c0,c1,c2,ix,iy,ic;
   int n,k,l;
   stdint_t lsize[4];
   double nplaq,su3plaq0,su3plaq1,u1plaq0,u1plaq1,eps;
   MPI_Status stat;
   FILE *fin=NULL;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   cnfg_type=which_cnfg_type(in,0);

   check_machine();

   udb=NULL;
   if ((cnfg_type&1)!=0)
   {
      udb=udfld();
      if (ubuf==NULL)
         alloc_ubuf(my_rank);
   }

   adb=NULL;
   if ((cnfg_type&2)!=0)
   {
      adb=adfld();
      if (abuf==NULL)
         alloc_abuf(my_rank);
   }


   dmy=1;
   tag0=mpi_tag();
   tag1=mpi_tag();

   if (my_rank==0)
   {
      fin=fopen(in,"rb");
      error_root(fin==NULL,1,"import_cnfg [archive.c]",
      "Unable to open input file");

      ir=fread(lsize,sizeof(stdint_t),4,fin);
      error_root(ir!=4,1,"import_cnfg [archive.c]","Incorrect read count");

      if ((cnfg_type&1)!=0)
      {
         ir=fread(&su3plaq0,sizeof(double),1,fin);
         error_root(ir!=1,1,"import_cnfg [archive.c]","Incorrect read count");
      }
      else
      {
         su3plaq0=0.0;
         u1plaq0=0.0;
      }

      if (endian==BIG_ENDIAN)
      {
         bswap_int(4,lsize);
         bswap_double(1,&su3plaq0);
      }

      np[0]=(int)(lsize[0]);
      np[1]=(int)(lsize[1]);
      np[2]=(int)(lsize[2]);
      np[3]=(int)(lsize[3]);

      error_root((np[0]<1)||((N0%np[0])!=0)||
                 (np[1]<1)||((N1%np[1])!=0)||
                 (np[2]<1)||((N2%np[2])!=0)||
                 (np[3]<1)||((N3%np[3])!=0),1,"import_cnfg [archive.c]",
                 "Unexpected or incompatible lattice sizes");

      if (bc_cstar()!=0)
      error_root((((N1/2)%np[1])!=0),1,"import_cnfg [archive.c]",
                 "Incompatible lattice sizes for cstar bc's");

      error_root((np[0]!=N0)&&(bc_type()!=3),1,"import_cnfg [archive.c]",
                 "Periodic extension in time is only possible when\n"
                 "periodic boundary conditions are chosen");
   }
   else
   {
      np[0]=0;
      np[1]=0;
      np[2]=0;
      np[3]=0;
   }

   MPI_Bcast(np,4,MPI_INT,0,MPI_COMM_WORLD);

   if ((cnfg_type&1)!=0)
      MPI_Bcast(&su3plaq0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


   n0=np[0];
   n1=np[1];
   n2=np[2];
   n3=np[3];

   nc0=N0/n0;
   nc1=N1/n1;
   nc2=N2/n2;
   nc3=N3/n3;
   ira=0;

   if ((cnfg_type&1)!=0)
   {
      for (ix=0;ix<(n0*n1*n2);ix++)
      {
         x0=ix/(n1*n2);
         x1=(ix/n2)%n1;
         x2=ix%n2;

         if (my_rank==0)
         {
            n=4*n3;
            ir=fread(vbuf,sizeof(su3_dble),n,fin);
            ira|=(ir!=n);

            if (endian==BIG_ENDIAN)
               bswap_double(n*18,vbuf);

            for (k=1;k<nc3;k++)
            {
               for (l=0;l<n;l++)
               vbuf[k*n+l]=vbuf[l];
            }
         }

         for (ic=0;ic<(nc0*nc1*nc2);ic++)
         {
            c0=ic/(nc1*nc2);
            c1=(ic/nc2)%nc1;
            c2=ic%nc2;

            y0=x0+c0*n0;
            y1=x1+c1*n1;
            y2=x2+c2*n2;
            iy=(y2%L2)+L2*(y1%L1)+L1*L2*(y0%L0);

            np[0]=y0/L0;
            np[1]=y1/L1;
            np[2]=y2/L2;

            for (y3=0;y3<N3;y3+=L3)
            {
               np[3]=y3/L3;
               n=ipr_global(np);

               if (n>0)
               {
                  if (my_rank==0)
                  {
                     MPI_Send(vbuf+4*y3,4*L3*18,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD);
                     MPI_Recv(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD,&stat);
                  }
                  else if (my_rank==n)
                  {
                     MPI_Recv(ubuf,4*L3*18,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD,&stat);
                     MPI_Send(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD);
                  }
               }
               else if (my_rank==0)
                  for (l=0;l<(4*L3);l++)
                     ubuf[l]=vbuf[4*y3+l];

               if (my_rank==n)
               {
                  set_su3_links(iy);
               }
            }
         }
      }
      set_flags(UPDATED_UD);
   }

   if ((cnfg_type&2)!=0)
   {
      if (my_rank==0)
      {
         ir=fread(&u1plaq0,sizeof(double),1,fin);
         error_root(ir!=1,1,"import_cnfg [archive.c]","Incorrect read count");
         if (endian==BIG_ENDIAN)
            bswap_double(1,&u1plaq0);
      }

      MPI_Bcast(&u1plaq0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      ira=0;
      for (ix=0;ix<(n0*n1*n2);ix++)
      {
         x0=ix/(n1*n2);
         x1=(ix/n2)%n1;
         x2=ix%n2;

         if (my_rank==0)
         {
            n=4*n3;
            ir=fread(vabuf,sizeof(double),n,fin);
            ira|=(ir!=n);

            if (endian==BIG_ENDIAN)
               bswap_double(n,vabuf);

            for (k=1;k<nc3;k++)
            {
               for (l=0;l<n;l++)
               vabuf[k*n+l]=vabuf[l];
            }
         }

         for (ic=0;ic<(nc0*nc1*nc2);ic++)
         {
            c0=ic/(nc1*nc2);
            c1=(ic/nc2)%nc1;
            c2=ic%nc2;

            y0=x0+c0*n0;
            y1=x1+c1*n1;
            y2=x2+c2*n2;
            iy=(y2%L2)+L2*(y1%L1)+L1*L2*(y0%L0);

            np[0]=y0/L0;
            np[1]=y1/L1;
            np[2]=y2/L2;

            for (y3=0;y3<N3;y3+=L3)
            {
               np[3]=y3/L3;
               n=ipr_global(np);

               if (n>0)
               {
                  if (my_rank==0)
                  {
                     MPI_Send(vabuf+4*y3,4*L3,MPI_DOUBLE,n,tag1,MPI_COMM_WORLD);
                     MPI_Recv(&dmy,1,MPI_INT,n,tag0,MPI_COMM_WORLD,&stat);
                  }
                  else if (my_rank==n)
                  {
                     MPI_Recv(abuf,4*L3,MPI_DOUBLE,0,tag1,MPI_COMM_WORLD,&stat);
                     MPI_Send(&dmy,1,MPI_INT,0,tag0,MPI_COMM_WORLD);
                  }
               }
               else if (my_rank==0)
                  for (l=0;l<(4*L3);l++)
                     abuf[l]=vabuf[4*y3+l];

               if (my_rank==n)
               {
                  set_u1_links(iy);
               }
            }
         }
      }
      set_flags(UPDATED_AD);
   }


   if (my_rank==0)
   {
      error_root(ira!=0,1,"import_cnfg [archive.c]","Incorrect read count");
      fclose(fin);
   }

   nplaq=(double)(6*N0*N1)*(double)(N2*N3);
   eps=sqrt(nplaq)*DBL_EPSILON;

   if (((cnfg_type&1)!=0)&&((gauge()&1)!=0))
   {
      ie=check_bc(64.0*DBL_EPSILON);
      error_root(ie!=1,1,"import_cnfg [archive.c]",
                 "Incompatible SU(3) boundary conditions");
      ie=0;
      su3plaq1=plaq_sum_dble(1)/nplaq;
      ie|=(fabs(su3plaq1-su3plaq0)>eps);
      set_bc();
      su3plaq1=plaq_sum_dble(1)/nplaq;
      ie|=(fabs(su3plaq1-su3plaq0)>eps);
      error_root(ie!=0,1,"import_cnfg [archive.c]",
                 "Incorrect average plaquette for the SU(3) gauge field");
   }

   if (((cnfg_type&2)!=0)&&((gauge()&2)!=0))
   {
      ie=check_ad_bc(64.0*DBL_EPSILON);
      error_root(ie!=1,1,"import_cnfg [archive.c]",
                 "Incompatible U(1) boundary conditions");
      ie=0;
      u1plaq1=u1_plaq_sum_dble(1)/nplaq;
      ie|=(fabs(u1plaq1-u1plaq0)>eps);
      set_ad_bc();
      u1plaq1=u1_plaq_sum_dble(1)/nplaq;
      ie|=(fabs(u1plaq1-u1plaq0)>eps);
      error_root(ie!=0,1,"import_cnfg [archive.c]",
                "Incorrect average plaquette for the U(1) gauge field");
   }
   
   if ((bc_cstar()!=0))
   {
      if ((cnfg_type&1)!=0)
      {
         if (cpr[1]>=NPROC1/2)
            cstar_su3_dble(4*VOLUME,udb);
         set_flags(UPDATED_UD);
      }
      if ((cnfg_type&2)!=0)
      {
         if (cpr[1]>=NPROC1/2)
            cstar_double(4*VOLUME,adb);
         set_flags(UPDATED_AD);
      }
   }

   return cnfg_type;
}
