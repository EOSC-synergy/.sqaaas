
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2005, 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Consistency checks on the global index arrays cpr,...,map
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "global.h"

#define NPROC_BLK (NPROC0_BLK*NPROC1_BLK*NPROC2_BLK*NPROC3_BLK)

static int ip_test[NPROC];
static int ix_test[VOLUME];
static int ia[2][9];

static void set_ia(void)
{
   int ifc;

   ia[0][0]=0;
   ia[0][1]=ia[0][0]+(FACE0/2);
   ia[0][2]=ia[0][1]+(FACE0/2);
   ia[0][3]=ia[0][2]+(FACE1/2);
   ia[0][4]=ia[0][3]+(FACE1/2);
   ia[0][5]=ia[0][4]+(FACE2/2);
   ia[0][6]=ia[0][5]+(FACE2/2);
   ia[0][7]=ia[0][6]+(FACE3/2);
   ia[0][8]=ia[0][7]+(FACE3/2);

   for (ifc=0;ifc<9;ifc++)
      ia[1][ifc]=ia[0][ifc]+(BNDRY/2);
}


int main(int argc,char *argv[])
{
   double phi[2];
   int my_rank,itest,cs;
   int in,ir,is,n[4],NMIRR;
   int mu,ix,x0,x1,x2,x3;
   int iy0,iy1,iy2,iy3,iz0,iz1,iz2,iz3;
   int g[4],gx[VOLUME+BNDRY][4],iy;
   int rbuf[5*BNDRY+5],sbuf[5*BNDRY+5],rlen,slen;
   int pc[NPROC][4];
   int tag;
   MPI_Status stat;
   FILE *flog=NULL;
   const int ns=48;
   const int sh[48][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},
                        {-1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1},
                        {2,0,0,0},{0,2,0,0},{0,0,2,0},{0,0,0,2},
                        {-2,0,0,0},{0,-2,0,0},{0,0,-2,0},{0,0,0,-2},
                        {3,0,0,0},{0,3,0,0},{0,0,3,0},{0,0,0,3},
                        {-3,0,0,0},{0,-3,0,0},{0,0,-3,0},{0,0,0,-3},
                        {-6,-9,1,-9},{-2,0,7,-3},{-5,-4,-1,6},{9,0,6,-5},
                        {-2,0,5,5},{-9,-2,4,-1},{-4,-2,7,-2},{4,-9,-3,-4},
                        {1,-9,-1,0},{-4,-7,7,-4},{-5,8,8,-3},{7,3,3,-3},
                        {3,0,2,-3},{-2,-7,-5,-3},{5,6,0,6},{-5,0,-2,-9},
                        {1,-7,-7,1},{1,-3,-4,-6},{-5,-5,-1,-1},{-4,-3,0,0},
                        {-2,0,9,-2},{0,-6,1,-1},{8,0,7,-7},{1,-5,2,0}};

   
   
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Consistency checks on the global index arrays cpr,...,map\n");
      printf("---------------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d grid blocks\n\n",
             NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);

      cs=find_opt(argc,argv,"-cs");

      if (cs!=0)
         error_root(sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-cs <cstar>]");
   }

   MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
   phi[0]=0.0;
   phi[1]=0.0;
   set_bc_parms(3,cs,phi,phi,0.0,0.0);
   print_bc_parms();

   geometry();
   set_ia();

   error(my_rank!=ipr_global(cpr),1,
         "main [check1.c]","Processor coordinates are incorrect");
   if (my_rank==0)
      printf("Process rank and coordinates are consistent with ipr_global() ... OK\n");
   

   if (my_rank==0)
   {
      for (in=0;in<NPROC;in++)
      {
         ir=in;
         n[0]=ir%NPROC0;
         ir/=NPROC0;
         n[1]=ir%NPROC1;
         ir/=NPROC1;
         n[2]=ir%NPROC2;
         ir/=NPROC2;
         n[3]=ir%NPROC3;

         ip_test[in]=ipr_global(n);
      }
   }

   MPI_Bcast(ip_test,NPROC,MPI_INT,0,MPI_COMM_WORLD);

   itest=0;

   for (in=0;in<NPROC;in++)
   {
      ir=in;
      n[0]=ir%NPROC0;
      ir/=NPROC0;
      n[1]=ir%NPROC1;
      ir/=NPROC1;
      n[2]=ir%NPROC2;
      ir/=NPROC2;
      n[3]=ir%NPROC3;

      if (ip_test[in]!=ipr_global(n))
         itest=1;

      n[0]-=(n[0]%NPROC0_BLK);
      n[1]-=(n[1]%NPROC1_BLK);
      n[2]-=(n[2]%NPROC2_BLK);
      n[3]-=(n[3]%NPROC3_BLK);
      ir=ipr_global(n);

      NMIRR=1;
      if ((bc_cstar()!=0)&&((NPROC1/NPROC1_BLK)%2==0))
      {
         NMIRR=2;
         if (n[1]>=NPROC1/2)
            n[1]-=NPROC1/2;
      }


      if ((ip_test[in]<ir)||(ip_test[in]>=(ir+NPROC_BLK*NMIRR)))
         itest=2;
   }

   error(itest==1,1,
         "main [check1.c]","ipr_global is process dependent");
   if (my_rank==0)
      printf("ipr_global() returns the same values on all processors ... OK\n");

   error(itest==2,1,
         "main [check1.c]","Processes are not properly blocked");
   if (my_rank==0)
      printf("Processes are properly blocked ... OK\n");

   itest=0;
   for (n[0]=0;n[0]<NPROC0;n[0]++)
   {
      for (n[1]=0;n[1]<NPROC1;n[1]++)
      {
         for (n[2]=0;n[2]<NPROC2;n[2]++)
         {
            for (n[3]=0;n[3]<NPROC3;n[3]++)
            {
               for (is=0;is<ns;is++)
               {
                  in=ipr_global(n);
                  if (((bc_cstar()==2)&&(safe_mod(sh[is][2],2)==1))||
                      ((bc_cstar()==3)&&(safe_mod(sh[is][2]+sh[is][3],2)==1)))
                  {
                     n[1]=(n[1]+NPROC1/2)%NPROC1;
                     in=ipr_global(n);
                     n[1]=(n[1]+NPROC1/2)%NPROC1;
                  }

                  n[0]+=sh[is][0]*NPROC0;
                  n[1]+=sh[is][1]*NPROC1;
                  n[2]+=sh[is][2]*NPROC2;
                  n[3]+=sh[is][3]*NPROC3;
                  itest|=(in!=ipr_global(n));
                  n[0]-=sh[is][0]*NPROC0;
                  n[1]-=sh[is][1]*NPROC1;
                  n[2]-=sh[is][2]*NPROC2;
                  n[3]-=sh[is][3]*NPROC3;
               }
            }
         }
      }
   }
   
   error(itest!=0,1,
         "main [check1.c]","Global topology is incorrect");
   if (my_rank==0)
      printf("Global topology is correctly built ... OK\n");


   n[0]=cpr[0];
   n[1]=cpr[1];
   n[2]=cpr[2];
   n[3]=cpr[3];

   for (mu=0;mu<4;mu++)
   {
      n[mu]-=1;
      if (npr[2*mu]!=ipr_global(n))
         itest=1;
      n[mu]+=2;
      if (npr[2*mu+1]!=ipr_global(n))
         itest=1;
      n[mu]-=1;
   }

   error(itest==1,1,
         "main [check1.c]","npr is inconsistent with ipr_global");
   if (my_rank==0)
      printf("npr is consistent with ipr_global() ... OK\n");

   for (ix=0;ix<VOLUME;ix++)
      ix_test[ix]=0;

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               if ((ix<0)||(ix>=VOLUME))
                  itest=1;
               else
                  ix_test[ix]+=1;
            }
         }
      }
   }

   error(itest==1,1,
         "main [check1.c]","The index ipt is out of range");
   if (my_rank==0)
      printf("The index ipt takes integer values in [0,VOLUME] on the local lattice... OK\n");

   for (ix=0;ix<VOLUME;ix++)
   {
      if (ix_test[ix]!=1)
         itest=1;
   }

   error(itest==1,1,
         "main [check1.c]","The index ipt is not one-to-one");
   if (my_rank==0)
      printf("The map ipt is one-to-one on the local lattice ... OK\n");

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
               ir=(x0+x1+x2+x3)%2;

               if (((ir==0)&&(ix>=(VOLUME/2)))||((ir==1)&&(ix<(VOLUME/2))))
                  itest=1;

               ir=(ir+1)%2;
               iy0=iup[ix][0];
               iz0=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*((x0+1)%L0)];

               if ((x0==(L0-1))&&(NPROC0>1))
               {
                  iy0-=VOLUME;
                  if ((iy0<ia[ir][1])||(iy0>=ia[ir][2]))
                     itest=2;
                  else
                     iy0=map[iy0];
               }

               iy1=iup[ix][1];
               iz1=ipt[x3+L3*x2+L2*L3*((x1+1)%L1)+L1*L2*L3*x0];

               if ((x1==(L1-1))&&(NPROC1>1))
               {
                  iy1-=VOLUME;
                  if ((iy1<ia[ir][3])||(iy1>=ia[ir][4]))
                     itest=2;
                  else
                     iy1=map[iy1];
               }

               iy2=iup[ix][2];
               iz2=ipt[x3+L3*((x2+1)%L2)+L2*L3*x1+L1*L2*L3*x0];

               if ((x2==(L2-1))&&(NPROC2>1))
               {
                  iy2-=VOLUME;
                  if ((iy2<ia[ir][5])||(iy2>=ia[ir][6]))
                     itest=2;
                  else
                     iy2=map[iy2];
               }

               iy3=iup[ix][3];
               iz3=ipt[((x3+1)%L3)+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               if ((x3==(L3-1))&&(NPROC3>1))
               {
                  iy3-=VOLUME;
                  if ((iy3<ia[ir][7])||(iy3>=ia[ir][8]))
                     itest=2;
                  else
                     iy3=map[iy3];
               }

               if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3))
                  itest=3;

               iy0=idn[ix][0];
               iz0=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*((x0+L0-1)%L0)];

               if ((x0==0)&&(NPROC0>1))
               {
                  iy0-=VOLUME;
                  if ((iy0<ia[ir][0])||(iy0>=ia[ir][1]))
                     itest=4;
                  else
                     iy0=map[iy0];
               }

               iy1=idn[ix][1];
               iz1=ipt[x3+L3*x2+L2*L3*((x1+L1-1)%L1)+L1*L2*L3*x0];

               if ((x1==0)&&(NPROC1>1))
               {
                  iy1-=VOLUME;
                  if ((iy1<ia[ir][2])||(iy1>=ia[ir][3]))
                     itest=4;
                  else
                     iy1=map[iy1];
               }

               iy2=idn[ix][2];
               iz2=ipt[x3+L3*((x2+L2-1)%L2)+L2*L3*x1+L1*L2*L3*x0];

               if ((x2==0)&&(NPROC2>1))
               {
                  iy2-=VOLUME;
                  if ((iy2<ia[ir][4])||(iy2>=ia[ir][5]))
                     itest=4;
                  else
                     iy2=map[iy2];
               }

               iy3=idn[ix][3];
               iz3=ipt[((x3+L3-1)%L3)+L3*x2+L2*L3*x1+L1*L2*L3*x0];

               if ((x3==0)&&(NPROC3>1))
               {
                  iy3-=VOLUME;
                  if ((iy3<ia[ir][6])||(iy3>=ia[ir][7]))
                     itest=4;
                  else
                     iy3=map[iy3];
               }

               if ((iy0!=iz0)||(iy1!=iz1)||(iy2!=iz2)||(iy3!=iz3))
                  itest=5;
            }
         }
      }
   }

   error(itest==1,1,
         "main [check1.c]","The index ipt does not respect eo ordering");
   if (my_rank==0)
      printf("The index ipt respects eo ordering ... OK\n");
   error(itest==2,1,
         "main [check1.c]","The index iup is out of range at the boundaries");
   if (my_rank==0)
      printf("The index iup takes values in the correct range at the boundaries ... OK\n");
   error(itest==3,1,
         "main [check1.c]","The index iup (combined with map) is incorrect");
   if (my_rank==0)
      printf("The index iup is consistent with the array map ... OK\n");
   error(itest==4,1,
         "main [check1.c]","The index idn is out of range at the boundaries");
   if (my_rank==0)
      printf("The index idn takes values in the correct range at the boundaries ... OK\n");
   error(itest==5,1,
         "main [check1.c]","The index idn (combined with map) is incorrect");
   if (my_rank==0)
      printf("The index idn is consistent with the array map ... OK\n");

   for (ix=0;ix<VOLUME+BNDRY;ix++)
      gx[ix][0]=-1;

   itest=0;
   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
               
               g[0]=cpr[0]*L0+x0;
               g[1]=cpr[1]*L1+x1;
               g[2]=cpr[2]*L2+x2;
               g[3]=cpr[3]*L3+x3;
               
               if(gx[ix][0]==-1)
               {
                  gx[ix][0]=g[0];
                  gx[ix][1]=g[1];
                  gx[ix][2]=g[2];
                  gx[ix][3]=g[3];
               }
               else
               {
                  if((gx[ix][0]!=g[0])||(gx[ix][1]!=g[1])||
                     (gx[ix][2]!=g[2])||(gx[ix][3]!=g[3]))
                  {
                     itest=1;
                  }
               }
            }
         }
      }
   }
   
   error(itest==1,1,
         "main [check1.c]","Global coordinates are not uniquely assigned (1)");
   
   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               for (mu=0;mu<4;mu++)
               {
                  ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
                  iy=iup[ix][mu];
                  
                  g[0]=gx[ix][0];
                  g[1]=gx[ix][1];
                  g[2]=gx[ix][2];
                  g[3]=gx[ix][3];
                  g[mu]+=1;
                  
                  g[0]=safe_mod(g[0],L0*NPROC0);
                  if (bc_cstar()>=2)
                  {
                     g[2]=safe_mod(g[2],2*L2*NPROC2);
                     if (g[2]>=L2*NPROC2)
                     {
                        g[2]-=L2*NPROC2;
                        g[1]+=L1*NPROC1/2;
                     }
                  }
                  else
                     g[2]=safe_mod(g[2],L2*NPROC2);
                  if (bc_cstar()>=3)
                  {
                     g[3]=safe_mod(g[3],2*L3*NPROC3);
                     if (g[3]>=L3*NPROC3)
                     {
                        g[3]-=L3*NPROC3;
                        g[1]+=L1*NPROC1/2;
                     }
                  }
                  else
                     g[3]=safe_mod(g[3],L3*NPROC3);
                  g[1]=safe_mod(g[1],L1*NPROC1);
                  
                  if(gx[iy][0]==-1)
                  {
                     gx[iy][0]=g[0];
                     gx[iy][1]=g[1];
                     gx[iy][2]=g[2];
                     gx[iy][3]=g[3];
                  }
                  else
                  {
                     if((gx[iy][0]!=g[0])||(gx[iy][1]!=g[1])||
                        (gx[iy][2]!=g[2])||(gx[iy][3]!=g[3]))
                     {
                        itest=1;
                     }
                  }
               }
            }
         }
      }
   }
   
   error(itest==1,1,
         "main [check1.c]","Global coordinates are not uniquely assigned (2)");

   for (x0=0;x0<L0;x0++)
   {
      for (x1=0;x1<L1;x1++)
      {
         for (x2=0;x2<L2;x2++)
         {
            for (x3=0;x3<L3;x3++)
            {
               for (mu=0;mu<4;mu++)
               {
                  ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];
                  iy=idn[ix][mu];
                  
                  g[0]=gx[ix][0];
                  g[1]=gx[ix][1];
                  g[2]=gx[ix][2];
                  g[3]=gx[ix][3];
                  g[mu]-=1;
                  
                  g[0]=safe_mod(g[0],L0*NPROC0);
                  if (bc_cstar()>=2)
                  {
                     g[2]=safe_mod(g[2],2*L2*NPROC2);
                     if (g[2]>=L2*NPROC2)
                     {
                        g[2]-=L2*NPROC2;
                        g[1]+=L1*NPROC1/2;
                     }
                  }
                  else
                     g[2]=safe_mod(g[2],L2*NPROC2);
                  if (bc_cstar()>=3)
                  {
                     g[3]=safe_mod(g[3],2*L3*NPROC3);
                     if (g[3]>=L3*NPROC3)
                     {
                        g[3]-=L3*NPROC3;
                        g[1]+=L1*NPROC1/2;
                     }
                  }
                  else
                     g[3]=safe_mod(g[3],L3*NPROC3);
                  g[1]=safe_mod(g[1],L1*NPROC1);
                  
                  if(gx[iy][0]==-1)
                  {
                     gx[iy][0]=g[0];
                     gx[iy][1]=g[1];
                     gx[iy][2]=g[2];
                     gx[iy][3]=g[3];
                  }
                  else
                  {
                     if((gx[iy][0]!=g[0])||(gx[iy][1]!=g[1])||
                        (gx[iy][2]!=g[2])||(gx[iy][3]!=g[3]))
                     {
                        itest=1;
                     }
                  }
               }
            }
         }
      }
   }
   
   error(itest==1,1,
         "main [check1.c]","Global coordinates are not uniquely assigned (3)");
   if (my_rank==0)
      printf("Global coordinates are uniquely assigned ... OK\n");
   
   for (ix=0;ix<VOLUME+BNDRY;ix++)
   {
      if(gx[ix][0]==-1) itest=2;
   }
   
   error(itest==2,1,
         "main [check1.c]","Some points have no global coordinates");
   if (my_rank==0)
      printf("Global coordinates are assigned to all points ... OK\n");
   
   if(NPROC>1)
   {
      for (mu=0;mu<4;mu++)
      {
         slen=0;
         for (ix=0;ix<VOLUME;ix++)
         {
            iy=iup[ix][mu];
            if(iy>=VOLUME)
            {
               sbuf[5*slen+0]=gx[iy][0];
               sbuf[5*slen+1]=gx[iy][1];
               sbuf[5*slen+2]=gx[iy][2];
               sbuf[5*slen+3]=gx[iy][3];
               sbuf[5*slen+4]=map[iy-VOLUME];
               slen++;
            }
         }
         
         tag=mpi_tag();
         MPI_Sendrecv(&slen,1,MPI_INT,npr[2*mu+1],tag,
                      &rlen,1,MPI_INT,npr[2*mu],tag,
                      MPI_COMM_WORLD,&stat);
         tag=mpi_tag();
         MPI_Sendrecv(sbuf,5*slen,MPI_INT,npr[2*mu+1],tag,
                      rbuf,5*rlen,MPI_INT,npr[2*mu],tag,
                      MPI_COMM_WORLD,&stat);
         
         for (ix=0;ix<rlen;ix++)
         {
            iy=rbuf[5*ix+4];
            if((gx[iy][0]!=rbuf[5*ix+0])||(gx[iy][1]!=rbuf[5*ix+1])||
               (gx[iy][2]!=rbuf[5*ix+2])||(gx[iy][3]!=rbuf[5*ix+3]))
            {
               itest=3;
            }
         }
      }

      for (mu=0;mu<4;mu++)
      {
         slen=0;
         for (ix=0;ix<VOLUME;ix++)
         {
            iy=idn[ix][mu];
            if(iy>=VOLUME)
            {
               sbuf[5*slen+0]=gx[iy][0];
               sbuf[5*slen+1]=gx[iy][1];
               sbuf[5*slen+2]=gx[iy][2];
               sbuf[5*slen+3]=gx[iy][3];
               sbuf[5*slen+4]=map[iy-VOLUME];
               slen++;
            }
         }
         
         tag=mpi_tag();
         MPI_Sendrecv(&slen,1,MPI_INT,npr[2*mu],tag,
                      &rlen,1,MPI_INT,npr[2*mu+1],tag,
                      MPI_COMM_WORLD,&stat);
         tag=mpi_tag();
         MPI_Sendrecv(sbuf,5*slen,MPI_INT,npr[2*mu],tag,
                      rbuf,5*rlen,MPI_INT,npr[2*mu+1],tag,
                      MPI_COMM_WORLD,&stat);
         
         for (ix=0;ix<rlen;ix++)
         {
            iy=rbuf[5*ix+4];
            if((gx[iy][0]!=rbuf[5*ix+0])||(gx[iy][1]!=rbuf[5*ix+1])||
               (gx[iy][2]!=rbuf[5*ix+2])||(gx[iy][3]!=rbuf[5*ix+3]))
            {
               itest=3;
            }
         }
      }

      error(itest==3,1,
            "main [check1.c]","Global coordinates are not consistent with iup/idn/map");
      if (my_rank==0)
         printf("Global coordinates are consistent with iup/idn/map ... OK\n");
   }

   
   if (bc_cstar()!=0)
   {
      sbuf[0]=gx[0][0];
      sbuf[1]=gx[0][1];
      sbuf[2]=gx[0][2];
      sbuf[3]=gx[0][3];
      
      tag=mpi_tag();
      MPI_Sendrecv(sbuf,4,MPI_INT,get_mirror_rank(),tag,
                   rbuf,4,MPI_INT,get_mirror_rank(),tag,
                   MPI_COMM_WORLD,&stat);
      
      if ((sbuf[0]!=rbuf[0])||(sbuf[2]!=rbuf[2])||(sbuf[3]!=rbuf[3]))
         itest=4;
      if ((cpr[1]<NPROC1/2)&&(sbuf[1]!=rbuf[1]-NPROC1*L1/2))
         itest=4;
      if ((cpr[1]>=NPROC1/2)&&(sbuf[1]!=rbuf[1]+NPROC1*L1/2))
         itest=4;

      error(itest==4,1,
            "main [check1.c]","Global coordinates are not consistent with get_mirror_rank");
      if (my_rank==0)
         printf("Global coordinates are consistent with get_mirror_rank ... OK\n");
   }
   
   if(my_rank==0)
   {
      pc[0][0]=cpr[0];
      pc[0][1]=cpr[1];
      pc[0][2]=cpr[2];
      pc[0][3]=cpr[3];
   }
   for(in=1;in<NPROC;in++)
   {
      tag=mpi_tag();
      if(in==my_rank)
         MPI_Send(cpr,4,MPI_INT,0,tag,MPI_COMM_WORLD);
      if(my_rank==0)
         MPI_Recv(pc[in],4,MPI_INT,in,tag,MPI_COMM_WORLD,&stat);
   }
   
   if(my_rank==0)
   {
      printf("\n    Rank    Process coords\n");
      printf("------------------------------------\n");
      for(in=0;in<NPROC;in++)
      {
         printf("%8d    %d %d %d %d\n",ipr_global(pc[in]),pc[in][0],pc[in][1],pc[in][2],pc[in][3]);
      }
   }

   if (my_rank==0)
   {
      printf("\nThe lattice is correctly mapped by the global arrays\n\n");
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
