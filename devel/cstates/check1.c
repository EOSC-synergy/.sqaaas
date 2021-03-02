
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2017 Nazario Tantalo
*               2020 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the cfactor for a particular choice of the U(1) field
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "random.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "u1flds.h"
#include "cstates.h"
#include "global.h"


#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1/2)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int np[4],bo[4];
static double ac[16];

static double afld(int *x,int mu)
{
	int nu;
	double xt[4],phi;

	xt[0]=(double)(safe_mod(x[0],N0));
	xt[1]=(double)(safe_mod(x[1],N1));
	xt[2]=(double)(safe_mod(x[2],N2));
	xt[3]=(double)(safe_mod(x[3],N3));

	phi=0.0;
	for (nu=0;nu<4;nu++)
		phi+=ac[nu+4*mu]*xt[nu];

	if( bc_cstar()>=1 && (x[1]<0 || x[1]>=N1) )
		phi*=-1.0;

	if( bc_cstar()>=2 && (x[2]<0 || x[2]>=N2) )
		phi*=-1.0;

	if( bc_cstar()>=3 && (x[3]<0 || x[3]>=N3) )
		phi*=-1.0;

	return phi;
}

static void choose_ac(void)
{
	int mu;
	double r[16];

	ranlxd(r,16);
	MPI_Bcast(r,16,MPI_DOUBLE,0,MPI_COMM_WORLD);

	for(mu=0;mu<16;++mu)
		ac[mu]=(double)((int)(3.0*r[mu])-1);
}


static void set_ad(void)
{
	int x[4];
	int x0,x1,x2,x3;
	int ix,ifc;
	double phi;
	double *adb,*a;

	adb=adfld();

	for (x0=0;x0<L0;x0++)
	for (x1=0;x1<L1;x1++)
	for (x2=0;x2<L2;x2++)
	for (x3=0;x3<L3;x3++)
	{
		ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

		if (ix>=(VOLUME/2))
		{
			x[0]=bo[0]+x0;
			x[1]=bo[1]+x1;
			x[2]=bo[2]+x2;
			x[3]=bo[3]+x3;

			a=adb+8*(ix-(VOLUME/2));

			for (ifc=0;ifc<8;ifc++)
			{
				if (ifc&0x1)
					x[ifc/2]-=1;

				phi=afld(x,ifc/2);

				if (ifc&0x1)
					x[ifc/2]+=1;

				(*a)=phi;
				a+=1;
			}		
		}
	}
	set_ad_bc();
	orbi_cpy_ad();

	set_flags(UPDATED_AD);
}


static void check_cstring_ad(int mu)
{
	int ix,x[4],nu;
	int x0,x1,x2,x3;
	double diff,sum,term;
	double gdiff,gsum;
	double *cstring;

	cstring=get_cfactor_phase(0,mu);

	diff=0.0;
	sum=0.0;
	for (x0=0;x0<L0;x0++)
	for (x2=0;x2<L2;x2++)
	for (x3=0;x3<L3;x3++)
	for (x1=0;x1<L1;x1++)
	{
		ix=ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];

		x[0]=bo[0]+x0;
		x[1]=bo[1]+x1;
		x[2]=bo[2]+x2;
		x[3]=bo[3]+x3;

		term=0.0;
		for(nu=0;nu<4;++nu)
			term+=((double)(2*x[mu]-np[mu]))*ac[nu+4*mu]*x[nu];
		term-=((double)(2*x[mu]-np[mu]))*ac[mu+4*mu]*x[mu];
		term+=((double)(x[mu]*(x[mu]-1)))*ac[mu+4*mu];
		term-=((double)(np[mu]*(np[mu]-1)))*ac[mu+4*mu]*0.5;

		if(cpr[1]>=NPROC1/2)
			term*=-1.0;

		if (bc_type()==1 && global_time(ix)==0 )
			term=0.0;

		diff+=fabs(term-cstring[ix]);
		sum+=fabs(cstring[ix]);
	}

	MPI_Reduce(&sum,&gsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&gsum,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Reduce(&diff,&gdiff,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&gdiff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	message("mu=%2d, sum= %.1e, dev= %.1e, rel_dev= %.1e\n",
				mu,gsum,gdiff,gdiff/gsum);
}


int main(int argc,char *argv[])
{
	int my_rank,bc,cs,mpr1,mu;
	double phi[2],phi_prime[2];
	FILE *flog=NULL;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	if (my_rank==0)
	{
		flog=freopen("check1.log","w",stdout);

		printf("\n");
		printf("Check of the cfactor for a particular choice of the U(1) field\n");
		printf("--------------------------------------------------------------\n\n");

		printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
		printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
		printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);

		bc=find_opt(argc,argv,"-bc");

		if (bc!=0)
			error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
							"Syntax: check1 -cs <cstar> [-bc <type>]");

			cs=find_opt(argc,argv,"-cs");
			error_root(cs==0 || sscanf(argv[cs+1],"%d",&cs)!=1,1,"main [check1.c]",
							"Syntax: check1 -cs <cstar> [-bc <type>]");
	}

	set_flds_parms(2,0);
	print_flds_parms();

	MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&cs,1,MPI_INT,0,MPI_COMM_WORLD);
	phi[0]=0.123;
	phi[1]=-0.534;
	phi_prime[0]=0.912;
	phi_prime[1]=0.078;
	set_bc_parms(bc,cs,phi,phi_prime,0.573,-1.827);
	print_bc_parms();

	start_ranlux(0,123456);
	geometry();

	np[0]=N0;
	np[1]=N1;
	np[2]=N2;
	np[3]=N3;

	mpr1=cpr[1];
	if(mpr1>=NPROC1/2)
		mpr1=(cpr[1]+NPROC1/2)%NPROC1;

	bo[0]=cpr[0]*L0;
	bo[1]=mpr1*L1;
	bo[2]=cpr[2]*L2;
	bo[3]=cpr[3]*L3;

	choose_ac();
	set_ad();

	for(mu=1;mu<4;++mu)
		if(bc_cstar()>=mu)
			check_cstring_ad(mu);

	message("\n");

	if (my_rank==0)
		fclose(flog);

	MPI_Finalize();
	exit(0);
}
