
/*******************************************************************************
*
* File minmax.c
*
* Copyright (C) 2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* MinMax rational approximation of ( x + mu^2 )^{p/q}
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpf2mpfr.h"
#include "utils.h"
#include "remes.h"

remes_t rem;
char logname[NAME_SIZE];
FILE *flog=NULL;


static void parse_input(int argc,char *argv[])
{
   int iopt,slen,n;
   int p,q,s,f,m,verbose;
   int precision,relflag;
   double ra,rb,mu,goal,alpha,beta;
   char syntax[NAME_SIZE];
   char syscmd[NAME_SIZE];   

   sprintf(syntax,"sintax: minmax -p <int> -q <int> -ra <float> -rb <float> -goal <float> [-mu <float>] [-abs] [-prec <int>] [-outpath <string>] [-nopen] [-nstart <int> <float> <float>]");
      
   precision=200;
   iopt=find_opt(argc,argv,"-prec");
   error(iopt==argc-1,1,"main [minmax.c]",syntax);
   if(iopt>0)
   {
      sscanf(argv[iopt+1],"%d",&precision);
      error(precision<64,1,"main [minmax.c]",
	    "invalid goal: the working precision is negative or too small");
   }
   rem.precision=precision;
   alloc_remes(&rem);

      
   rem.nstart=-1;
   iopt=find_opt(argc,argv,"-nstart");
   error(iopt==argc-3,1,"main [minmax.c]",syntax);
   if(iopt>0)
   {
      sscanf(argv[iopt+1],"%d",&n);
      error(n<1,1,"main [minmax.c]",
	    "nstart must be >1");

      sscanf(argv[iopt+2],"%lf",&alpha);
      error(alpha<0.0,1,"main [minmax.c]",
	    "alpha must be >0.0");

      sscanf(argv[iopt+3],"%lf",&beta);

      rem.nstart=n;
      mpfr_set_d(rem.alpha,alpha,ROUNDING);
      mpfr_set_d(rem.beta,beta,ROUNDING);
   }

   iopt=find_opt(argc,argv,"-p");
   error(iopt==0 || iopt==argc-1,1,"main [minmax.c]",syntax);	 
   sscanf(argv[iopt+1],"%d",&p);
   error(p==0,1,"main [minmax.c]",
	 "invalid power: the program approximates ( x + mu^2 )^{p/q} where p and q are integers with p!=0 and q>0");

   iopt=find_opt(argc,argv,"-q");
   error(iopt==0 || iopt==argc-1,1,"main [minmax.c]",syntax);	 
   sscanf(argv[iopt+1],"%d",&q);
   error(q<0,1,"main [minmax.c]",
	 "invalid power: the program approximates ( x + mu^2 )^{p/q} where p and q are integers with p!=0 and q>0");
   rem.p=p;
   rem.q=q;
   f=p/q;
   m=p%q;
   if(f!=0 && m==0)
   {
      fprintf(stderr,"\nthe approximation you are searching is x^{%d}\n\n",p/q);
      exit(0);
   }
   else if(f!=0 && m!=0)
   {
      s=1;
      if(p<0)
	 s=-1;
      if(q%m==0)
      {
	 q=s*(q/m);
	 p=s;
      }
      else
	 p=m;
      fprintf(stderr,"\nyou are trying to approximate ( x + mu^2 )^{%d} * ( x + mu^2 )^{%d/%d}, ",f,p,q);
      fprintf(stderr,"set p= %d, q= %d and try again\n\n",p,q);
      exit(0);
   }
   
      
   mu=0.0;
   iopt=find_opt(argc,argv,"-mu");
   error(iopt==argc-1,1,"main [minmax.c]",syntax);
   if(iopt>0)
      sscanf(argv[iopt+1],"%lf",&mu);
   mpfr_set_d(rem.mu,mu,ROUNDING);


   iopt=find_opt(argc,argv,"-ra");
   error(iopt==0 || iopt==argc-1,1,"main [minmax.c]",syntax);
   sscanf(argv[iopt+1],"%lf",&ra);
   error(ra*ra+mu*mu==0.0,1,"main [minmax.c]",
	 "invalid range: the program approximates ( x + mu^2 )^{p/q} for x in the range [ra^2,rb^2] with rb > ra and ra^2+mu^2 > 0");

   iopt=find_opt(argc,argv,"-rb");
   error(iopt==0 || iopt==argc-1,1,"main [minmax.c]",syntax);
   sscanf(argv[iopt+1],"%lf",&rb);
   error(rb<=ra,1,"main [minmax.c]",
	 "invalid range: the program approximates ( x + mu^2 )^{p/q} for x in the range [ra^2,rb^2] with rb > ra and ra^2+mu^2 > 0");
   mpfr_set_d(rem.ra,ra,ROUNDING);
   mpfr_set_d(rem.rb,rb,ROUNDING);


   iopt=find_opt(argc,argv,"-goal");
   error(iopt==0 || iopt==argc-1,1,"main [minmax.c]",syntax);
   sscanf(argv[iopt+1],"%lf",&goal);
   error(goal<=1.0e-30,1,"main [minmax.c]",
	 "invalid goal: the target error is negative or too small");
   mpfr_set_d(rem.goal,goal,ROUNDING);


   relflag=1;
   iopt=find_opt(argc,argv,"-abs");
   error(iopt==argc-1,1,"main [minmax.c]",syntax);
   if(iopt>0)
      relflag=0;
   rem.relflag=relflag;

   verbose=1;
   iopt=find_opt(argc,argv,"-nopen");
   if(iopt>0)
      verbose=0;
   rem.verbose=verbose;

   sprintf(rem.path,"./p%dq%dmu%.8era%.8erb%.8e",p,q,mu,ra,rb);
   iopt=find_opt(argc,argv,"-outpath");
   error(iopt==argc-1,1,"main [minmax.c]",syntax);
   if(iopt>0)
   {
      slen=strlen(argv[iopt+1]);
      error(slen>NAME_SIZE/2,1,"main [minmax.c]","path too long");
      sscanf(argv[iopt+1],"%s",rem.path);
      sprintf(rem.path,"%s/p%dq%dra%.8erb%.8e",rem.path,p,q,ra,rb);
   }
   sprintf(syscmd,"/bin/rm -rf %s",rem.path);      
   system(syscmd);
   sprintf(syscmd,"mkdir %s",rem.path);      
   system(syscmd);

   sprintf(logname,"%s/minmax.log",rem.path);         
   flog=freopen(logname,"w",stdout);
}


int main(int argc,char *argv[])
{

   parse_input(argc,argv);

   remes(&rem);

   fclose(flog);      
   
   exit(0);
}
