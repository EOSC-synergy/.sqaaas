
/*******************************************************************************
*
* File remes.h
*
* Copyright (C) 2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef REMES_H
#define REMES_H

#include "gmp.h"
#include "mpfr.h"
#include "mpf2mpfr.h"

#define ROUNDING MPFR_RNDN

typedef struct
{
   char path[NAME_SIZE];
   int p,q,nfit,nstart;
   int relflag,sflag,verbose;
   int n,d,nd,nd2,nx;
   int precision;
   mpfr_t sdx,alpha,beta;
   mpfr_t eeps,eeps2;
   mpfr_t reeps,reeps2;
   mpfr_t ra,rb,mu,mu2,a,b,bpa,bma;
   mpfr_t ee,emax,emin,eprev;
   mpfr_t *cn,*cd;
   mpfr_t *x,*y,*dx;
   mpfr_t *e,*se;
   mpfr_t *rn,*rd;
   mpfr_t power,goal;
} remes_t;


/* REMES_C */
extern void alloc_remes(remes_t *rem);
extern void remes(remes_t *rem);


#endif
