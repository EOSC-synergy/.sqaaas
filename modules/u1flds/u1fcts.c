
/*******************************************************************************
*
* File u1fcts.c
*
* Copyright (C) 2016 Marina Marinkovic
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Various functions for complex-number multiplication.
*
* The externally accessible functions are
*
*  void u1xu1(complex_dble *u,complex_dble *v,complex_dble *w)
*    Computes w=u*v assuming that w is different from u.
*
*  void u1dagxu1(complex_dble *u,complex_dble *v,complex_dble *w)
*    Computes w=u^dag*v assuming that w is different from u.
*
*  void u1xu1dag(complex_dble *u,complex_dble *v,complex_dble *w)
*    Computes w=u*v^dag assuming that w is different from u.
*
*  void u1dagxu1dag(complex_dble *u,complex_dble *v,complex_dble *w)
*    Computes w=u^dag*v^dag assuming that w is different from u.
*
*  void cm1x1_retr(complex_dble *u,complex_dble *v,double *tr)
*    Assigns the real part of (*u)*(*v) to (*tr)
*
*  void cm1x1_imtr(complex_dble *u,complex_dble *v,double *tr)
*    Assigns the imaginary part of (*u)*(*v) to (*tr)
*
*
* Notes:
*
* The programs from u1fcts module do not perform any communications and can be
* called locally.
*
*******************************************************************************/

#define U1FCTS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "su3.h"


static complex_dble psi;


void u1xu1(complex_dble *u,complex_dble *v,complex_dble *w)
{ 
   psi.re=(*v).re;
   psi.im=(*v).im;

   (*w).re=(*u).re*psi.re-(*u).im*psi.im; 
   (*w).im=(*u).re*psi.im+(*u).im*psi.re; 
}

void u1dagxu1(complex_dble *u,complex_dble *v,complex_dble *w)
{
    psi.re=(*v).re;
    psi.im=(*v).im;
    
    (*w).re= (*u).re*psi.re+(*u).im*psi.im;  
    (*w).im= (*u).re*psi.im-(*u).im*psi.re;
}

void u1xu1dag(complex_dble *u,complex_dble *v,complex_dble *w)
{
    psi.re= (*v).re;
    psi.im=-(*v).im;
    
    (*w).re=(*u).re*psi.re-(*u).im*psi.im; 
    (*w).im=(*u).re*psi.im+(*u).im*psi.re; 
}

void u1dagxu1dag(complex_dble *u,complex_dble *v,complex_dble *w)
{
    psi.re= (*v).re;
    psi.im=-(*v).im;

    (*w).re= (*u).re*psi.re+(*u).im*psi.im;
    (*w).im= (*u).re*psi.im-(*u).im*psi.re;
}

void cm1x1_retr(complex_dble *u,complex_dble *v,double *tr)
{
    (*tr) =(*u).re*(*v).re-(*u).im*(*v).im;
}

void cm1x1_imtr(complex_dble *u,complex_dble *v,double *tr)
{
    (*tr) =(*u).re*(*v).im+(*u).im*(*v).re;
}
