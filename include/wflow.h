/*******************************************************************************
*
* File wflow.h
*
* Copyright (C) 2009, 2010, 2011, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef WFLOW_H
#define WFLOW_H

/* WFLOW_C */
extern void fwd_su3_euler(int n,double eps);
extern void fwd_su3_rk2(int n,double eps);
extern void fwd_su3_rk3(int n,double eps);

/* WFLOW_U1_C */
extern void fwd_u1_euler(int n,double eps);
extern void fwd_u1_rk2(int n,double eps);
extern void fwd_u1_rk3(int n,double eps);

#endif
