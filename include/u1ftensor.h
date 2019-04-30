
/*******************************************************************************
*
* File u1ftensor.h
*
* Copyright (C) 2016 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef U1FTENSOR_H
#define U1FTENSOR_H

/* U1FTCOM_C */
extern void copy_bnd_u1ft(int n,double *ft);
extern void add_bnd_u1ft(int n,double *ft);

/* U1FTENSOR_C */
extern double **u1ftensor(void);

/* U1FLUXES_C */
extern double u1fluxes(int munu);
extern double u1fluxes_slices(int munu,double *Phisl);

/* MXW_ACTION_C */
extern double mxw_action(void);
extern double mxw_action_slices(double *asl);

#endif
