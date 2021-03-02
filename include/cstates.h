
/*******************************************************************************
*
* File cstates.h
*
* Copyright (C) 2017 Nazario Tantalo
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef CSTATES_H
#define CSTATES_H

/* CFACTOR_C */
extern double* get_cfactor_phase(int coulomb,int mu);
extern void mul_cfactor(int iflag,int coulomb,int mu,spinor_dble *pk,spinor_dble *pl);
extern void mul_cfactor_muaverage(int iflag,int coulomb,spinor_dble *pk,spinor_dble *pl);


/* CCOULOMB_C */
extern void nabla_sq_dvec(double *s, double *r);
extern void div_sym_dvec(double *s1, double *s2, double *s3, double *r);
extern double inv_nabla_sq(double *out, double *in);
extern double* get_coulomb_amu(void);

#endif
