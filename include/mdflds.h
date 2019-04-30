
/*******************************************************************************
*
* File mdflds.h
*
* Copyright (C) 2017 Nazario Tantalo
* Copyright (C) 2016 Agostino Patella
* Copyright (C) 2011, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef MDFLDS_H
#define MDFLDS_H

#ifndef SU3_H
#include "su3.h"
#endif

typedef struct
{
   int npf;
   su3_alg_dble *su3mom,*su3frc;
   double *u1mom,*u1frc;
   spinor_dble **pf;
} mdflds_t;

/* FCOM_C */
extern void copy_bnd_su3frc(void);
extern void add_bnd_su3frc(void);
extern void copy_bnd_u1frc(void);
extern void add_bnd_u1frc(void);

/* MDFLDS_C */
extern mdflds_t *mdflds(void);
extern void set_su3frc2zero(void);
extern void bnd_su3mom2zero(void);
extern void random_su3mom(void);
extern double su3momentum_action(int icom);
extern void set_u1frc2zero(void);
extern void bnd_u1mom2zero(void);
extern void random_u1mom(void);
extern double u1momentum_action(int icom);

/* MD_CSTAR_C */
extern void orbi_cpy_su3mom(void);
extern void orbi_cpy_u1mom(void);

/* U1MOM_MAP_C */
extern void gather_u1mom(double* inmom,complex_dble **mom);
extern void scatter_u1mom(complex_dble **mom,double* outmom);

/* U1MOM_FACC_C */
extern void u1mom_Delta_no0(int inv,double* inmom,double* outmom);

#endif
