
/*******************************************************************************
*
* File u1flds.h
*
* Copyright (C) 2016 Marina Marinkovic
                2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef U1FLDS_H
#define U1FLDS_H


/* ADCOM_C */
extern void copy_bnd_ad(void);

/* ADFLDS_C */
extern double *adfld(void);
extern void random_ad(void);
extern void renormalize_ad(void);

/* AD_BCNDS_C */
extern void set_ad_bc(void);
extern int check_ad_bc(double tol);

/* AD_CSTAR_C */
extern void orbi_cpy_ad(void);

/* AD_SHIFT_C */
extern int shift_ad(int *s);

/* U1FCTS_C */
extern void u1xu1(complex_dble *u,complex_dble *v,complex_dble *w);
extern void u1dagxu1(complex_dble *u,complex_dble *v,complex_dble *w);
extern void u1xu1dag(complex_dble *u,complex_dble *v,complex_dble *w);
extern void u1dagxu1dag(complex_dble *u,complex_dble *v,complex_dble *w);
extern void cm1x1_retr(complex_dble *u,complex_dble *v,double *tr);
extern void cm1x1_imtr(complex_dble *u,complex_dble *v,double *tr);

/* U1FLDS_C */
typedef enum
{
   LOC=0,
   EXT
} u1lat_t;
extern complex_dble *u1dfld(u1lat_t lat);

/* U1_BSTAP_C */
extern complex_dble *u1_bstap(void);
extern void set_u1_bstap(void);

/* U1_PLAQ_SUM_C */
extern double u1_plaq_sum_dble(int icom);
extern double u1_plaq_wsum_dble(int icom);
extern double u1_plaq_action_slices(double *asl);

#endif
