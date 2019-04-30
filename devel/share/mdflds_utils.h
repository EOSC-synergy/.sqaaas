#ifndef MDFLDS_UTILS_H

#define MDFLDS_UTILS_H

#ifdef MDFLDS_UTILS_C
   #define EXTERN
#else /* MDFLDS_UTILS_C */
   #define EXTERN extern
#endif /* MDFLDS_UTILS_C */


/* UPDATE */
EXTERN void rot_ud(double eps);
EXTERN void rot_ad(double eps);


/* CHECKS */
EXTERN int check_bnd_su3frc(su3_alg_dble *frc);
EXTERN int check_zero_su3frc(su3_alg_dble *frc);
EXTERN int check_bnd_u1frc(double *frc);
EXTERN int check_zero_u1frc(double *frc);


#undef EXTERN

#endif /* MDFLDS_UTILS_H */
