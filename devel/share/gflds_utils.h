#ifndef GFLDS_UTILS_H

#define GFLDS_UTILS_H

#ifdef GFLDS_UTILS_C
   #define EXTERN
#else /* GFLDS_UTILS_C */
   #define EXTERN extern
#endif /* GFLDS_UTILS_C */


/* INITIALIZATION */
EXTERN void unit_gflds(void);
EXTERN void random_gflds(void);


/* TRANSLATION */
EXTERN int shift_gflds(int *s);


/* GAUGE TRANSFORMATIONS */
EXTERN su3_dble *g3tr(void);
EXTERN double *g1tr(void);
EXTERN void random_g(void);
EXTERN void transform_gflds(void);


#undef EXTERN

#endif /* GFLDS_UTILS_H */
