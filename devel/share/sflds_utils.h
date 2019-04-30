#ifndef SFLDS_UTILS_H

#define SFLDS_UTILS_H

#ifdef SFLDS_UTILS_C
   #define EXTERN
#else /* SFLDS_UTILS_C */
   #define EXTERN extern
#endif /* SFLDS_UTILS_C */


EXTERN void transform_sd(spinor_dble *pk,spinor_dble *pl);


#undef EXTERN

#endif /* SFLDS_UTILS_H */
