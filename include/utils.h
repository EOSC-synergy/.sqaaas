
/*******************************************************************************
*
* File utils.h
*
* Copyright (C) 2011, 2016 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include <limits.h>
#include <float.h>

#ifndef SU3_H
#include "su3.h"
#endif

#if ((DBL_MANT_DIG!=53)||(DBL_MIN_EXP!=-1021)||(DBL_MAX_EXP!=1024))
#error : Machine is not compliant with the IEEE-754 standard
#endif

#if (SHRT_MAX==0x7fffffff)
typedef short int stdint_t;
typedef unsigned short int stduint_t;
#elif (INT_MAX==0x7fffffff)
typedef int stdint_t;
typedef unsigned int stduint_t;
#elif (LONG_MAX==0x7fffffff)
typedef long int stdint_t;
typedef unsigned long int stduint_t;
#else
#error : There is no four-byte integer type on this machine
#endif

#undef UNKNOWN_ENDIAN
#undef LITTLE_ENDIAN
#undef BIG_ENDIAN

#define UNKNOWN_ENDIAN 0
#define LITTLE_ENDIAN 1
#define BIG_ENDIAN 2

#undef IMAX
#define IMAX(n,m) ((n)+((m)-(n))*((m)>(n)))

typedef enum
{
   ALL_PTS,EVEN_PTS,ODD_PTS,NO_PTS,PT_SETS
} ptset_t;

/* ENDIAN_C */
extern int endianness(void);
extern void bswap_int(int n,void *a);
extern void bswap_double(int n,void *a);

/* ERROR_C */
extern void set_error_file(char *path,int loc_flag);
extern void error(int test,int no,char *name,char *format,...);
extern void error_root(int test,int no,char *name,char *format,...);
extern void error_loc(int test,int no,char *name,char *format,...);

/* HSUM_C */
extern int init_hsum(int n);
extern void reset_hsum(int id);
extern void add_to_hsum(int id,double *x);
extern void local_hsum(int id,double *sx);
extern void global_hsum(int id,double *sx);

/* IOUTILS_C */
extern int write_little_dble(int blocking,FILE *fdat,int nargs,...);
extern int write_little_int(int blocking,FILE *fdat,int nargs,...);
extern int write_little_dblearray(int blocking,FILE *fdat,int size,double *data);
extern int write_little_intarray(int blocking,FILE *fdat,int size,int *data);
extern int read_little_dble(int blocking,FILE *fdat,int nargs,...);
extern int read_little_int(int blocking,FILE *fdat,int nargs,...);
extern int read_little_dblearray(int blocking,FILE *fdat,int size,double *data);
extern int read_little_intarray(int blocking,FILE *fdat,int size,int *data);
extern void check_little_dble(const char fnm[256],FILE *fdat,int nargs,...);
extern void check_little_int(const char fnm[256],FILE *fdat,int nargs,...);
extern void check_little_dblearray(const char fnm[256],FILE *fdat,int size,double *data);
extern void check_little_intarray(const char fnm[256],FILE *fdat,int size,int *data);

/* MUTILS_C */
extern int find_opt(int argc,char *argv[],char *opt);
extern int fdigits(double x);
extern void check_dir(char* dir);
extern void check_dir_root(char* dir);
extern int name_size(char *format,...);
extern long find_section(char *title);
extern long read_line(char *tag,char *format,...);
extern int count_tokens(char *tag);
extern void read_iprms(char *tag,int n,int *iprms);
extern void read_dprms(char *tag,int n,double *dprms);
extern void copy_file(char *in,char *out);

/* UTILS_C */
extern int safe_mod(int x,int y);
extern void *amalloc(size_t size,int p);
extern void afree(void *addr);
extern int mpi_permanent_tag(void);
extern int mpi_tag(void);
extern void message(char *format,...);
extern void check_global_dble(const char fnm[256],int nargs,...);
extern void check_global_dblearray(const char fnm[256],int size,double *data);
extern void check_global_int(const char fnm[256],int nargs,...);
extern void check_global_intarray(const char fnm[256],int size,int *data);

/* WSPACE_C */
extern void alloc_wud(int n);
extern su3_dble **reserve_wud(int n);
extern int release_wud(void);
extern int wud_size(void);
extern void alloc_wad(int n);
extern double **reserve_wad(int n);
extern int release_wad(void);
extern int wad_size(void);
extern void alloc_wf3d(int n);
extern su3_alg_dble **reserve_wf3d(int n);
extern int release_wf3d(void);
extern int wf3d_size(void);
extern void alloc_wf1d(int n);
extern double **reserve_wf1d(int n);
extern int release_wf1d(void);
extern int wf1d_size(void);
extern void alloc_ws(int n);
extern spinor **reserve_ws(int n);
extern int release_ws(void);
extern int ws_size(void);
extern void alloc_wsd(int n);
extern spinor_dble **reserve_wsd(int n);
extern int release_wsd(void);
extern int wsd_size(void);
extern void alloc_wv(int n);
extern complex **reserve_wv(int n);
extern int release_wv(void);
extern int wv_size(void);
extern void alloc_wvd(int n);
extern complex_dble **reserve_wvd(int n);
extern int release_wvd(void);
extern int wvd_size(void);

#endif
