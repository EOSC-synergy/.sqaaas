#ifndef LIB_MAIN_H
#define LIB_MAIN_H


/* lib/dat_utils.c */
extern void write_dat(FILE *dat);
extern int read_dat(FILE *dat,int *nt);
extern void set_dat(int nt,int iac,double dH,double avpl3,double avpl1);
extern void print_dat(void);


/* lib/ms3_utils.c */
extern void write_ms3dat(FILE *dat);
extern int read_ms3dat(FILE *dat,int *nt);
extern void write_ms3dat_head(FILE *fdat);
extern void check_ms3dat_head(FILE *fdat);
extern void set_ms3dat(int nt);
extern void print_ms3dat(void);


/* lib/ms5_utils.c */
extern void write_ms5dat(FILE *dat);
extern int read_ms5dat(FILE *dat,int *nt);
extern void write_ms5dat_head(FILE *fdat);
extern void check_ms5dat_head(FILE *fdat);
extern void set_ms5dat(int nt);
extern void print_ms5dat(void);


/* lib/wflow_parms.c */
typedef struct
{
   int flint,dn,nn,tmax;
   double eps;
} wflow_parms_t;
extern wflow_parms_t wflow_parms(void);
extern void read_wflow_parms(void);
extern void check_wflow_parms(FILE *fpar);
extern void write_wflow_parms(FILE *fpar);
extern void print_wflow_parms(void);


#endif /* LIB_MAIN_H */
