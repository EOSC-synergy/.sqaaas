
/*******************************************************************************
*
* File flags.h
*
* Copyright (C) 2009-2014, 2016 Martin Luescher, Isabel Campos
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef FLAGS_H
#define FLAGS_H

#include <stdio.h>

#ifndef BLOCK_H
#include "block.h"
#endif

typedef enum
{
   UPDATED_U,UPDATED_UD,ASSIGNED_UD2U,
   COPIED_BND_UD,SET_BSTAP,SHIFTED_UD,COMPUTED_FTS,
   ERASED_SW,ERASED_SWD,COMPUTED_SWD,ASSIGNED_SWD2SW,
   INVERTED_SW_E,INVERTED_SW_O,
   INVERTED_SWD_E,INVERTED_SWD_O,
   ASSIGNED_SWD2SWBGR,ASSIGNED_SWD2SWDBGR,
   ERASED_AW,ERASED_AWHAT,COMPUTED_AW,COMPUTED_AWHAT,
   UPDATED_AD, COMPUTED_U1D,
   COMPUTED_BND_U1D, COPIED_BND_AD,COMPUTED_U1FTS, SET_U1BSTAP,
   COMPUTED_HD, ASSIGNED_HD2H, ERASED_H,ERASED_HD,
   ASSIGNED_HD2HBGR, ERASED_HBGR,
   EVENTS
} event_t;

typedef enum
{
   U_MATCH_UD,UDBUF_UP2DATE,BSTAP_UP2DATE,
   FTS_UP2DATE,
   SW_UP2DATE,SW_E_INVERTED,SW_O_INVERTED,
   SWD_UP2DATE,SWD_E_INVERTED,SWD_O_INVERTED,
   AW_UP2DATE,AWHAT_UP2DATE,
   U1D_UP2DATE, ADBUF_UP2DATE, U1DBUF_UP2DATE,
   U1FTS_UP2DATE, U1BSTAP_UP2DATE,
   HD_UP2DATE,H_UP2DATE,HBGR_UP2DATE,
   QUERIES
} query_t;

typedef enum
{
   ACG_SU3,ACF_TM1,ACF_TM1_EO,ACF_TM1_EO_SDET,
   ACF_TM2,ACF_TM2_EO,ACF_RAT,ACF_RAT_SDET,ACG_U1,
   ACTIONS
} action_t;

typedef enum
{
   LPFR,OMF2,OMF4,
   INTEGRATORS
} integrator_t;

typedef enum
{
   FRG_SU3,FRF_TM1,FRF_TM1_EO,FRF_TM1_EO_SDET,
   FRF_TM2,FRF_TM2_EO,FRF_RAT,FRF_RAT_SDET,FRG_U1,
   FORCES
} force_t;

typedef enum
{
   RWTM1,RWTM1_EO,RWTM2,RWTM2_EO,RWRAT,
   RWFACTS
} rwfact_t;

typedef enum
{
   CGNE,MSCG,SAP_GCR,DFL_SAP_GCR,
   SOLVERS
} solver_t;

typedef struct
{
   action_t action;
   int ipf,im0;
   int irat[3],imu[4];
   int isp[4];
} action_parms_t;

typedef struct
{
   int type;
   int SFtype;
   int cstar;
   double phi[2][3];
   double ad[2];
} bc_parms_t;

typedef struct
{
   int qhat;
   double m0,su3csw,u1csw,cF[2],theta[3];
} dirac_parms_t;

typedef struct
{
   int bs[4];
   int Ns;
} dfl_parms_t;

typedef struct
{
   int nkv,nmx;
   double res;
} dfl_pro_parms_t;

typedef struct
{
   int ninv,nmr,ncy;
   double kappa,mu;
   dirac_parms_t dp;
} dfl_gen_parms_t;

typedef struct
{
   int nsm;
   double dtau;
} dfl_upd_parms_t;

typedef struct
{
   int gauge;
   int nfl;
} flds_parms_t;

typedef struct
{
   force_t force;
   int ipf,im0;
   int irat[3],imu[4];
   int isp[4];
   int ncr[4],icr[4];
} force_parms_t;

typedef struct
{
   int npf,nlv;
   int nact,nmu;
   int *iact;
   double tau,*mu;
} hmc_parms_t;

typedef struct
{
   integrator_t integrator;
   double lambda;
   int nstep,nfr;
   int *ifr;
} mdint_parms_t;

typedef struct
{
   int power[2];
   int degree;
   double range[2];
   double mu0,A,delta;
   double *nu,*mu;
} rat_parms_t;

typedef struct
{
   rwfact_t rwfact;
   int im0,nsrc;
   int irp,nfct;
   double *mu;
   int *np,*isp;
} rw_parms_t;

typedef struct
{
   int bs[4];
   int isolv;
   int nmr,ncy;
} sap_parms_t;

typedef struct
{
   solver_t solver;
   int nmx,nkv;
   int isolv,nmr,ncy;
   double res;
} solver_parms_t;

typedef struct
{
   double beta,c0,c1,cG[2];
} su3lat_parms_t;

typedef struct
{
   int eoflg;
} tm_parms_t;

typedef struct
{
   int type;
   double alpha,invqel,beta,lambda,c0,c1,cG[2];
} u1lat_parms_t;

typedef struct
{
   int n;
   double eps;
} wflow_parms_t;

/* FLAGS_C */
extern void set_flags(event_t event);
extern void set_grid_flags(blk_grid_t grid,event_t event);
extern int query_flags(query_t query);
extern int query_grid_flags(blk_grid_t grid,query_t query);
extern void print_flags(void);
extern void print_grid_flags(blk_grid_t grid);

/* ACTION_PARMS_C */
extern action_parms_t set_action_parms(int iact,action_t action,int ipf,
                                       int im0,int *irat,int *imu,int *isp);
extern action_parms_t action_parms(int iact);
extern void read_action_parms(int iact);
extern void print_action_parms(void);
extern void write_action_parms(FILE *fdat);
extern void check_action_parms(FILE *fdat);

/* DIRAC_PARMS_C */
extern dirac_parms_t set_dirac_parms9(int qhat,double m0,
                                      double su3csw,double u1csw,
                                      double cF,double cF_prime,
                                      double th1,double th2,double th3);
extern dirac_parms_t set_dirac_parms1(dirac_parms_t *par);
extern void print_dirac_parms(void);
extern dirac_parms_t dirac_parms(void);
extern tm_parms_t set_tm_parms(int eoflg);
extern tm_parms_t tm_parms(void);

/* DFL_PARMS_C */
extern dfl_parms_t set_dfl_parms(int *bs,int Ns);
extern dfl_parms_t dfl_parms(void);
extern dfl_pro_parms_t set_dfl_pro_parms(int nkv,int nmx,double res);
extern dfl_pro_parms_t dfl_pro_parms(void);
extern dfl_gen_parms_t set_dfl_gen_parms(double kappa,double mu,
                                    int qhat,double su3csw,double u1csw,
                                    double cF,double cF_prime,
                                    double th1,double th2,double th3,
                                    int ninv,int nmr,int ncy);
extern dfl_gen_parms_t dfl_gen_parms(void);
extern dfl_upd_parms_t set_dfl_upd_parms(double dtau,int nsm);
extern dfl_upd_parms_t dfl_upd_parms(void);
extern void print_dfl_parms(int ipr);
extern void write_dfl_parms(FILE *fdat);
extern void check_dfl_parms(FILE *fdat);

/* FORCE_PARMS_C */
extern force_parms_t set_force_parms(int ifr,force_t force,int ipf,int im0,
                                     int *irat,int *imu,int *isp,int *ncr);
extern force_parms_t force_parms(int ifr);
extern void read_force_parms(int ifr);
extern void read_force_parms2(int ifr);
extern void print_force_parms(void);
extern void print_force_parms2(void);
extern void write_force_parms(FILE *fdat);
extern void check_force_parms(FILE *fdat);

/* HMC_PARMS_C */
extern hmc_parms_t set_hmc_parms(int nact,int *iact,int npf,
                                 int nmu,double *mu,int nlv,double tau);
extern hmc_parms_t hmc_parms(void);
extern void print_hmc_parms(void);
extern void write_hmc_parms(FILE *fdat);
extern void check_hmc_parms(FILE *fdat);

/* LAT_PARMS_C */
extern flds_parms_t set_flds_parms(int gauge,int nfl);
extern bc_parms_t set_bc_parms(int type,int SFtype,int cstar,
                               double *phi,double *phi_prime);
extern su3lat_parms_t set_su3lat_parms(double beta,double c0,
                      double cG,double cG_prime);
extern u1lat_parms_t set_u1lat_parms(int type,double alpha,double invqel,
                           double lambda,double c0,double cG,double cG_prime);
extern dirac_parms_t set_qlat_parms(int ifl,double kappa,int qhat,
                     double su3csw,double u1csw,double cF,double cF_prime,
                     double th1,double th2,double th3);
extern flds_parms_t flds_parms(void);
extern int gauge(void);
extern bc_parms_t bc_parms(void);
extern int bc_type(void);
extern int bc_cstar(void);
extern su3lat_parms_t su3lat_parms(void);
extern u1lat_parms_t u1lat_parms(void);
extern dirac_parms_t qlat_parms(int ifl);
extern void print_flds_parms(void);
extern void print_bc_parms(void);
extern void print_lat_parms(void);
extern void print_flds_bc_lat_parms(void);
extern void write_flds_bc_lat_parms(FILE *fdat);
extern void check_flds_bc_lat_parms(FILE *fdat);

/* MDINT_PARMS_C */
extern mdint_parms_t set_mdint_parms(int ilv,
                                     integrator_t integrator,double lambda,
                                     int nstep,int nfr,int *ifr);
extern mdint_parms_t mdint_parms(int ilv);
extern void read_mdint_parms(int ilv);
extern void print_mdint_parms(void);
extern void write_mdint_parms(FILE *fdat);
extern void check_mdint_parms(FILE *fdat);

/* RAT_PARMS_C */
extern rat_parms_t set_rat_parms(int irp,int num,int den,int degree,
                                 double *range,double mu0,double A,double *nu,
                                 double *mu,double delta,double *x);
extern rat_parms_t rat_parms(int irp);
extern void read_rat_parms(int irp);
extern void print_rat_parms(void);
extern void write_rat_parms(FILE *fdat);
extern void check_rat_parms(FILE *fdat);

/* RW_PARMS_C */
extern rw_parms_t set_rw_parms(int irw,rwfact_t rwfact,int im0,int nsrc,
                               int irp,int nfct,double *mu,int *np,int *isp);
extern rw_parms_t rw_parms(int irw);
extern void read_rw_parms(int irw);
extern void print_rw_parms(void);
extern void write_rw_parms(FILE *fdat);
extern void check_rw_parms(FILE *fdat);

/* SAP_PARMS_C */
extern sap_parms_t set_sap_parms(int *bs,int isolv,int nmr,int ncy);
extern sap_parms_t sap_parms(void);
extern void print_sap_parms(int ipr);
extern void write_sap_parms(FILE *fdat);
extern void check_sap_parms(FILE *fdat);

/* SOLVER_PARMS_C */
extern solver_parms_t set_solver_parms(int isp,solver_t solver,
                                       int nkv,int isolv,int nmr,int ncy,
                                       int nmx,double res);
extern solver_parms_t solver_parms(int isp);
extern void read_solver_parms(int isp);
extern void print_solver_parms(int *isap,int *idfl);
extern void write_solver_parms(FILE *fdat);
extern void check_solver_parms(FILE *fdat);

/* WFLOW_PARMS_C */
extern wflow_parms_t set_wflow_parms(int n,double eps);
extern wflow_parms_t wflow_parms(void);

/* PARMS_UTILS */
extern void check_global_dble(const char fnm[256],int nargs,...);
extern void check_global_dblearray(const char fnm[256],int size,double *data);
extern void check_global_int(const char fnm[256],int nargs,...);
extern void check_global_intarray(const char fnm[256],int size,int *data);
extern void write_little_dble(FILE *fdat,int nargs,...);
extern void write_little_int(FILE *fdat,int nargs,...);
extern void write_little_dblearray(FILE *fdat,int size,double *data);
extern void write_little_intarray(FILE *fdat,int size,int *data);
extern void check_fpar_dble(const char fnm[256],FILE *fdat,int nargs,...);
extern void check_fpar_int(const char fnm[256],FILE *fdat,int nargs,...);
extern void check_fpar_dblearray(const char fnm[256],FILE *fdat,int size,double *data);
extern void check_fpar_intarray(const char fnm[256],FILE *fdat,int size,int *data);


#endif
