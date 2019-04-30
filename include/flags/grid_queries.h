
/*******************************************************************************
*
* File grid_queries.h
*
* Copyright (C) 2011 Martin Luescher, 2016 Isabel Campos
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Block grid queries modified to include QED
*
*******************************************************************************/

#define GRID_QUERIES_H

#if (defined FLAGS_C)

static int (*grid_query_fcts[(int)(QUERIES)+1])(void)={NULL};

static int GridQueryHbgrUp2Date(void)
{
   if ((*gf).shf&0x1)
   {
      error_loc(1,1,"GridQueryHbgrMatchHd [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
      return (((*gf).h[0]==lat.ud)&&((*gf).h[1]==lat.ad));
}

static int GridQuerySwUp2date(void)
{
   if ((*gf).shf&0x1)
   {
      error_loc(1,1,"GridQuerySwUp2date [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
     return (((*gf).sw[0]==lat.ud)&&((*gf).sw[3]==lat.ad));
}

static int GridQuerySwEInverted(void)
{
   if ((*gf).shf&0x1)
   {
      error_loc(1,1,"GridQuerySwEInverted [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
      return ((*gf).sw[1]==1);
}

static int GridQuerySwOInverted(void)
{
   if ((*gf).shf&0x1)
   {
      error_loc(1,1,"GridQuerySwOInverted [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
      return ((*gf).sw[2]==1);
}

static int GridQuerySwdUp2date(void)
{
   if ((*gf).shf&0x2)
   {
      error_loc(1,1,"GridQuerySwdUp2date [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
     return (((*gf).swd[0]==lat.ud)&&((*gf).swd[3]==lat.ad));
}

static int GridQuerySwdEInverted(void)
{
   if ((*gf).shf&0x2)
   {
      error_loc(1,1,"GridQuerySwdEInverted [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
      return ((*gf).swd[1]==1);
}

static int GridQuerySwdOInverted(void)
{
   if ((*gf).shf&0x2)
   {
      error_loc(1,1,"GridQuerySwdOInverted [grid_queries.h]",
                "Query involving shared fields");
      return -1;
   }
   else
      return ((*gf).swd[2]==1);
}

static void set_grid_queries(void)
{
   grid_query_fcts[(int)(HBGR_UP2DATE)]=GridQueryHbgrUp2Date;
   grid_query_fcts[(int)(SW_UP2DATE)]=GridQuerySwUp2date;
   grid_query_fcts[(int)(SW_E_INVERTED)]=GridQuerySwEInverted;
   grid_query_fcts[(int)(SW_O_INVERTED)]=GridQuerySwOInverted;   
   grid_query_fcts[(int)(SWD_UP2DATE)]=GridQuerySwdUp2date;
   grid_query_fcts[(int)(SWD_E_INVERTED)]=GridQuerySwdEInverted;
   grid_query_fcts[(int)(SWD_O_INVERTED)]=GridQuerySwdOInverted;
}

#endif
