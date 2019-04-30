
/*******************************************************************************
*
* File flags/queries.h
*
* Copyright (C) 2009-2012, 2016 Martin Luescher, Isabel Campos
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Query descriptions including QED
*
*******************************************************************************/

#define QUERIES_H

#if (defined FLAGS_C)

static int (*query_fcts[(int)(QUERIES)+1])(void)={NULL};

static int QueryUMatchUd(void)
{
   return (lat.u==lat.ud);
}

static int QueryUdbufUp2date(void)
{
   return ((lat.ud>0)&&(lat.udbuf==lat.ud));
}

static int QueryBstapUp2date(void)
{
   return ((lat.ud>0)&&(lat.bstap==lat.ud));
}

static int QueryFtsUp2date(void)
{
   return ((lat.ud>0)&&(lat.fts==lat.ud));
}

static int QuerySwUp2date(void)
{
  return (((lat.ad>0)||(lat.ud>0))&&(lat.sw[0]==lat.ud)&&(lat.sw[3]==lat.ad));
}

static int QuerySwEInverted(void)
{
   return (lat.sw[1]==1);
}

static int QuerySwOInverted(void)
{
   return (lat.sw[2]==1);
}

static int QuerySwdUp2date(void)
{
  return (((lat.ad>0)||(lat.ud>0))&&(lat.swd[0]==lat.ud)&&(lat.swd[3]==lat.ad));
}

static int QuerySwdEInverted(void)
{
   return (lat.swd[1]==1);
}

static int QuerySwdOInverted(void)
{
   return (lat.swd[2]==1);
}

static int QueryAwUp2date(void)
{
  return (((lat.ad>0)||(lat.ud>0))&&(lat.aw[0]==lat.ud)&&(lat.aw[1]==lat.ad));
}

static int QueryAwhatUp2date(void)
{
  return (((lat.ad>0)||(lat.ud>0))&&(lat.awh[0]==lat.ud)&&(lat.awh[1]==lat.ad));
}

/************ QED queries *******/

static int QueryU1dUp2date(void)
{
   return ((lat.ad>0)&&(lat.u1d==lat.ad));
}

static int QueryAdbufUp2date(void)
{
   return ((lat.ad>0)&&(lat.adbuf==lat.ad));
}

static int QueryU1dbufUp2date(void)
{
   return ((lat.ad>0)&&(lat.u1dbuf==lat.ad));
}

static int QueryU1BstapUp2date(void)
{
  return ((lat.ad>0)&&(lat.u1bstap==lat.ad));
}

static int QueryU1FtsUp2date(void)
{
  return ((lat.ad>0)&&(lat.u1fts==lat.ad));
}

static int QueryHdUp2date(void)
{
   return (((lat.ad>0)||(lat.ud>0))&&(lat.hd[0]==lat.ud)&&(lat.hd[1]==lat.ad));
}

static int QueryHUp2date(void)
{
  return (((lat.ad>0)||(lat.ud>0))&&(lat.h[0]==lat.ud)&&(lat.h[1]==lat.ad));
}


static void set_queries(void)
{
   query_fcts[(int)(U_MATCH_UD)]=QueryUMatchUd;
   query_fcts[(int)(UDBUF_UP2DATE)]=QueryUdbufUp2date;
   query_fcts[(int)(BSTAP_UP2DATE)]=QueryBstapUp2date;
   query_fcts[(int)(FTS_UP2DATE)]=QueryFtsUp2date;
   query_fcts[(int)(SW_UP2DATE)]=QuerySwUp2date;
   query_fcts[(int)(SW_E_INVERTED)]=QuerySwEInverted;
   query_fcts[(int)(SW_O_INVERTED)]=QuerySwOInverted;
   query_fcts[(int)(SWD_UP2DATE)]=QuerySwdUp2date;
   query_fcts[(int)(SWD_E_INVERTED)]=QuerySwdEInverted;
   query_fcts[(int)(SWD_O_INVERTED)]=QuerySwdOInverted;
   query_fcts[(int)(AW_UP2DATE)]=QueryAwUp2date;
   query_fcts[(int)(AWHAT_UP2DATE)]=QueryAwhatUp2date;
   query_fcts[(int)(U1D_UP2DATE)]=QueryU1dUp2date;
   query_fcts[(int)(ADBUF_UP2DATE)]=QueryAdbufUp2date;
   query_fcts[(int)(U1DBUF_UP2DATE)]=QueryU1dbufUp2date;
   query_fcts[(int)(U1BSTAP_UP2DATE)]=QueryU1BstapUp2date;
   query_fcts[(int)(U1FTS_UP2DATE)]=QueryU1FtsUp2date;
   query_fcts[(int)(HD_UP2DATE)]=QueryHdUp2date;
   query_fcts[(int)(H_UP2DATE)]=QueryHUp2date;

}

#endif
