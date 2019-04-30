
/*******************************************************************************
*
* File grid_events.h
*
* Copyright (C) 2011 Martin Luescher, 2016 Isabel Campos
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Block grid events modified to include QED (switching u/ud <-> h/hd)
*
*******************************************************************************/

#define GRID_EVENTS_H

#if (defined FLAGS_C)

static void (*grid_event_fcts[(int)(EVENTS)+1])(void)={NULL};

static void GridAssignedHd2hbgr(void)
{
   if ((*gf).shf&0x1)
      error_root(1,1,"GridAssignedHd2hbgr [grid_events.h]",
                 "Event involving shared fields");
   else
     {
      (*gf).h[0]=lat.hd[0];
      (*gf).h[1]=lat.hd[1];
     }
}

static void GridErasedHbgr(void)
{
   if ((*gf).shf&0x1)
      error_root(1,1,"GridErasedHbgr [grid_events.h]",
                 "Event involving shared fields");
   else
     {
      (*gf).h[0]=0;
      (*gf).h[1]=0;
     }
}

static void GridAssignedSwd2swbgr(void)
{
   if ((*gf).shf&0x1)
      error_root(1,1,"GridAssignedSwd2swbgr [grid_events.h]",
                 "Event involving shared fields");
   else
   {
      (*gf).sw[0]=lat.swd[0];
      (*gf).sw[1]=lat.swd[1];
      (*gf).sw[2]=lat.swd[2];      
      (*gf).sw[3]=lat.swd[3];      
   }
}

static void GridAssignedSwd2swdbgr(void)
{
   if ((*gf).shf&0x2)
      error_root(1,1,"GridAssignedSwd2swdbgr [grid_events.h]",
                 "Event involving shared fields");
   else
   {
      (*gf).swd[0]=lat.swd[0];
      (*gf).swd[1]=lat.swd[1];
      (*gf).swd[2]=lat.swd[2];      
      (*gf).swd[3]=lat.swd[3];      
   }
}

static void GridInvertedSwdE(void)
{
   if ((*gf).shf&0x2)
      error_root(1,1,"GridInvertedSwdE [grid_events.h]",
                 "Event involving shared fields");
   else
      (*gf).swd[1]^=0x1;
}

static void GridInvertedSwdO(void)
{
   if ((*gf).shf&0x2)
      error_root(1,1,"GridInvertedSwdO [grid_events.h]",
                 "Event involving shared fields");
   else
      (*gf).swd[2]^=0x1;
}

static void GridInvertedSwE(void)
{
   if ((*gf).shf&0x1)
      error_root(1,1,"GridInvertedSwE [grid_events.h]",
                 "Event involving shared fields");
   else
      (*gf).sw[1]^=0x1;
}

static void GridInvertedSwO(void)
{
   if ((*gf).shf&0x1)
      error_root(1,1,"GridInvertedSwO [grid_events.h]",
                 "Event involving shared fields");
   else
      (*gf).sw[2]^=0x1;
}

static void GridErasedSw(void)
{
   if ((*gf).shf&0x1)
      error_root(1,1,"GridErasedSw [grid_events.h]",
                 "Event involving shared fields");
   else
   {
      (*gf).sw[0]=0;
      (*gf).sw[1]=0;
      (*gf).sw[2]=0;
      (*gf).sw[3]=0;
   }
}

static void GridErasedSwd(void)
{
   if ((*gf).shf&0x2)
      error_root(1,1,"GridErasedSwd [grid_events.h]",
                 "Event involving shared fields");
   else
   {
      (*gf).swd[0]=0;
      (*gf).swd[1]=0;
      (*gf).swd[2]=0;
      (*gf).swd[3]=0;
   }
}

static void set_grid_events(void)
{
   grid_event_fcts[(int)(ASSIGNED_HD2HBGR)]=GridAssignedHd2hbgr;
   grid_event_fcts[(int)(ERASED_HBGR)]=GridErasedHbgr;
   grid_event_fcts[(int)(ASSIGNED_SWD2SWBGR)]=GridAssignedSwd2swbgr;
   grid_event_fcts[(int)(ASSIGNED_SWD2SWDBGR)]=GridAssignedSwd2swdbgr;
   grid_event_fcts[(int)(INVERTED_SWD_E)]=GridInvertedSwdE;   
   grid_event_fcts[(int)(INVERTED_SWD_O)]=GridInvertedSwdO;
   grid_event_fcts[(int)(INVERTED_SW_E)]=GridInvertedSwE;   
   grid_event_fcts[(int)(INVERTED_SW_O)]=GridInvertedSwO;   
   grid_event_fcts[(int)(ERASED_SW)]=GridErasedSw;
   grid_event_fcts[(int)(ERASED_SWD)]=GridErasedSwd;
}

#endif
