
/*******************************************************************************
* QED MODIFIED: Isabel Campos
*
* File flags/events.h
*
* Copyright (C) 2009, 2010, 2012, 2016 Martin Luescher, Isabel Campos
*               2017 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Full-lattice events, including QED 
*
*******************************************************************************/

#define EVENTS_H

#if (defined FLAGS_C)

static void (*event_fcts[(int)(EVENTS)+1])(void)={NULL};


static void LatUpdatedU(void)
{
   lat.u=next_tag();
}

static void LatUpdatedUd(void)
{
   lat.ud=next_tag();
}

static void LatAssignedUd2u(void)
{
   lat.u=lat.ud;
}

static void LatCopiedBndUd(void)
{
   lat.udbuf=lat.ud;
}

static void LatSetBstap(void)
{
   lat.bstap=lat.ud;
}

static void LatShiftedUd(void)
{
   lat.ud=next_tag();
   lat.udbuf=0;
}

static void LatComputedFts(void)
{
   lat.fts=lat.ud;
}

static void LatErasedSw(void)
{
   lat.sw[0]=0;
   lat.sw[1]=0;
   lat.sw[2]=0;
   lat.sw[3]=0;
}

static void LatErasedSwd(void)
{
   lat.swd[0]=0;
   lat.swd[1]=0;
   lat.swd[2]=0;
   lat.swd[3]=0;
}

static void LatComputedSwd(void)
{
   lat.swd[0]=lat.ud;
   lat.swd[1]=0;
   lat.swd[2]=0;
   lat.swd[3]=lat.ad;
}

static void LatAssignedSwd2sw(void)
{
   lat.sw[0]=lat.swd[0];
   lat.sw[1]=lat.swd[1];
   lat.sw[2]=lat.swd[2];
   lat.sw[3]=lat.swd[3];
}

static void LatInvertedSwdE(void)
{
   lat.swd[1]^=0x1;
}

static void LatInvertedSwdO(void)
{
   lat.swd[2]^=0x1;
}

static void LatErasedAw(void)
{
   lat.aw[0]=0;
   lat.aw[1]=0;
}

static void LatErasedAwhat(void)
{
   lat.awh[0]=0;
   lat.awh[1]=0;
}

static void LatComputedAw(void)
{
   lat.aw[0]=lat.ud;
   lat.aw[1]=lat.ad;
}

static void LatComputedAwhat(void)
{
   lat.awh[0]=lat.ud;
   lat.awh[1]=lat.ad;
}

/******************** NEW QED EVENTS TO CHECK *********************/

static void LatUpdatedAd(void)
{
   lat.ad=next_tag();
}

static void LatComputedU1d(void)
{
   lat.u1d=lat.ad;
}

static void LatComputedBndU1d(void)
{
   lat.u1dbuf=lat.u1d;
}

static void LatCopiedBndAd(void)
{
   lat.adbuf=lat.ad;
}

static void LatSetU1Bstap(void)
{
   lat.u1bstap=lat.ad;
}

static void LatComputedU1Fts(void)
{
   lat.u1fts=lat.ad;
}

static void LatComputedHd(void)
{
   lat.hd[0]=lat.ud;
   lat.hd[1]=lat.ad;
}

static void LatAssignedHd2h(void)
{
   lat.h[0]=lat.hd[0];
   lat.h[1]=lat.hd[1];
}

static void LatErasedHd(void)
{
   lat.hd[0]=0;
   lat.hd[1]=0;
}
static void LatErasedH(void)
{
   lat.h[0]=0;
   lat.h[1]=0;
}

static void set_events(void)
{
   event_fcts[(int)(UPDATED_U)]=LatUpdatedU;
   event_fcts[(int)(UPDATED_UD)]=LatUpdatedUd;
   event_fcts[(int)(ASSIGNED_UD2U)]=LatAssignedUd2u;
   event_fcts[(int)(COPIED_BND_UD)]=LatCopiedBndUd;
   event_fcts[(int)(SET_BSTAP)]=LatSetBstap;
   event_fcts[(int)(SHIFTED_UD)]=LatShiftedUd;
   event_fcts[(int)(COMPUTED_FTS)]=LatComputedFts;
   event_fcts[(int)(ERASED_SW)]=LatErasedSw;
   event_fcts[(int)(ERASED_SWD)]=LatErasedSwd;
   event_fcts[(int)(COMPUTED_SWD)]=LatComputedSwd;
   event_fcts[(int)(ASSIGNED_SWD2SW)]=LatAssignedSwd2sw;
   event_fcts[(int)(INVERTED_SWD_E)]=LatInvertedSwdE;
   event_fcts[(int)(INVERTED_SWD_O)]=LatInvertedSwdO;
   event_fcts[(int)(ERASED_AW)]=LatErasedAw;
   event_fcts[(int)(ERASED_AWHAT)]=LatErasedAwhat;
   event_fcts[(int)(COMPUTED_AW)]=LatComputedAw;
   event_fcts[(int)(COMPUTED_AWHAT)]=LatComputedAwhat;
   event_fcts[(int)(UPDATED_AD)]=LatUpdatedAd;
   event_fcts[(int)(COMPUTED_U1D)]=LatComputedU1d;
   event_fcts[(int)(COPIED_BND_AD)]=LatCopiedBndAd;
   event_fcts[(int)(COMPUTED_BND_U1D)]=LatComputedBndU1d;
   event_fcts[(int)(SET_U1BSTAP)]=LatSetU1Bstap;
   event_fcts[(int)(COMPUTED_U1FTS)]=LatComputedU1Fts;
   event_fcts[(int)(COMPUTED_HD)]=LatComputedHd;
   event_fcts[(int)(ERASED_HD)]=LatErasedHd;
   event_fcts[(int)(ASSIGNED_HD2H)]=LatAssignedHd2h;
   event_fcts[(int)(ERASED_H)]=LatErasedH;

}

#endif
