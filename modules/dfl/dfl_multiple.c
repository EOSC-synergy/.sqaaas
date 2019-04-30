
/*******************************************************************************
*
* File dfl_mutiple.c
*
* Copyright (C) 2019 Agostino Patella
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Utilities to manage multiple deflation spaces.
*
* The externally accessible functions are
*
*   void use_dfl_subspace(int idfl)
*     It saves the currently used (or active) deflation subspace, and recalls
*     the deflation subspace with index idfl, which becomes the new active
*     deflation subspace. The deflation subspace is defined by the single- and
*     double-precision spinors in the DFL_BLOCKS grid, and by the vector fields
*     referenced by vflds() and vdflds().
* 
* The programs in this module perform global operations and must be called
* simultaneously on all MPI processes. The required workspaces are
*
*******************************************************************************/

#define DFL_MULTIPLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "flags.h"
#include "utils.h"
#include "block.h"
#include "dfl.h"
#include "vflds.h"
#include "global.h"

static int dflmax=-1,Ns,vol,nb,nv,adfl=-1;
static int vs_size, vds_size, sb_size;
static block_t *b;
static complex **vs=NULL;
static complex_dble **vds=NULL;
static spinor **sb=NULL;
static spinor_dble **sdb=NULL;

#ifdef DFL_MODES_DBG
   int my_rank;
#endif


static void alloc_storage(void)
{
   complex *vsptr=NULL;
   complex_dble *vdsptr=NULL;
   spinor *sbptr=NULL;
   spinor_dble *sdbptr=NULL;
   int ndfl,idfl;
   int isw;
   int *bs;
   dfl_parms_t dfl;
   dflst_t status;

#ifdef DFL_MODES_DBG
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#endif

   error_root(sizeof(complex)!=(2*sizeof(float)),1,
              "alloc_storage [dfl_multiple.c]",
              "The complex structures are not properly packed");
   error_root(sizeof(complex_dble)!=(2*sizeof(double)),1,
              "alloc_storage [dfl_multiple.c]",
              "The complex_dble structures are not properly packed");

   dfl=dfl_parms();
   bs=dfl.bs;
   Ns=dfl.Ns;

   dflmax=0;
   ndfl=0;
   while(1)
   {
      status=dfl_gen_parms(dflmax).status;
      if(status==DFL_OUTOFRANGE) break;
      if(status==DFL_DEF) ndfl++;
      dflmax++;
   }

   if(ndfl==0)
      return;

#ifdef DFL_MODES_DBG
   if (my_rank==0)
   {
      printf("dfl_multiple.c: Allocate %d pointers and %d subspaces\n",dflmax,ndfl);
      fflush(stdout);
   }
#endif

   error_root(dfl.Ns==0,1,"alloc_storage [dfl_multiple.c]",
         "The deflation subspace parameters are not set");
   
   vol=bs[0]*bs[1]*bs[2]*bs[3];
   nb=VOLUME/vol;
   nv=Ns*nb;

   b=blk_list(DFL_BLOCKS,&nb,&isw);

   if (nb==0)
   {
      alloc_bgr(DFL_BLOCKS);
      b=blk_list(DFL_BLOCKS,&nb,&isw);
   }

   error_root((nb!=(VOLUME/vol)),1,"alloc_storage [dfl_multiple.c]",
         "The number of blocks in DFL_BLOCKS does not match (VOLUME/vol)");
   error_root(((*b).vol!=vol),1,"alloc_storage [dfl_multiple.c]",
         "The size of spinors in DFL_BLOCKS does not match vol");
   error_root(((*b).ns!=(Ns+1))||((*b).nsd!=(Ns+1)),1,"alloc_storage [dfl_multiple.c]",
         "The number of spinors in DFL_BLOCKS does not match Ns");


   vs=malloc(dflmax*sizeof(*vs));
   vds=malloc(dflmax*sizeof(*vds));
   sb=malloc(dflmax*sizeof(*sb));
   sdb=malloc(dflmax*sizeof(*sdb));

   vs_size=2*Ns*nv;
   vsptr=amalloc(ndfl*vs_size*sizeof(*vsptr),ALIGN);

   vds_size=Ns*nv;
   vdsptr=amalloc(ndfl*vds_size*sizeof(*vdsptr),ALIGN);
   
   error((vsptr==NULL)||(vdsptr==NULL),1,"alloc_storage [dfl_multiple.c]",
         "Unable to allocate vector fields");
      
   sb_size=(Ns+1)*(vol+1);
   sbptr=amalloc(ndfl*nb*sb_size*sizeof(*sbptr),ALIGN);
   sdbptr=amalloc(ndfl*nb*sb_size*sizeof(*sdbptr),ALIGN);

   error((sbptr==NULL)||(sdbptr==NULL),1,"alloc_storage [dfl_multiple.c]",
         "Unable to allocate spinor fields");

   for (idfl=0;idfl<dflmax;idfl++)
   {
      vs[idfl]=NULL;
      vds[idfl]=NULL;
      sb[idfl]=NULL;
      sdb[idfl]=NULL;
      
      if(dfl_gen_parms(idfl).status==DFL_DEF)
      {
         vs[idfl]=vsptr;
         vds[idfl]=vdsptr;
         sb[idfl]=sbptr;
         sdb[idfl]=sdbptr;
         
         vsptr+=vs_size;
         vdsptr+=vds_size;
         sbptr+=nb*sb_size;
         sdbptr+=nb*sb_size;
      }
#ifdef DFL_MODES_DBG
      if (my_rank==0)
      {
         printf("dfl_multiple.c: Set up pointer for subspace %d\n",idfl);
         fflush(stdout);
      }
#endif
   }
}


void use_dfl_subspace(int idfl)
{
   int n;
   
   if(dflmax==-1)
      alloc_storage();
   
   error((idfl<0)||(idfl>=dflmax),1,"use_dfl_subspace [dfl_multiple.c]",
         "idfl is out of range");

   error(vs[idfl]==NULL,1,"use_dfl_subspace [dfl_multiple.c]",
         "The deflation subspace idfl=%d is not defined",idfl);

   if(idfl==adfl)
      return;

#ifdef DFL_MODES_DBG
   if (my_rank==0)
   {
      printf("dfl_multiple.c: Switch to dfl subspace %d from currently active dfl subspace %d\n",idfl,adfl);
      fflush(stdout);
   }
#endif
   
   if(adfl==-1)
   {
      adfl=idfl;
      set_flags(ERASED_AW);
      set_flags(ERASED_AWHAT);
      return;
   }
   
   /*
   vs_size=2*Ns*nv;
   vds_size=Ns*nv;
   */

   memcpy(vs[adfl],vflds()[0],vs_size*sizeof(**vs));
   memcpy(vds[adfl],vdflds()[0],vds_size*sizeof(**vds));

   /* sb_size=(Ns+1)*(vol+1); */

   for (n=0;n<nb;n++)
   {
      memcpy(sb[adfl]+n*sb_size,b[n].s[0],sb_size*sizeof(**sb));
      memcpy(sdb[adfl]+n*sb_size,b[n].sd[0],sb_size*sizeof(**sdb));
   }

   adfl=idfl;
   
   memcpy(vflds()[0],vs[adfl],vs_size*sizeof(**vs));
   memcpy(vdflds()[0],vds[adfl],vds_size*sizeof(**vds));

   for (n=0;n<nb;n++)
   {
      memcpy(b[n].s[0],sb[adfl]+n*sb_size,sb_size*sizeof(**sb));
      memcpy(b[n].sd[0],sdb[adfl]+n*sb_size,sb_size*sizeof(**sdb));
   }
   
   set_flags(ERASED_AW);
   set_flags(ERASED_AWHAT);
}
