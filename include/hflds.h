
/*******************************************************************************
*
* File hflds.h
*
* Copyright (C) 2016 Agostino Patella
* 
* Based on openQCD-1.4/include/uflds.h
* Copyright (C) 2011, 2012, 2013 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef HFLDS_H
#define HFLDS_H

#ifndef SU3_H
#include "su3.h"
#endif

extern su3 *hfld(void);
extern su3_dble *hdfld(void);

#endif
