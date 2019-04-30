
/*******************************************************************************
*
* File utils.h
*
* Copyright (C) 2017 Nazario Tantalo
*
* based on the openQCD version of this file written by
* Copyright (C) 2009, 2010, 2011 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#define NAME_SIZE (10*1024)


/* UTILS_C */
extern int find_opt(int argc,char *argv[],char *opt);
extern void error(int test,int no,char *name,char *format,...);
extern void message(char *format,...);

#endif
