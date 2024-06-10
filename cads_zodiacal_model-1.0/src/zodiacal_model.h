/***************************************************************************
 * This File is a part of the CADS Zodiacal light model software           *
 *   Copyright (C) 2012         by The CADS Software Team,                 *
 *                              Indian Institute of Astrophysics,          *
 *                              Bangalore 560 034                          *
 *                              cads_AT_iiap.res.in                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/*--------------------------------------------------------------------------
  zodiacal_model.h: Contains various definitions for zodiacal_model.c

  Revision history:       ----
  Reks, 25th September, 2007
  1. Created this file
  2. Shifted all #define statements from main to here
  Reks, 22 June 2012
  1. Rebranded to cads software group
 *--------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Maximum number of elements in a text line                     */
#define MAX_TEXT 150

/* Define the parameters for reading the zodiacal file           */
#define ZOD_XSIZE 19                            /*Number of rows */
#define ZOD_YSIZE 10                        /* Number of Columns */
#define ZOD_SPEC 25763                     /* Number of elements */
#define RAD_PER_DEG 0.0174532925           /* radians per degree */
#define FLEN_FILENAME 1025          /* max length of a filename  */
#define PARAM_FILE "zodiacal_initparams.txt"     /* input params */
#define DEFAULT_ZOD_FILE "leinert_dist.txt"   /* Leinert's model */
#define DEFAULT_SPEC_FILE "zodiacal_spec.txt"  /* Solar Spectrum */
#define DEFAULT_OUT_FILE "zodiacal_output.txt"    /* Output file */

/* path to input data files, if we were GNU-built by autotools   */
#ifndef CADSDATA
#define CADSDATA ' '
#endif

/* this is something else */
#define INITVAL -400
