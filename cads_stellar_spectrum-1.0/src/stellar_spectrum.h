/***************************************************************************
 * This File is a part of the CADS/UV Stellar Spectrum model software      *
 *   Copyright (C) 2012 by CADS/UV Software Team,                          *
 *                         Indian Institute of Astrophysics                *
 *                         Bangalore 560034                                *
 *                         cads_AT_iiap.res.in                             *
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
/********************************************************************
*Modelling of the UV Sky
*Sujatha N V & Jayant Murthy
********************************************************************
*file :        sptmagflux.h
*
*Description:     Header file for the program(s)
*                 "SpTMAGFlux.c"
*
*Last modified date: Aug. 04, 2004

Reks, 07th July, 2007
 Renamed file to mag2flux.h, as a part of re-organizing project
to mag2counts.

Reks, 01 July 2012
 Another rename to stellar_spectrum.h as a part of rebranding the
 whole thing to cads.
*********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef CADSDATA
#define CADSDATA ' '
#endif

#define PROGRAM "UVS_stellar_spectrum"

#define N_Tmax 4000
#define N_FLUX 4000
#define TAB_LEN 150
#define TAB_DELIM "|"
/*Maximum number of elements in a data file*/
#define MAX_ELEMENTS 10000
/*Maximum number of elements in a text line*/
#define MAX_TEXT 150

/* The input parameter file */
#define paramfile "stellarspec_initparams.txt"

/* directory containing spectral data files
   (inside SPDATA_DIR, read from parameter file) */
#define MKSPDATA_DIR "MKSpDATA"

/* Root name for the spectral data files */
#define MKSPDATA_FROOT "MKSpType"

/* File containing wavelengths for which flux will be calculated */
#define WAVELENGTH_FILE "wave.dat"

/* File with cross sections at wavelengths specified in WAVELENGTH_FILE */
#define SIGMA_FILE "crossec0.dat"

/*Definition of the data structure of the Table Spec-Type, Temp & g */
struct S_TABLE
{
    char           star_name[N_Tmax][13];     /* Spectral Type.                     */
    char           TG_name[N_Tmax][MAX_TEXT]; /* Kurucz flux File name              */
    float          Bflux[N_Tmax];             /* Kurucz Flux in B filter            */
    float          Vflux[N_Tmax];             /* Kurucz Flux for V filter (5500Å).  */
    float          Oflux[N_Tmax];             /* Kurucz Flux under investigation.   */
    float          OWlength;                  /* Wavelength of investigation.       */
    int            N_tab;                     /* Number of TABLE items.             */
    char           KFluxdir0[MAX_TEXT];       /* Directory of the Kurucz flux Files */
    char           w_filename[MAX_TEXT];      /* String to hold MKSpTypewwww.in     */
};

/*Definition of the data structure for Hipparcos stars*/
struct STARS
{
    char   sp_type[22];        /* Spectral Type.                   */
    char   tg_type[16];        /* Kurucz flux  FILE NAME           */
    float  Wlength[N_FLUX];    /* Wavelength of investigation.     */
    double C_Section[N_FLUX];  /* Extinction Cross-Section         */
    float  mv;                 /* Visual Magnitude of the Star.    */
    double Bflux;              /* B - V expected                   */
    double Vflux;              /* V flux form Kurucz flux          */
    double Oflux;              /* O flux form Kurucz flux          */
    float  alpha;              /* Scale value for Kurucz flux      */
    float  Ebv;                /* Jhonson V mag                    */
    char SigmaFile[MAX_TEXT];  /* Crosssection file                */
    char WaveFile[MAX_TEXT];   /* Wavelength File                  */
    char dir_root[MAX_TEXT];   /*Location of the spectral data     */
};
