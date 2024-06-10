/*
  sky_model.h
  sky_model

  Created by Jayant Murthy on 01/10/12.
  Copyright (c) 2012 Jayant Murthy. All rights reserved.
*/

#ifndef sky_model_sky_model_h
#define sky_model_sky_model_h


/***************************************************************************
 * This File is a part of the UVS Zodiacal light model software            *
 *   Copyright (C) 2007         by The Tauvex Software Team,               *
 *                              Indian Institute of Astrophysics,          *
 *                              Bangalore 560 034                          *
 *                              tauvex_AT_iiap.res.in                      *
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
 uvs_zodiacal_model.h: Contains various definitions for
 uvs_zodiacal_model.c
 URL: tauvex.iiap.res.in/tauwiki/

 Revision history:       ----
 Reks, 25th September, 2007
 1. Created this file
 2. Shifted all #define statements from main to here
 *---------------------------------------------------------------------------*/

/*Include files*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include "parse_par.h"


/*Define program specific parameters*/
#define BIGERR 1e32

/*Physical Constants*/
#define PI            3.141592653589793238462643
#define RADEG         0.017453292	/* Conversion from degrees to radians */
#define DEGRAD(x)    ((x)*PI/180.)
#define RADDEG(x)    ((x)*180./PI)
#define PCinAU       206265
#define PCinCM(x)    (3.085678E+18 * (x)) /* Converts from Parsec to cm */
#define eRStar       0.0       /* Average Stellar radius (in cm ) */
#define rd           0        /* distance to the reference plane from earth(in PC)*/
#define ERG_TO_PHOT  50306871.92

/*Program specific variables*/
#define HIP_DELIM           "|"
#define HIP_MAIN_REC_LEN    1000
#define MAX_FILE_LENGTH     80
#define N_CASTELLI_SPECTRA  1121
#define N_CASTELLI_MODELS   76
#define ZOD_XSIZE           19              /*Number of rows from Leinert */
#define ZOD_YSIZE           11              /* Number of Columns */
#define MAX_BKGD_ELEM       100000          /*Number of elements in bkgd file*/
#define NSTARS              120000          /*Number of stars in Catalog*/
#define MAX_SPEC_ELEM       100000          /*Maximum number of spectral elements*/

/*Default file names*/
#define PARAM_FILE "diffuse_initparams.txt"         /* input params */
#define DEFAULT_ZOD_FILE "leinert_dist.txt"         /* Leinert's model */
#define DEFAULT_ZOD_SPEC_FILE "zodiacal_spec.txt"   /* Solar Spectrum */
#define DEFAULT_OUT_FILE "diffuse_output.txt"       /* Output file */

/*Structures*/
struct WCS_HEAD {
	double xrefval;
	double yrefval;
	double xrefpix;
	double yrefpix;
	double xinc;
	double yinc;
	double rot;
    double angsize;
	char coordtype[5];
};

struct BKGD_DATA{
    double *ra;
    double *dec;
    double *x;
    double *y;
    double *z;
    double *bkgd;
    double *err;
    int npoints;
};

struct ZOD_DATA{
    float *zod_dist;
    float *table_ecl; /*heliocliptic longitude*/
    float *table_beta;/*ecliptic latitude*/
};
struct INP_PAR {
    char par_file[MAX_FILE_LENGTH];
    char out_file[MAX_FILE_LENGTH];
    char wave_file[MAX_FILE_LENGTH];
    char sigma_file[MAX_FILE_LENGTH];
    char location_file[MAX_FILE_LENGTH];
    char hipparcos_file[MAX_FILE_LENGTH];
    char HIPFile[MAX_FILE_LENGTH];
    char MKSpFile[MAX_FILE_LENGTH];
    char castelli_file[MAX_FILE_LENGTH];
    char out_dir[MAX_FILE_LENGTH];
    char zod_file[MAX_FILE_LENGTH];
    char zod_spec[MAX_FILE_LENGTH];
    char bkgd_model[MAX_FILE_LENGTH];
    char bkgd_spec[MAX_FILE_LENGTH];
    char filter_file[MAX_FILE_LENGTH];
    float hour;
    float day;
    float month;
    float year;
    float ang_limit;
    float csc_law;
    float dc_level;
    int zod_inc;
    int gal_inc;
    int star_inc;

};

/*Definition of the data structure for Hipparcos stars*/
struct STARS
{
    long   HIP_NO;             /* Hipparcos Catalog No             */
    long   HD_NO;              /* HD Number                        */
    char   sp_type[12];        /* Spectral Type.                   */
    float  distance;           /* Distance to the Star (in pc)     */
    float  ra;                 /* Right Ascension (2000)           */
    float  dec;                /* Declination                      */
    double x;                  /* X Cartesian Coordinate.          */
    double y;                  /* Y Cartesian Coordinate.          */
    double z;                  /* Z Cartesian Coordinate.          */
    float  B_mag;              /* B mag of the star                */
    float  V_mag;              /* V mag of the star                */
    float  gl;                 /* Galactic longitude l             */
    float  gb;                 /* Galactic latitude  b             */
    int  temperature;          /* index of temperature of the star */
    int nstars;                /*Number of elements*/
    float scale;               /*Scale factor for stellar bightness*/
    float ebv;                 /*E(B-V)*/
    float flux_at_earth;       /*Counts from star at the Earth*/
};

struct SPECTRA
{
    char filename[MAX_FILE_LENGTH];
    float *spectrum;
    float *wavelength;
    float nelements;
};


/*Function Definitions*/
void EC_to_EQ(double coord[2]);
void EQ_to_EC(double coord[2]);
void GA_to_EQ(double coord[2]);
void GA_to_EC(double coord[2]);
void EC_to_GA(double coord[2]);
void EQ_to_GA(double coord[2]);
void coord_to_string(double coord[2], char s_coord[2][40]);
int CALC_BKGD_FLUX(struct WCS_HEAD wcs_in, long naxes[2], float *data,
                   float *var, float ang_limit, float csc_law,
                   struct BKGD_DATA bkgd);
void SUN_RA_DEC(float t_ut, float day, float month,
                float year, double *ra, double *dec);
int CALC_ZOD_FLUX(struct ZOD_DATA zodi,
						  float hour, float day, float month, float year,
						  long *naxes,
						  float *zod_data, struct WCS_HEAD wcs_in);
int SPECT_READ(float *wave, float *spec, char *spec_file);
struct ZOD_DATA ZOD_DIST_READ(char *zod_file);
int READ_PARAMS(struct INP_PAR *inp_par, struct WCS_HEAD *wcs_in);
void FITS_WRITE_FILE(char file_out[FLEN_FILENAME], struct WCS_HEAD wcs_out,
					 long naxes[2], float *data, char gal[3]);
void GEN_DEFAULT_PARAMS();
struct BKGD_DATA READ_GAL_BKGD(char bkgd_file[FLEN_FILENAME]);
int CONV_EQ_TO_GAL(double ra, double dec, double *gl, double *gb);
int READ_CROSS_SEC(struct SPECTRA line, struct INP_PAR inp_par);
struct INP_PAR READ_INPUT_PARAMETERS(int argc, char *argv[]);
int HIP_READ_LINE(FILE *star_file, struct STARS *line);
int READ_CASTELLI_SPECTRA(char *, struct SPECTRA *);
int GET_STAR_TEMP(struct STARS *);
int GET_SCALE_FACTOR(struct STARS *hipstars,
                     struct SPECTRA *stellar_spectra);
float CALC_FLUX(struct STARS *hipstars, struct INP_PAR files,
                struct SPECTRA cross_sec,
                struct SPECTRA *stellar_spectra,
                struct SPECTRA filters,
                float gl, float gb, float dist, float lambda,
                float ave_sigma);
float CALC_SCALE_FACTOR(struct SPECTRA x, struct SPECTRA y);
int READ_HIPPARCOS_CAT(struct INP_PAR inp_par, struct SPECTRA *stellar_spectra,
                       struct STARS *hipstars, struct SPECTRA filter_wave);




#endif

#ifndef CADSDATA
#define CADSDATA "./"
#endif