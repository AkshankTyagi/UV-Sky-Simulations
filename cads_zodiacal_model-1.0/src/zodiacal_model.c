/***************************************************************************
 * This File is a part of the CADS/UV Zodiacal light model software        *
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
/*-------------------------------------------------------------------------*
 Program : cads_zodiacal_model.c
 This program will write a 2 dimensional file with the zodiacal light
 spectrum. The spectrum is based on the solar spectrum from Colina et al
 while the distribution is from Leinert et al. Documentation is in the
 cads website at http://cads.iiap.res.in/software The Sun position is
 calculated using libnova (http://libnova.sf.net)

 The output is in units of photons cm-2 s-1 sr-1 A-1

 Date  : Sep. 8, 2007
 Author: Jayant Murthy

  Modification history:

 Reks, 17th Sept 2007:
 1. Improved formatting and overall beautification :-)
 2. configure, make and added into the simulation project
 3. Removed parse_par. Inputs are now read from a parameter file,
    including the location of leinart model and zodiacal spectrum files.
 4. output written to a file instead of terminal.

JM: Sep. 20, 2007:
 1. helioecliptic longitude should be absolute value.
 2. latitude should also be fabs because we are symmetric around the equator.

Reks: 25th Sept 2007: v1.0 ready for release!
 1. Put back all the options for sending in input parameters - paramfile
    or commandline or interactive (in the order of preference).

JM: Sep. 27, 2007
 1. Changed input methods. Will read either from a parameter file or
    from a command line but not interactively.
 2. Corrected an error where I need to look only at the difference
    between solar longitude and look direction.
 3. Header of output file was wrong.

 Reks: 23 June 2012
 1. Rebranded to cads software team
 2. All coordinate transformations and calculations using libnova
 3. Removed parse_par => no command line options
 4. Can handle comments in spectrum file and removed the 'rather strict'
    ZOD_SPEC number of elements in spectrum :)
 5. More gnu-autotools friendly. Path to input data files are taken care of
*--------------------------------------------------------------------------*/
#include <libnova/solar.h>
#include <libnova/transform.h>
#include <libnova/julian_day.h>
#include <libnova/utility.h>
#include "zodiacal_model.h"

#define PROGRAM "UVS_Zodiacal_Model"

/* ------------------------------------------------------------------------- */

void PRINT_BANNER()
{
    printf("-----------------------------------------------------------\n");
    printf("%s %s\n", PROGRAM, VERSION);
    printf(" Calculates the Zodiacal light spectrum towards a given\n");
    printf(" direction, on a given date and time.\n");
    printf("-----------------------------------------------------------\n");
}/* End PRINT_BANNER */

/* ------------------------------------------------------------------------- */

void usage()
{

    printf (" Usage: uvs_zodiacal_model \n\n");
    printf ("        If no parameter file exists, new one will be created\n");
    printf ("        which may be edited as required\n\n");
    printf ("-----------------------------------------------------------\n");
    printf ("Copyright (c) 2012 The CADS Team, IIAp. [GPL v3.0 or later]\n");
    printf ("Report Bugs to: http://cads.iiap.res.in/bugzilla \n");
    printf ("Documentation : http://cads.iiap.res.in/software\n\n");
}/*End usage */

/* ------------------------------------------------------------------------- */

void get_sun_hecl(float t_ut, int day, int month, int year,
                  float ra, float dec, float *sun_hecl, float *sun_beta)
{
    double jd = 0.0;
    double f_temp = 0.0;
    struct ln_lnlat_posn sun_ecl;
    struct ln_date date;
    struct ln_equ_posn eqcoords;
    struct ln_lnlat_posn eclcoords;

    /* calculate Julian day */
    date.years   = year;
    date.months  = month;
    date.days    = day;
    date.hours   = (int) t_ut;
    date.minutes = 0;
    date.seconds = 0;
    jd = ln_get_julian_day(&date);

    /* calculate sun's ecliptic coordinates */
    ln_get_solar_ecl_coords(jd, &sun_ecl);

    eqcoords.ra  = (double) ra;
    eqcoords.dec = (double) dec;

    ln_get_ecl_from_equ(&eqcoords, jd, &eclcoords);

    /* things are symmetric around 0 latitude here..  */
    f_temp    = fabs(eclcoords.lng - sun_ecl.lng);

    /* wrap coordinates around 180 deg.*/
    if (f_temp > 180) f_temp = 360 - f_temp;

    *sun_hecl = (float) f_temp;
    *sun_beta = (float) fabs(eclcoords.lat);
}/*End calculate solar coordinates*/

/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[])
{

    char param_file[MAX_TEXT] = "", out_file[MAX_TEXT] = "";
    char zod_file[MAX_TEXT] = "", spec_file[MAX_TEXT] = "";

    float *zod_dist, *zod_spec, *zod_wave;
    float *table_ecl, *table_beta;
    float ra = 0.0, dec = 0.0, hour = 0.0;
    int   day = 0, month = 0, year = 0;
    float helio_beta = 0.0, helio_ecl = 0.0;
    int i_arg = 0;
    long hecl_index = 0, hbeta_index = 0, i = 0, wave_index = 0;
    float zod_intensity = 0.0;
    float zod_scale = 252;      /*1 unit in the table in phot units at 5000 A*/
    int ndat = 0;
    FILE *fout = NULL;

    /* --------------------------------------------------------------------- */
    /*                 Before we begin... set the signals and print a banner */

    set_signals();                            /* set signals before starting */
    PRINT_BANNER();

    /* start with some default values. Use it
       if nobody re-defined them in paramfile */
    strcpy(spec_file, DEFAULT_SPEC_FILE);
    strcpy(zod_file, DEFAULT_ZOD_FILE);
    strcpy(out_file, DEFAULT_OUT_FILE);

    /*Read parameters if present*/
    /* Get all the input parameters from parameter file */
    READ_PARAMS(&hour, &day, &month, &year, &ra, &dec, zod_file, spec_file,
                out_file);

    /* Tell the user about the input values */
    printf("Info   : Time & date       : %04.2f Hrs %02d/%02d/%04d \n",
           hour, day, month, year);
    printf("Info   : Sky coordinates   : (RA = %3.4f, Decl. = %3.4f)\n",
           ra, dec);

    /*Zodiacal light related variables*/
    zod_dist   = (float *) malloc(sizeof(float) * ZOD_XSIZE * ZOD_YSIZE);
    zod_spec   = (float *) malloc(sizeof(float) * ZOD_SPEC);
    zod_wave   = (float *) malloc(sizeof(float) * ZOD_SPEC);
    table_ecl  = (float *) malloc(sizeof(float) * ZOD_XSIZE);
    table_beta = (float *) malloc(sizeof(float) * ZOD_YSIZE);

    /*Read zodiacal light data*/
    ZOD_DIST_READ(zod_dist, table_ecl, table_beta, zod_file);
    ZOD_SPECT_READ(zod_wave, zod_spec, &ndat, spec_file);

    /* Get helio-ecliptic coordinates of Sun */
    get_sun_hecl(hour, day, month, year, ra, dec, &helio_ecl, &helio_beta);

    /* Print ecliptic coordinates*/
    printf("Info   : Helioecliptic coordinates are: %3.4f %3.4f\n",
           helio_ecl, helio_beta);

    /*Now lookup zodiacal intensity*/
    hecl_index = 0;
    while (table_ecl[hecl_index] < helio_ecl)
        ++hecl_index;
    hbeta_index = 0;
    while (table_beta[hbeta_index] < fabs(helio_beta))
        ++hbeta_index;

    zod_intensity = zod_dist[hecl_index + hbeta_index*ZOD_XSIZE];

    /*Scale zodiacal light spectrum*/
    wave_index = 0;
    while (zod_wave[wave_index] < 5000)
        ++wave_index;

    /*Scale factor for zodiacal light*/
    zod_scale *= (zod_intensity / zod_spec[wave_index]);
    for (i = 0; i < ndat; ++i)
        zod_spec[i] *= zod_scale;



    /*Print out zodiacal light*/

    fout = fopen(out_file, "w");
    if (fout == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to create output file, %s\n",
                out_file);
        exit(EXIT_FAILURE);
    }

    /* Write a header with some info on developers and all input values */
    fprintf(fout, "#------------------------------------------------#\n");
    fprintf(fout, "# Zodiacal light Spectrum generated by\n");
    fprintf(fout, "# %s %s\n", PROGRAM, VERSION);
    fprintf(fout, "# (based on Leinert's model) \n");
    fprintf(fout, "# Copyright (c) by the CADS/UV Software Team\n");
    fprintf(fout, "# URL: http://cads.iiap.res.in/software\n");
    fprintf(fout, "#------------------------------------------------#\n");
    fprintf(fout, "# user inputs\n");
    fprintf(fout, "#  Time      : %f Hrs %02d/%02d/%04d \n",
            hour, day, month, year);
    fprintf(fout, "#  Sky Coord : (RA = %3.4f, Decl. = %3.4f)\n", ra, dec);
    fprintf(fout, "#\n");
    fprintf(fout, "#  Solar spectrum : %s\n", spec_file);
    fprintf(fout, "#  Zodiacal model : %s\n", zod_file);
    fprintf(fout, "#------------------------------------------------#\n");
    fprintf(fout, "# Column 1: Wavelength (Angstroms)                \n");
    fprintf(fout, "# Column 2: Flux density (photons/cm^2/s/sr/A)    \n");
    fprintf(fout, "#------------------------------------------------#\n");

    for (i = 0; i < ndat; ++i)
        fprintf(fout, "%f %f\n", zod_wave[i], zod_spec[i]);
    close(fout);

    printf("Info   : Output written to %s\n", out_file);
    return(EXIT_SUCCESS);
}
