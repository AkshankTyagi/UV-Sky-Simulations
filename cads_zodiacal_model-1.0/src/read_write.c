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

#include "zodiacal_model.h"

void ZOD_DIST_READ(float *arr, float *hecl_lon, float *hecl_lat,
                   char *zod_file)
{
    FILE *ZODIAC=NULL;
    int ix, iy;
    float f_tmp;
    char s_tmp[FLEN_FILENAME];

    if ((ZODIAC = fopen(zod_file, "r")) == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to open Zodiacal light file %s \n",
                zod_file);
        exit(EXIT_FAILURE);
    }/*Open zodiacal light file if it exists*/

    fgets(s_tmp, FLEN_FILENAME, ZODIAC);
    fgets(s_tmp, FLEN_FILENAME, ZODIAC);

    fscanf(ZODIAC, "%f", &f_tmp);         /* Dummy read of first element */
    for (iy = 0; iy < ZOD_YSIZE; ++iy)
    {                                /* Now read the remaining latitudes */
        fscanf(ZODIAC, "%f", &f_tmp);
        *(hecl_lat + iy) = f_tmp;
    }

    for (ix = 0; ix < ZOD_XSIZE; ++ix)
    {
        fscanf(ZODIAC, "%f", &f_tmp);
        *(hecl_lon + ix) = f_tmp;
        for (iy = 0; iy < ZOD_YSIZE; ++iy)
        {
            fscanf(ZODIAC, "%f", &f_tmp);
            *(arr + ix + iy*ZOD_XSIZE) = f_tmp;
        }
    }/*Read over zodiacal file*/

    fclose(ZODIAC);
}/*End ZOD_DIST*/
/*********************************************************************/

void ZOD_SPECT_READ(float *wave, float *spec, int *ndat, char *spec_file)
{
    FILE *SPECTRUM = NULL;
    int ispec = 0;
    char s_tmp[FLEN_FILENAME], s_tmp1[FLEN_FILENAME], s_tmp2[FLEN_FILENAME];

    if ((SPECTRUM = fopen(spec_file, "r")) == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to open spectrum file %s \n",
                spec_file);
        exit(EXIT_FAILURE);
    }/*Open zodiacal spectrum if it is exists*/

    while (!feof(SPECTRUM))
    {
        fgets(s_tmp, FLEN_FILENAME, SPECTRUM);
        if((s_tmp[0] != '#') && (s_tmp[0] != '\n'))
        {
            sscanf(s_tmp, "%s %s", s_tmp1, s_tmp2);
            *(wave + ispec) = atof(s_tmp1);
            *(spec + ispec) = atof(s_tmp2);
            ispec++;
        }
    }
    fclose(SPECTRUM);
    ispec -= 1;
    *ndat = ispec;
}/*END ZOD_SPECT*/
/*********************************************************************/

void gen_default_params()
{

    FILE *fparam = NULL;
    usage();
    printf("Info   : Creating default parameter file %s\n", PARAM_FILE);

    fparam = fopen(PARAM_FILE, "w");
    if (fparam == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to create %s\n", PARAM_FILE);
        exit(EXIT_FAILURE);
    }

    fprintf(fparam,
            "#-----------------------------------------------------#\n"
            "# Parameter FILE for UVS Zodiacal light model          \n"
            "#                                                      \n"
            "# Purpose: Input parameters for running the Zodiacal   \n"
            "#          light modelling software will be read from  \n"
            "#          this file.                                  \n"
           );

    fprintf(fparam,
            "# NOTE:    Any line beginning with a hash (\"#\") or   \n"
            "#          those enclosed within \"/* -- */ \" are     \n"
            "#          comments and will be ignored by the software\n"
            "#                                                      \n"
            "# format of this file is: KEY = VALUE                  \n"
            "# value is expected to be float/string                 \n"
            "#                                                      \n"
            "#-----------------------------------------------------#\n\n"
           );

    fprintf(fparam,
            "# TIME & DATE Information:\n"
            "# -----------------------\n"
            "# 1. Time of the day (UT, hours. range: 0.0 - 24.0)\n"
            "TIME_UT      = 12.00\n\n"
            "# 2. Day of the month (unitless, range: 0 - 31 (or 28/29/30))\n"
            "DAY_OF_MONTH = 21\n\n"
            "# 3. Month of the year (unitless, range: 1 - 12) \n"
            "MONTH        = 12\n\n"
            "# 4. Year (unitless, range: 1900 - 2100 for max. accuracy) \n"
            "YEAR         = 2012\n\n\n"
           );

    fprintf(fparam,
            "# COORDINATE Information:\n"
            "# -----------------------\n"
            "# 1. RA (deg. range: 0.0 - 360.0)\n"
            "RA           = 266.50\n\n"
            "# 2. Declination (deg. range: -90 - +90)\n"
            "DECLINATION  = -29.0\n\n\n"
           );

    /*The files */
    fprintf(fparam,
            "# OUTPUT FILE:\n"
            "# -----------------------\n"
            "# The calculated spectrum is written to this file. \n"
            "OUTPUT_FILE     = zodiacal_output.txt\n\n\n"
           );

    fprintf(fparam,
            "# ZODIACAL SPECTRUM FILE:\n"
            "# -----------------------\n"
            "# zodiacal light spectrum is read from this file. \n"
            "# Data is based on the solar spectrum from Colina et al\n"
            "ZODIAC_SPECTRUM = %szodiacal_spec.txt\n\n", CADSDATA
           );

    fprintf(fparam,
            "# ZODIACAL LIGHT MODEL FILE:\n"
            "# --------------------------\n"
            "# zodiacal light distribution data (Leinert et al.) \n"
            "ZODIAC_MODEL    = %sleinert_dist.txt\n\n", CADSDATA
           );

    fprintf(fparam,
            "# --------------------------\n"
            "#\n"
           );

    printf("Info   : Default Parameter file has been generated.\n");
    fclose(fparam);

}/* END GEN_DEFAULT_PARAMS */
/*********************************************************************/

int READ_PARAMS(float *hour, int *day, int *month, int *year,
                float *ra, float *dec,
                char *zod_file, char *spec_file, char *out_file)
{

    FILE *fp;
    char line[MAX_TEXT], name[MAX_TEXT], equal[MAX_TEXT], data[MAX_TEXT];

    if ((fp = fopen(PARAM_FILE,"r")) == NULL)
    {
        gen_default_params();
        if ((fp = fopen(PARAM_FILE,"r")) == NULL) exit(EXIT_FAILURE);
    }

    /* Read input parameters */
    while (!feof(fp))
    {
        if (fgets(line, MAX_TEXT, fp) != NULL);
        if ((line[0] != '#') && (line[0] != '\n'))
        {
            sscanf(line, "%s %s %s ", name, equal, data);
            if (!strcmp("TIME_UT", name))       *hour = atof(data);
            if (!strcmp("DAY_OF_MONTH", name))   *day = atoi(data);
            if (!strcmp("MONTH", name))        *month = atoi(data);
            if (!strcmp("YEAR", name))          *year = atoi(data);
            if (!strcmp("RA", name))              *ra = atof(data);
            if (!strcmp("DECLINATION", name))    *dec = atof(data);
            if (!strcmp("OUTPUT_FILE", name))     strcpy(out_file, data);
            if (!strcmp("ZODIAC_MODEL", name))    strcpy(zod_file, data);
            if (!strcmp("ZODIAC_SPECTRUM", name)) strcpy(spec_file, data);
            strcpy(line, "End of File");
            /* This is to take care of the end of file */
        }
    }
    fclose(fp);
}/* END READ_PARAMS */
/*********************************************************************/
