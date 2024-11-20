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
#include "stellar_spectrum.h"

void gen_stspec_params(char *ParamFile)
{
    FILE *fparam = NULL;
    printf("Info   : Creating default parameter file %s\n", ParamFile);

    fparam = fopen(ParamFile,"w");
    if (fparam == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to create %s\n", ParamFile);
        exit(EXIT_FAILURE);
    }

    fprintf(fparam,
            "#-----------------------------------------------------#\n"
            "# CADS/UVS Stellar Spectrum Calculator - Parameter FILE\n"
            "# \n"
            "# Purpose: Input parameters for generating a stellar \n"
            "#          spectrum using KURUCZ models. \n"
            "# \n"
            "# NOTE:    Any line beginning with a hash (\"#\") or\n"
            "#          those enclosed within \"/* -- */ \" are \n"
            "#          comments and will be ignored by the software. \n"
            "# \n"
            "# format of this file is: KEY = VALUE \n"
            "#-----------------------------------------------------#\n\n\n"
           );

    fprintf(fparam,
            "# Source Information:  \n"
            "# -------------------  \n"
            "# 1. Source V magnitude\n"
            "V_MAGNITUDE   = 0.0    \n\n"
            "# 2. Spectral type     \n"
            "SPECTRAL_TYPE = G2V    \n\n"
            "# 3. E(B - V)          \n"
            "E(B-V)        = 0.0    \n\n\n"
           );

    fprintf(fparam,
            "# Output file \n"
            "# -------------------       \n"
            "OUTPUT_FILE   = spectrum.dat\n\n\n"
           );

    fprintf(fparam,
            "# Stellar Data files information:       \n"
            "# ----------------------------------    \n"
            "# 1. Directory containing Spectral data \n"
            "SPDATA_DIR        = %sspectral_data    \n\n", CADSDATA
           );

    fprintf(fparam,
            "#************************************** \n"
            "# NO MORE CHANGES ARE REQUIRED IN THIS  \n"
            "# SECTION UNDER NORMAL CIRCUMSTANCES    \n"
            "#************************************** \n"
            "# 2. The Cross Section File. Default is \n"
            "#    <SPDATA_DIR>/crossec0.dat          \n"
            "# CROSS_SECTION_FILE = crossec0.dat     \n\n"
            "# 3. File containing list of wavelengths\n"
            "#     for which flux will be calculated \n"
            "#     Default is <SPDATA_DIR>/wave.dat  \n"
            "# WAVELENGTH_FILE    = wave.dat         \n\n"
           );

    printf("Info   : Modify the file as per your taste and re-run\n\n");
    fclose(fparam);
}

/****************************************************************************/
void PRINT_BANNER()
{
    printf("-----------------------------------------------------------\n");
    printf(" %s %s  \n", PROGRAM, VERSION);
    printf(" Calculates stellar flux using KURUCZ models\n");
    printf("-----------------------------------------------------------\n");
}

void usage()
{
    printf (" Usage: uvs_stellar_spectrum [options]\n\n");
    printf ("    options are:\n");
    printf ("    [-g] Generates the default parameter file \n");
    printf ("    [-h] help \n\n" );
    printf (" In the absence of command line options, the program would\n");
    printf (" read inputs on stellar  magnitude, reddening and spectral \n");
    printf (" type from  the parameter file and generates  the  stellar \n");
    printf (" spectrum using  kurucz stellar models.                    \n");
    printf ("-----------------------------------------------------------\n");
    printf ("Copyright (c) 2012 CADS Team, IIAp. [GPL v3.0 or later]\n");
    printf ("Report Bugs to: cads.iiap.res.in/bugzilla\n" );
    printf ("Documentation : cads.iiap.res.in/software\n\n");
    exit(1);
}
/***************************************************************/

void READ_PARAMS(char *ParamFile, struct STARS *stars, char *OutFile)
{
    FILE *file_ptr = NULL;
    char   line[MAX_TEXT], name[MAX_TEXT], equal[MAX_TEXT], data[MAX_TEXT];
    char   tmpstr[MAX_TEXT];  /* for temperory usage */

    if ((file_ptr = fopen(ParamFile,"r")) == NULL)
    {
        fprintf(stderr,"ERROR  : Unable to open the input parameter file\n");
        gen_stspec_params(ParamFile);
        exit(EXIT_FAILURE);
    }

    /* read the parameter file */
    while ((fgets(line, MAX_TEXT, file_ptr)) != NULL)
    {
        if ((line[0] != '#') && (line[0] != '\n'))
        {
            sscanf(line, "%s %s %s ", name, equal, data);
            if (!strcmp("V_MAGNITUDE", name)) stars->mv = atof(data);
            if (!strcmp("E(B-V)", name)) stars->Ebv = atof(data);
            if (!strcmp("SPECTRAL_TYPE", name)) strcpy(stars->sp_type, data);
            if (!strcmp("OUTPUT_FILE", name)) strcpy(OutFile, data);
            if (!strcmp("SPDATA_DIR", name)) strcpy(stars->dir_root, data);

            /* Now prepare default values for Wavefile and Sigmafile..
                                   (will use these if they are not defined in paramfile ) */
            strcpy(tmpstr, stars->dir_root);
            strcat(tmpstr, "/");
            strcpy(stars->WaveFile, tmpstr);
            strcpy(stars->SigmaFile, tmpstr);

            if (!strcmp("CROSS_SECTION_FILE", name)) strcpy(stars->SigmaFile, data);
            if (!strcmp("WAVELENGTH_FILE", name)) strcpy(stars->WaveFile, data);
            strcpy(line, "End of File");
        }
    }/*End File Read*/
    fclose(file_ptr);/* close the parameter file */

    /* Tell the user about the input values */
    printf("Info   : V Magnitude        = %f\n", stars->mv);
    printf("Info   : Spectral Type      = %s\n", stars->sp_type);
    printf("Info   : Reddening          = %f\n", stars->Ebv);
    printf("Info   : Spectral Data dir  = %s\n", stars->dir_root);

    /* Print out cross section filename only if user specified the location.
              Otherwise, load default value and keep quiet */
    if (strcmp(stars->SigmaFile, tmpstr) != 0)
        printf("Info   : Cross section file = %s\n", stars->SigmaFile);
    else strcat(stars->SigmaFile, SIGMA_FILE);


    /* Same for the Wavelength filename */
    if (strcmp(stars->WaveFile, tmpstr) != 0)
        printf("Info   : Wavelength file    = %s\n", stars->WaveFile);
    else strcat(stars->WaveFile, WAVELENGTH_FILE);

    printf("Info   : Output file        = %s\n", OutFile);

}/*End READ_PARAM*/
/****************************************************************************/
int GET_WAVE_NAME(FILE *wave_ptr, char *dir_root, struct S_TABLE *mkstTG)
{
    float W0;
    char wwww[6];
    char filename[MAX_TEXT];
    int tst;

    tst = fscanf(wave_ptr, "%f", &W0);
    if (tst < 0) return(EXIT_FAILURE);
    mkstTG->OWlength=W0;

    /* make the file name as MKSpTypewwww.in */
    /* first copy the dir/root */
    strcpy(filename, dir_root);
    strcat(filename, "/");
    strcat(filename, MKSPDATA_DIR);
    strcat(filename, "/");
    strcat(filename, MKSPDATA_FROOT);
    sprintf(wwww,"%.0f",mkstTG->OWlength);
    strcat(filename,wwww);
    strcat(filename,".in");

    /* copy the file name to structure field */
    strcpy(mkstTG->w_filename, filename);
    return(EXIT_SUCCESS);
}
/***************************************************************************/

void WRITE_INFO(FILE *out_ptr, struct STARS stars)
{

    /* Write a header with some info on developers and all input values */
    fprintf(out_ptr, "#------------------------------------------------#\n");
    fprintf(out_ptr, "# Stellar Spectrum generated by %s %s\n", PROGRAM, VERSION);
    fprintf(out_ptr, "# (based on KURUCZ model) \n");
    fprintf(out_ptr, "# Copyright (c) by the CADS Team, IIAp\n");
    fprintf(out_ptr, "# Bug reports: http://cads.iiap.res.in/bugzilla/\n");
    fprintf(out_ptr, "# Project Website: http://cads.iiap.res.in\n");
    fprintf(out_ptr, "#------------------------------------------------#\n");
    fprintf(out_ptr, "# user inputs:\n");
    fprintf(out_ptr, "#             V magnitude   = %f\n", stars.mv);
    fprintf(out_ptr, "#             Spectral Type = %s\n", stars.sp_type);
    fprintf(out_ptr, "#             E(B - V)      = %f\n", stars.Ebv);
    fprintf(out_ptr, "#------------------------------------------------#\n");
    fprintf(out_ptr, "# Column 1: Wavelength (Angstroms)                \n");
    fprintf(out_ptr, "# Column 2: Flux density (ergs/cm^2/s/A)          \n");
    fprintf(out_ptr, "#------------------------------------------------#\n");
}

