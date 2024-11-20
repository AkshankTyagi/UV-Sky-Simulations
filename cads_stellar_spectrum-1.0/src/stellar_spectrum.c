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

/*--------------------------------------------------------------------------
File          : mag2flux.c
Author        : Sujatha N V
Date          : 12th August, 2004
Version       : 1.0

Program Description:
This program reads sp_type, magnitude & E(B-V) from an input
parameter file & compares the spectral type in the  TABLE in file :
MKSpTypewwww.in. Selects the corresponding KURUCZ flux file
coresponding to Temp. T and gravity g of the sp. type. Then
calculate scale value for Kurucz flux and calculate the flux at
Earth for different wavelengths read from 'wave.dat'
output file contains: Wavelength, Flux at Earth F(lambda)

---------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
Rivision            : Author, dd/mm/yyyy
Revision history    :     -----

Reks, 12th July, 2005.
1. Introduced signal handling (signal_handling.c and fun_prototypes.h)
2. Modified all printf for error messages to fprintf->STDERR
   Those two modifications are a must to be able to run this on web...

Fayaz, Sept., 2005.
1. Added UVS_msgs.c, to print error messages (supersedes #2 in previous
   revision) and other fancy things like a prologue() before starting.

Reks, 08th Dec. 2005
1. Can take care of comments in input parameter file.

Reks, 07th July, 2007
1. configure-make
2. renamed project to mag2counts
3. renamed file to mag2flux
4. Better error messages and overall readability
5. moved all functions into kurucz_functions.c

JM: 10 July, 2007
1. Corrected several memory leaks in text handling

JM: 15 July, 2007
1. Cleaned up code in main program; moved parameter reading etc to functions.
   (Still have to go through other functions.)

Reks: 21 July, 2007
1. Removed SPDATA_FROOT from paramfile.
2. SPDATA_DIR now points to the top level directory that contain all data
   files, subdirectories KURUCZ_FLUX and MkSpDATA
3. shifted wave.dat and crossec0.dat into SPDATA_DIR (default location)
   experienced users may over-ride this from parameter file.
4. added command line options to create default paramfile and help message
5. Still have to fix the memory leak problems with sptype inside GET_REF1!

JM: 21 July, 2007
1. fixed all memory leaks!

Reks: 01 July 2012:
1. Rebranded to CADS/UVS, reset version to 1.0 in new cads git repository
2. Renamed to stellar_spectrum
3. read/write ops shifted to read_write.c
----------------------------------------------------------------------------*/

#include "stellar_spectrum.h"
#include <sys/types.h>

/* Main program starts here*/
/****************************************************************************/
int main (int argc, char *argv[])
{

    FILE            *out_ptr = NULL, *wave_ptr = NULL;
    struct STARS     stars;
    struct S_TABLE   mkstTG;
    int              gr = 0, tst;
    char             OutFile[MAX_TEXT];
    char             ParamFile[80];
    int c;           /* For getopt */

    /* Before we begin... set the signals and print a banner */
    set_signals();
    PRINT_BANNER();

    strcpy(ParamFile, paramfile);/*paramfile is defined in header file*/

    /* Command line options*/
    while ( (c = my_getopt (argc, argv, "gh") ) != EOF)
        switch (c)
        {
        case 'g':
            gen_stspec_params(ParamFile);
            return(EXIT_SUCCESS);
            break;
        case 'h':
            usage();
            break;
        default:
            usage();
            break;
        }

    /* Open and read the parameter file */
    READ_PARAMS(paramfile, &stars, OutFile);
    /* try opening the wave file */
    if ((wave_ptr = fopen(stars.WaveFile,"r")) == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to open file %s\n", stars.WaveFile);
        exit(EXIT_FAILURE);
    }
    /*Output file*/
    if ((out_ptr = fopen(OutFile, "w")) == NULL)
    {
        fprintf(stderr, "ERROR: Unable to write to file %s\n", OutFile);
        exit(EXIT_FAILURE);
    }

    WRITE_INFO(out_ptr, stars);/*Writes program information to file*/

    /* for all wavelengths listed in this file */
    while ((feof(wave_ptr))==0)
    {

        tst = GET_WAVE_NAME(wave_ptr, stars.dir_root, &mkstTG);
        if (tst == EXIT_SUCCESS)
        {

            /* Read the wavelength file*/
            printf("\rWorking: |");
            READ_TABLE(&mkstTG);

            /* get the crossection */
            printf("\rWorking: /");
            RCross_Sec(&stars, mkstTG.w_filename);

            /*Reads ISM Absorption CrosSection*/
            printf("\rWorking: -");
            gr  = GET_REF1(&stars, &mkstTG);

            printf("\rWorking: \\");
            if (gr == 0)
            {
                ALPHA(&stars);
                if (stars.alpha != 0.0) WRITE_LINE(out_ptr, &stars);
            }

        }/*End tst*/
    } /* end of while loop */
    printf("\rInfo   : Done\n");
    fclose(wave_ptr);
    fclose(out_ptr);

    return(EXIT_SUCCESS);
}
