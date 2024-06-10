/*
 *  zodiacal_programs.c
 *  
 *  Created by Jayant Murthy on 21/04/11.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sky_model.h"
#include "fitsio.h"
#include "solar.h"

/********************** SUN_RA_DEC *****************************/
/*Using libnova*/

void SUN_RA_DEC(float t_ut, float day, float month,
                float year, double *ra, double *dec)
{
	double jd;
    struct ln_equ_posn sun_eq;
	
	/*Calculate Julian Day*/
	if (month<=2) {
		month=month+12; 
		year=year-1;
	}
	jd = (int) (365.25*year) + (int)(30.6001*(month+1)) - 15 + 1720996.5 + day + t_ut/24.0;
    ln_get_solar_equ_coords (jd, &sun_eq);
	*ra = sun_eq.ra;
	*dec= sun_eq.dec;
}/*End calculate solar coordinates*/
/*********************************************************************/

/*************** ZOD_DIST_READ ****************************/
/* Reads zodiacal light from Leinert table ***************/

struct ZOD_DATA ZOD_DIST_READ(char *zod_file)
{
	FILE *ZODIAC=NULL;
	int ix, iy;
	float f_tmp;
	char s_tmp[FLEN_FILENAME];
    struct ZOD_DATA zodi;

/*Check for existence of zodiacal light file*/	
	if ((ZODIAC = fopen(zod_file, "r")) == NULL){
		fprintf(stderr, "ERROR  : Unable to open Zodiacal light file %s \n", zod_file);
		exit(EXIT_FAILURE);
	}/*Open zodiacal light file if it exists*/

/*Zodiacal light related variables*/
	zodi.zod_dist   = (float *) malloc(sizeof(float) * ZOD_XSIZE * ZOD_YSIZE);
	zodi.table_ecl  = (float *) malloc(sizeof(float) * ZOD_XSIZE);
	zodi.table_beta = (float *) malloc(sizeof(float) * ZOD_YSIZE);

	
	/*The Zodiacal light file has two comment lines*/	
	fgets(s_tmp, FLEN_FILENAME, ZODIAC);
	fgets(s_tmp, FLEN_FILENAME, ZODIAC);
	
	fscanf(ZODIAC, "%f", &f_tmp);         /*Dummy read of first element*/
	for (iy = 0; iy < ZOD_YSIZE; ++iy){ /*Now read the remaining latitudes*/
		fscanf(ZODIAC, "%f", &f_tmp);
		zodi.table_beta[iy] = f_tmp;
	}
	
	for (ix = 0; ix < ZOD_XSIZE; ++ix){
		fscanf(ZODIAC, "%f", &f_tmp);
		zodi.table_ecl[ix] = f_tmp;/*Reading longitudes*/
		for (iy = 0; iy < ZOD_YSIZE; ++iy){
			fscanf(ZODIAC, "%f", &f_tmp);/*Read actual values*/
			zodi.zod_dist[ix + iy*ZOD_XSIZE] = f_tmp;            
		}
	}/*Read over zodiacal file*/
    
    return(zodi);
	fclose(ZODIAC);
}/*End ZOD_DIST*/
/*********************************************************************/

/*************** CALC_ZOD_SCALE_FACTOR *******************/
/* Calculates the zodiacal light over the entire sky based on
   helioecliptic coordinates.
 */

/* CALC_ZOD_SCALE_FACTOR*/

int CALC_ZOD_FLUX(struct ZOD_DATA zodi,
						  float hour, float day, float month, float year,
						  long *naxes,
						  float *zod_data, struct WCS_HEAD wcs_in)
{
	int hecl_index, hbeta_index, ix, iy;
	float zod_scale;
	float user_beta, helio_ecl;
	int status=0;
	double ra, dec, sun_ra, sun_dec, sun_ecl, sun_beta, coord[2];
		
	/*Read zodiacal light data*/	
	/*Get position of Sun*/
	SUN_RA_DEC(hour, day, month, year, &sun_ra, &sun_dec);
	
	/*We need the helioecliptic longitude and the beta angle for our model*/
	/*Convert celestial to ecliptic coordinates for Sun*/
	coord[0] = sun_ra;
	coord[1] = sun_dec;
	EQ_to_EC(coord);
	sun_ecl  = coord[0];
	sun_beta = coord[1];

/*Step through the requested field. The WCS coordinates are in the WCS structure.
  The axis size is given by naxes*/
	
	for (ix = 0; ix < naxes[0]; ++ix) {
		for (iy = 0; iy < naxes[1]; ++iy) {
			fits_pix_to_world(ix + 1, iy + 1, wcs_in.xrefval, wcs_in.yrefval,
							  wcs_in.xrefpix, wcs_in.yrefpix, wcs_in.xinc,
							  wcs_in.yinc, wcs_in.rot, wcs_in.coordtype,
							  &ra, &dec, &status);/*Returns galactic coordinates*/
			if (status == 0) {
				/*Convert galactic coordinates into helioecliptic*/
				coord[0] = ra;	coord[1] = dec;
				GA_to_EC(coord);
				helio_ecl = fabs(coord[0] - sun_ecl);
				if (helio_ecl > 180) helio_ecl = 360 - helio_ecl;
				user_beta = fabs(coord[1]);
				
				/*Now lookup zodiacal intensity*/
				hecl_index = 0;
				while (zodi.table_ecl[hecl_index] < helio_ecl)
					++hecl_index;
                if (hecl_index > 0) {
                    if ((zodi.table_ecl[hecl_index] - helio_ecl) >
                        (helio_ecl - zodi.table_ecl[hecl_index - 1]))
                        --hecl_index;
                }
				hbeta_index = 0;
				while (zodi.table_beta[hbeta_index] < user_beta)
					++hbeta_index;
                if (hbeta_index > 0) {
                    if ((zodi.table_beta[hbeta_index] - user_beta) >
                        (user_beta - zodi.table_beta[hbeta_index - 1]))
                        --hbeta_index;
                }
				zod_scale = zodi.zod_dist[hecl_index + hbeta_index*ZOD_XSIZE];
				zod_data[ix + iy*naxes[0]] = zod_scale;
			}
			status = 0;
		}
	}
	return(EXIT_SUCCESS);
}