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
/*
Coordinate conversion routines cannibalized from the SkyView routines of IPAC
Copyright (C) 1995, California Institute of Technology.
U.S. Government Sponsorship under NASA Contract NAS7-918 is acknowledged.
*/

#define M_PI 3.141592653589793
static double   dtr = M_PI/180;
static double   rtd = 180/M_PI;
#define debug 0 /*Do not debug*/
#define TRUE 1
#define FALSE 0
#define debug_file "COORDINATES.DEBUG"
#include "sky_model.h"

/*This program converts J2000 into galactic coordinates using a rotation matrix*/
int CONV_EQ_TO_GAL(double ra, double dec, double *gl, double *gb)
{
    int i, j;
    double x[3], xp[3], r, d, lat, lon;
    /* The following matrix comes from _The Hipparcos & Tycho */
    /* Catalogues:  Introduction & Guide to the Data_, p 92:  */
    double A[3][3] = {
        -.0548755604, -.8734370902, -.4838350155,
        .4941094279, -.4448296300,  .7469822445,
        -.8676661490, -.1980763734,  .4559837762 };
    
    r = DEGRAD(ra);
    d = DEGRAD(dec);
    x[0] = cos(r)*cos(d);
    x[1] = sin(r)*cos(d);
    x[2] = sin(d);
    
    for (i = 0; i<3; ++i) {
        xp[i] = 0;
        for (j = 0; j < 3; ++j) {
            xp[i] += A[i][j]*x[j];
        }
    }
    
    if (xp[2] > 1.)
        xp[2] = 1;
    if (xp[2] < -1.)
        xp[2] = -1;
    lat = asin(xp[2]);
    lon = atan2(xp[1], xp[0]);
    if (lon < 0)
        lon += 2.*PI;
    *gl = RADDEG(lon);
    *gb = RADDEG(lat);
    
    return(EXIT_SUCCESS);
}

/*********************************************************************/
/* 								     */
/* EC_to_EQ                                                          */
/* converts ecliptic coordinates to equatorial (both in degrees)     */
/* 								     */
/*********************************************************************/

void
EC_to_EQ(coord)
double coord[2];
{
    int             i, j;
    double          ra, dec, lat, lon, x[3], xp[3];
    double          dt1, dt2, dt0;

    static double   A[3][3] = {1.000000000, 0.000000000, 0.000000000,
			       0.000000000, 0.917436945, -0.397881203,
			       0.000000000, 0.397881203, 0.917436945};

    if (debug == TRUE)
	printf("EC_to_EQ: in %-g %-g\n", coord[0], coord[1]);

    lon = coord[0] * dtr;
    lat = coord[1] * dtr;

    x[0] = cos(lon) * cos(lat);
    x[1] = sin(lon) * cos(lat);
    x[2] = sin(lat);

    for (i = 0; i < 3; ++i)
    {
	xp[i] = 0.;
	for (j = 0; j < 3; ++j)
	{
	    xp[i] += A[i][j] * x[j];
	}
    }

    if (xp[2] > 1.)
	xp[2] = 1.;

    if (xp[2] < -1.)
	xp[2] = -1.;
    dt2 = xp[2];
    dt1 = xp[1];
    dt0 = xp[0];
    dec = asin(dt2);
    ra = atan2(dt1, dt0);

    if (ra < 0.)
	ra += 2.* M_PI;

    coord[0] = ra / dtr;
    coord[1] = dec / dtr;

    if (debug == TRUE)
	printf("EC_to_EQ: out %-g %-g\n", coord[0], coord[1]);

}



/*********************************************************************/
/* 								     */
/* GA_to_EQ                                                          */
/* converts galactic coordinates to equatorial (both in degrees)     */
/* 								     */
/*********************************************************************/

void
GA_to_EQ(coord)
double coord[2];
{
    int             i, j;
    double          ra, dec, lat, lon, x[3], xp[3];

    static double   A[3][3] = {-0.066988740, 0.492728470, -0.867600820,
			       -0.872755770, -0.450346960, -0.188374600,
			       -0.483538920, 0.744584640, 0.460199790};

    if (debug == TRUE)
	printf("GA_to_EQ: in %-g %-g\n", coord[0], coord[1]);

    lon = coord[0] * dtr;
    lat = coord[1] * dtr;

    x[0] = cos(lon) * cos(lat);
    x[1] = sin(lon) * cos(lat);
    x[2] = sin(lat);

    for (i = 0; i < 3; ++i)
    {
	xp[i] = 0.;
	for (j = 0; j < 3; ++j)
	{
	    xp[i] += A[i][j] * x[j];
	}
    }

    if (xp[2] > 1.)
	xp[2] = 1.;

    if (xp[2] < -1.)
	xp[2] = -1.;

    dec = asin(xp[2]);
    ra = atan2(xp[1], xp[0]);

    if (ra < 0.)
	ra += 2.* M_PI;

    coord[0] = ra * rtd;
    coord[1] = dec * rtd;

    if (debug == TRUE)
	printf("GA_to_EQ: out %-g %-g\n", coord[0], coord[1]);
}


/*********************************************************************/
/* 								     */
/* EQ_to_EC                                                          */
/* converts equatorial coordinates to ecliptic (both in degrees)     */
/* 								     */
/*********************************************************************/

void
EQ_to_EC(coord)
double coord[2];
{
    int             i, j;
    double          ra, dec, lat, lon, x[3], xp[3];

    static double   A[3][3] = {1.000000000, 0.000000000, 0.000000000,
			       0.000000000, 0.917436945, 0.397881203,
			       0.000000000, -0.397881203, 0.917436945};

    if (debug == TRUE)
	printf("EQ_to_EC: in %-g %-g\n", coord[0], coord[1]);

    ra = coord[0] * dtr;
    dec = coord[1] * dtr;

    x[0] = cos(ra) * cos(dec);
    x[1] = sin(ra) * cos(dec);
    x[2] = sin(dec);

    for (i = 0; i < 3; ++i)
    {
	xp[i] = 0.;
	for (j = 0; j < 3; ++j)
	{
	    xp[i] += A[i][j] * x[j];
	}
    }

    if (xp[2] > 1.)
	xp[2] = 1.;

    if (xp[2] < -1.)
	xp[2] = -1.;

    lat = asin(xp[2]);
    lon = atan2(xp[1], xp[0]);

    if (lon < 0.)
	lon += 2.* M_PI;

    coord[0] = lon * rtd;
    coord[1] = lat * rtd;

    if (debug == TRUE)
	printf("EQ_to_EC: out %-g %-g\n", coord[0], coord[1]);
}


/*********************************************************************/
/* 								     */
/* GA_to_EC                                                          */
/* converts galactic coordinates to ecliptic (both in degrees)       */
/* 								     */
/*********************************************************************/

void
GA_to_EC(coord)
double coord[2];
{
    if (debug == TRUE)
	printf("GA_to_EC: in %-g %-g\n", coord[0], coord[1]);

    GA_to_EQ(coord);
    EQ_to_EC(coord);

    if (debug == TRUE)
	printf("GA_to_EC: out %-g %-g\n", coord[0], coord[1]);
}



/*********************************************************************/
/* 								     */
/* EQ_to_GA                                                          */
/* converts equatorial coordinates to galactic (both in degrees)     */
/* 								     */
/*********************************************************************/

void
EQ_to_GA(coord)
double coord[2];
{
    int             i, j;
    double          ra, dec, lat, lon, x[3], xp[3];

    static double   A[3][3] = {-0.066988740, -0.872755770, -0.483538920,
			       0.492728470, -0.450346960, 0.744584640,
			       -0.867600820, -0.188374600, 0.460199790};

    ra = coord[0] * dtr;
    dec = coord[1] * dtr;

    x[0] = cos(ra) * cos(dec);
    x[1] = sin(ra) * cos(dec);
    x[2] = sin(dec);

    for (i = 0; i < 3; ++i)
    {
	xp[i] = 0.;
	for (j = 0; j < 3; ++j)
	{
	    xp[i] += A[i][j] * x[j];
	}
    }

    if (xp[2] > 1.)
	xp[2] = 1.;

    if (xp[2] < -1.)
	xp[2] = -1.;

    lat = asin(xp[2]);
    lon = atan2(xp[1], xp[0]);

    if (lon < 0.)
	lon += 2.* M_PI;

    coord[0] = lon * rtd;
    coord[1] = lat * rtd;

}



/*********************************************************************/
/* 								     */
/* EC_to_GA                                                          */
/* converts ecliptic coordinates to galactic (both in degrees)       */
/* 								     */
/*********************************************************************/

void
EC_to_GA(coord)
double coord[2];
{
    if (debug == TRUE)
	printf("EC_to_GA: in %-g %-g\n", coord[0], coord[1]);

    EC_to_EQ(coord);
    EQ_to_GA(coord);

    if (debug == TRUE)
	printf("EC_to_GA: in %-g %-g\n", coord[0], coord[1]);
}




/**********************************************************************/
/* 								     */
/* COORD_TO_STRING                                                    */
/* converts the coordinates to strings to be printed out.  The        */
/* variable 'sex' is used to decide whether to output in decimal      */
/* degrees or sexigesimal.                                            */
/* 								     */
/**********************************************************************/


void
coord_to_string(coord, s_coord)
double          coord[2];
char            s_coord[2][40];
{
    char            sign;
    int             hours, min, degs;
    double          dsec, c_hours, c_degs;
	int sex = TRUE;
	int s1_server_mode = FALSE;
	int output_coord_sys = 1;
	int EQ = 1;
	
    if (debug == TRUE)
	printf("coord_to_string: lon,lat %-.8g %-.8g \n", coord[0], coord[1]);

    if ((sex == TRUE) && (s1_server_mode == FALSE))  /* server gets decimal */
    {
	if (coord[0] < 0.)
	    sign = '-';
	else
	    sign = ' ';

	/* Convert longitude */
	if (output_coord_sys == EQ)
	{
	    c_hours = coord[0] / 15;	/* use hours, min, sec */
	    if (c_hours < 0)
		c_hours += 24;	/* no negative hours */
	    c_hours += .005 / 3600;   /* round up one-half of .01 sec */
	}
	else
	{
	    c_hours = fabs(coord[0]);	/* use deg, min, sec */
	    c_hours += .05 / 3600;   /* round up one-half .1 sec */
	}

	hours = c_hours;
	min = (c_hours - hours) * 60.;
	dsec = (c_hours - hours - min / 60.) * 3600.;

	if (output_coord_sys == EQ)
	{
	    /* remove the rounding (f format does it) */
	    dsec -= .005;
	    if (dsec < 0.0)
		dsec = 0.0;
	    sprintf(s_coord[0], "%dh%02dm%05.2fs", hours, min, dsec);
	}
	else
	{
	    /* remove the rounding (f format does it) */
	    dsec -= .05;
	    if (dsec < 0.0)
		dsec = 0.0;
	    sprintf(s_coord[0], "%c%dd%02dm%04.1fs", sign, hours, min, dsec);
	}



	/* Convert latitude */
	if (coord[1] < 0.)
	    sign = '-';
	else
	    sign = ' ';

	c_degs = fabs(coord[1]);/* Convert latitude     */
	c_degs += .05 / 3600;   /* round up one-half of .1 sec */

	degs = c_degs;
	min = (c_degs - degs) * 60.;
	dsec = (c_degs - degs - min / 60.) * 3600.;

	/* remove the rounding (f format does it) */
	dsec -= .05;
	if (dsec < 0.0)
	    dsec = 0.0;
	sprintf(s_coord[1], "%c%dd%02dm%04.1fs", sign, degs, min, dsec);
    }

    else			/* if sex == FALSE  */
    {
	sprintf(s_coord[0], "%f", coord[0]);
	sprintf(s_coord[1], "%f", coord[1]);
    }
}
