//
//  coordinates.c
//  ISRF
//
//  Created by Jayant Murthy on 11/09/12.
//  Copyright (c) 2012 Jayant Murthy. All rights reserved.
//

#include "diffuse_model.h"


#define PI 3.141592653589793238462643
#define RADEG         0.017453292	/* Conversion from degrees to radians */
#define DEGRAD(x)    ((x)*PI/180.)
#define RADDEG(x)    ((x)*180./PI)

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

/*****************   Begin Cart_to_celest  *******************************/
/*Convert from x, y, z to ra and dec*/
void CART_TO_CELEST(double xp, double yp, double zp, double *ra, double *dec, struct INP_PAR inp_par)
{
    double dist, dxy, cra, sra, sdec;
    double  x, y, z;
    
    x = (xp -inp_par.sun_x)*inp_par.dust_binsize;
    y = (yp -inp_par.sun_y)*inp_par.dust_binsize;
    z = (zp -inp_par.sun_z)*inp_par.dust_binsize;
    
    dist = sqrt(x*x + y*y + z*z);
    dxy  = sqrt(x*x + y*y);
    sdec = z/dist;
    cra  = x/dxy;
    sra  = y/dxy;
    
    if (sra >= 0)
        *ra = RADDEG(acos(cra));
    else *ra = 360 - RADDEG(acos(cra));
    *dec = RADDEG(asin(sdec));
 }

/********************  Begin CALC_DIST ****************************/
double CALC_DIST(double a, double b, double c, double x, double y, double z)
{
    double dist;
    dist = pow(a - x,2) + pow(b - y,2) + pow(c - z,2);
    return (dist);
}

/**********************WCS_READ********************************************/
struct WCS_HEAD WCS_READ(char wcs_file[MAX_FILE_LENGTH])
{
    FILE *fptr;
    struct WCS_HEAD wcs;
    double dtmp;
    char    coordtype[5];
    long    itmp;
    
    fptr = fopen(wcs_file, "r");
    fscanf(fptr, "%lf", &dtmp);
    wcs.xrefval = dtmp;
    fscanf(fptr, "%lf", &dtmp);
    wcs.yrefval = dtmp;
    fscanf(fptr, "%lf", &dtmp);
    wcs.xrefpix = dtmp;
    fscanf(fptr, "%lf", &dtmp);
    wcs.yrefpix = dtmp;
    fscanf(fptr, "%lf", &dtmp);
    wcs.xinc = dtmp;
    fscanf(fptr, "%lf", &dtmp);
    wcs.yinc = dtmp;
    fscanf(fptr, "%lf", &dtmp);
    wcs.rot = dtmp;
    fscanf(fptr, "%s", coordtype);
    strcpy(wcs.coordtype, coordtype);
    fscanf(fptr, "%li", &itmp);
    wcs.nx = itmp;
    fscanf(fptr, "%li", &itmp);
    wcs.ny = itmp;
    return (wcs);
}
/******************************************************************/