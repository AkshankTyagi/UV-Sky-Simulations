/*
  gal_bkgd.c
  sky_model

  Created by Jayant Murthy on 07/10/12.
  Copyright (c) 2012 Jayant Murthy. All rights reserved.
*/
#include "sky_model.h"

struct BKGD_DATA READ_GAL_BKGD(char bkgd_file[FLEN_FILENAME]) {
    
    FILE *file_ptr;
    int data_index;
    double val, eval, t1, t2;
    struct BKGD_DATA gal_bkgd;
    
    file_ptr = fopen(bkgd_file, "r"); /*Open input text file*/
    if (file_ptr == NULL){
        printf("Could not find the file %s\n", bkgd_file);
        exit(EXIT_FAILURE);
    }
    
	/*Read the data. The maximum number of elements is set by MAX_BKGD_ELEM*/
	data_index = 0;
	gal_bkgd.ra   = malloc(sizeof(double) * MAX_BKGD_ELEM);
	gal_bkgd.dec  = malloc(sizeof(double) * MAX_BKGD_ELEM);
	gal_bkgd.x	 = malloc(sizeof(double) * MAX_BKGD_ELEM);
	gal_bkgd.y	 = malloc(sizeof(double) * MAX_BKGD_ELEM);
	gal_bkgd.z	 = malloc(sizeof(double) * MAX_BKGD_ELEM);
	gal_bkgd.bkgd = malloc(sizeof(double) * MAX_BKGD_ELEM);
	gal_bkgd.err  = malloc(sizeof(double) * MAX_BKGD_ELEM);
    
    /*Read UV background file*/
    while (feof(file_ptr) == 0){ /*Continue until end of input file*/
		fscanf(file_ptr, "%lf %lf %lf %lf", &t1, &t2, &val, &eval);
        /*Galactic background is in galactic coordinates*/
		gal_bkgd.ra[data_index] = t1;
		gal_bkgd.dec[data_index] = t2;
        /*Convert to cartesian coordinates*/
		gal_bkgd.x[data_index] =
        cos(DEGRAD(gal_bkgd.ra[data_index]))*
        cos(DEGRAD(gal_bkgd.dec[data_index]));
		gal_bkgd.y[data_index] =
        sin(DEGRAD(gal_bkgd.ra[data_index]))*
        cos(DEGRAD(gal_bkgd.dec[data_index]));
		gal_bkgd.z[data_index] =
        sin(DEGRAD(gal_bkgd.dec[data_index]));
        /*Background value at that position and the error*/
		gal_bkgd.bkgd[data_index] = val;
		gal_bkgd.err[data_index] = eval;
		++data_index;
	}/*End read file*/
    gal_bkgd.npoints = data_index;
    return(gal_bkgd);
    
}

int CALC_BKGD_FLUX(struct WCS_HEAD wcs_in, long naxes[2], float *data, float *var,
				   float ang_limit, float csc_law,
                   struct BKGD_DATA bkgd)
{
	
	double ra0, dec0, x0, y0, z0;
	double val, max_data=10000;
	double cang_limit, dist, cm_angle;
	int status = 0, data_index, index;
	int i, j;
    int *inc_list, inc_index = 0;
	
	
	/*ang_limit is the maximum angle for binning data but it is easier to work with
	 the cosine*/
 	cang_limit = cos(DEGRAD(ang_limit));
    data_index = bkgd.npoints;
    
    /*I don't have to search the huge dataset each time 
      so I'll restrict the search range.*/
    inc_list = (int *) malloc(sizeof(int)*data_index);
    /*Go overboard on the limiting size*/
    cm_angle = cos(DEGRAD(wcs_in.angsize*1.5 + ang_limit));

    for (index = 0; index < data_index; ++index) {
        ra0 = wcs_in.xrefval;
        dec0 = wcs_in.yrefval;
        x0 = cos(DEGRAD(ra0))*cos(DEGRAD(dec0));
        y0 = sin(DEGRAD(ra0))*cos(DEGRAD(dec0));
        z0 = sin(DEGRAD(dec0));
        dist = bkgd.x[index]*x0 + bkgd.y[index]*y0 + bkgd.z[index]*z0;
        if (dist > cm_angle) {
            inc_list[inc_index] = index;
            ++inc_index;
        }
    }
  	for (i = 0; i < naxes[0]; ++i){
		for (j = 0; j < naxes[1]; ++j) {
			status = 0;
            /*Get the l and b for each pixel in the FITS file*/
			fits_pix_to_world(i + 1, j + 1, wcs_in.xrefval, wcs_in.yrefval,
							  wcs_in.xrefpix, wcs_in.yrefpix, wcs_in.xinc,
							  wcs_in.yinc, wcs_in.rot, wcs_in.coordtype,
							  &ra0, &dec0, &status);
			
			if (status == 0){/*If the point is valid*/
                /*Convert coordinates into cartesian coordinates*/
				x0 = cos(DEGRAD(ra0))*cos(DEGRAD(dec0));
				y0 = sin(DEGRAD(ra0))*cos(DEGRAD(dec0));
				z0 = sin(DEGRAD(dec0));
                /*Run through all the data points find the angular distance between them*/
				for (index = 0; index < inc_index; ++index) {
					dist = bkgd.x[inc_list[index]]*x0 +
                           bkgd.y[inc_list[index]]*y0 +
                           bkgd.z[inc_list[index]]*z0;
                    /*If the distance is less than our limit, we bin weighting with the square
                     of the cos(angle) Note that because we want nearer points to be weighted
                     higher we actually use 1/cos(angle)*/
					if (dist > cang_limit){
						data[i + j*naxes[0]] += bkgd.bkgd[inc_list[index]]*dist*dist;
						var[i + j*naxes[0]] += dist*dist;
					}/*Endif*/
				}/*End check through all data values*/
				if (var[i + j*naxes[0]] > 0){
                    /*Divide by the weights*/
					data[i + j*naxes[0]] /= var[i + j*naxes[0]];
					var[i + j*naxes[0]] = sqrt(var[i + j*naxes[0]]);
				} else {/*If there is no valid data*/
					data[i + j*naxes[0]] = -1;
					var[i + j*naxes[0]] = 1e9;
				}
				if (data[i + j*naxes[0]] > max_data)
                    max_data = data[i + j*naxes[0]];
			}/*End add data*/
		}/*End j loop*/
	}/*End i loop*/
	
    /*Fill in the empty areas using a cosecant law*/
	for (i = 0; i < naxes[0]; ++i){
		for (j = 0; j < naxes[1]; ++j) {
			if (var[i + j*naxes[0]] >= 1e9) {
				status = 0;
				fits_pix_to_world(i + 1, j + 1, wcs_in.xrefval, wcs_in.yrefval,
								  wcs_in.xrefpix, wcs_in.yrefpix, wcs_in.xinc,
								  wcs_in.yinc, wcs_in.rot, wcs_in.coordtype,
								  &ra0, &dec0, &status);
				if (fabs(dec0) > 0) val = csc_law/sin(fabs(dec0)*RADEG);
                if (val > max_data)
                    val = max_data;/*Cap the values*/
				data[i + j*naxes[0]] = max_data;
			}
		}
	}
	
	return(EXIT_SUCCESS);
	
}
