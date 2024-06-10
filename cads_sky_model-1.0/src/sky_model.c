/****************************************************************************
 This program will calculate the count rate per pixel
 in a given area of the sky.
 
 Author: Jayant Murthy (jmurthy@yahoo.com)
 Program Documentation: https://docs.google.com/document/d/1BIT1I2IkJfBDlnYJI9a-E5dtHK-9ZWWv5k4p8vAcmb4/edit
 Date: Oct. 7, 2012
 Oct. 18, 2012 (JM): Initialized variables, added integrated counts, DC level
 
 **************************************************************************/
#include "sky_model.h"

int main(int argc, char *argv[], char *envp[])
{
    /*Diffuse variables*/
    struct SPECTRA bkgd_spec;
    struct  BKGD_DATA bkgd;
    float bkgd_scale, solid_angle, tot_bkgd;
    float *data;

    /*Location*/
    double ra, dec, x, y;
    
    /*Input*/
    struct  INP_PAR  inp_par;
    
    /*Output*/
    struct  WCS_HEAD wcs_in;
    long    naxes[2];
	int ix, iy;
    float *diffuse_data;
    
    /*Zodiacal variables*/
	float   *zod_data, zod_scale, tot_zod;
    struct  ZOD_DATA zodi;
    struct SPECTRA zod_spec;
    
    /*Filters*/
    struct SPECTRA filter_spec;
    
    /*Stars*/
	struct  STARS   *hipstars;
    struct  SPECTRA *stellar_spectra;
    int istars;
    float tot_stars;
    
    /*Unclassified variables*/
    int status;
	int     i, tst;
	float   *err;


/*====================== READS OUT HEADERLINES===============================*/
    
    /*Begin Initialization*/
    /*Read parameters if present*/
    tst = READ_PARAMS(&inp_par, &wcs_in);
    
	naxes[0] = wcs_in.angsize/fabs(wcs_in.xinc) + 1;
	naxes[1] = wcs_in.angsize/wcs_in.yinc + 1;
    wcs_in.xrefpix = naxes[0]/2.;
    wcs_in.yrefpix = naxes[1]/2.;
    /*We have to define the solid angle for diffuse sources*/
    solid_angle = fabs(DEGRAD(wcs_in.xinc)*DEGRAD(wcs_in.yinc));
    
	/*Allocate arrays*/
    data  = calloc(sizeof(float), naxes[0]*naxes[1]);
    err   = calloc(sizeof(float), naxes[0]*naxes[1]);
	diffuse_data = calloc(sizeof(float), naxes[0]*naxes[1]);
    
    /*Zodiacal light*/
	zod_data   = calloc(sizeof(float), naxes[0]*naxes[1]);
    zod_spec.wavelength = (float *) malloc(sizeof(float) * MAX_SPEC_ELEM);
	zod_spec.spectrum   = (float *) malloc(sizeof(float) * MAX_SPEC_ELEM);
    
    /*Diffuse background*/
    bkgd_spec.wavelength = (float *) malloc(sizeof(float) * MAX_SPEC_ELEM);
    bkgd_spec.spectrum = (float *) malloc(sizeof(float) * MAX_SPEC_ELEM);
    
    /*Filters*/
    filter_spec.wavelength = (float *) malloc(sizeof(float) * MAX_SPEC_ELEM);
    filter_spec.spectrum = (float *) malloc(sizeof(float) * MAX_SPEC_ELEM);
    
    /*Stars*/
    stellar_spectra = (struct SPECTRA *)
                    malloc(sizeof(struct STARS)*N_CASTELLI_MODELS);
    for (i = 0; i < N_CASTELLI_MODELS; ++i){
        stellar_spectra[i].spectrum = malloc(sizeof(float) * N_CASTELLI_SPECTRA);
        stellar_spectra[i].wavelength = malloc(sizeof(float) * N_CASTELLI_SPECTRA);
        stellar_spectra[i].nelements = N_CASTELLI_MODELS;
    }
    hipstars = (struct STARS *) malloc(sizeof(struct STARS)*NSTARS);
 
    /*Read the filter*/
    filter_spec.nelements = SPECT_READ(filter_spec.wavelength,
                                       filter_spec.spectrum, inp_par.filter_file);

    /*Read Galactic background if needed*/
	if (inp_par.gal_inc > 0){
        bkgd = READ_GAL_BKGD(inp_par.bkgd_model);
        bkgd_spec.nelements = SPECT_READ(bkgd_spec.wavelength,
                                         bkgd_spec.spectrum, inp_par.bkgd_spec);
        bkgd_scale = CALC_SCALE_FACTOR(filter_spec, bkgd_spec);
		CALC_BKGD_FLUX(wcs_in, naxes, data, err, inp_par.ang_limit,
                       inp_par.csc_law, bkgd);
    }
	if (inp_par.zod_inc > 0){
        zodi = ZOD_DIST_READ(inp_par.zod_file);
        zod_spec.nelements = SPECT_READ(zod_spec.wavelength,
                                        zod_spec.spectrum, inp_par.zod_spec);
		tst = CALC_ZOD_FLUX(zodi, inp_par.hour, inp_par.day, inp_par.month,
                              inp_par.year, naxes, zod_data, wcs_in);
        zod_scale = CALC_SCALE_FACTOR(filter_spec, zod_spec);
    }
	tot_bkgd  =0;
    tot_zod = 0;
	for (iy = 0; iy < naxes[1]; ++iy) {
		for (ix = 0; ix < naxes[0]; ++ix) {
            diffuse_data[ix + iy*naxes[0]] =
            data[ix + iy*naxes[0]]*bkgd_scale*solid_angle +
			zod_data[ix + iy*naxes[0]]*zod_scale*solid_angle +
            inp_par.dc_level;
            tot_bkgd += data[ix + iy*naxes[0]]*bkgd_scale*solid_angle;
            tot_zod += zod_data[ix + iy*naxes[0]]*zod_scale*solid_angle;
		}
	}

    /*Do we include stars?*/
    if (inp_par.star_inc > 0) {
        tst = READ_HIPPARCOS_CAT(inp_par, stellar_spectra, hipstars,
                                 filter_spec);
        tot_stars = 0;
        for (istars = 0; istars < NSTARS; ++istars) {
            
            status = 0;
            ra = hipstars[istars].gl;
            dec = hipstars[istars].gb;
            status = 0;
            fits_world_to_pix(ra, dec, wcs_in.xrefval, wcs_in.yrefval,
                          wcs_in.xrefpix, wcs_in.yrefpix, wcs_in.xinc,
                          wcs_in.yinc, wcs_in.rot, wcs_in.coordtype,
                          &x, &y, &status);
            ix = (int) x - 1;
            iy = (int) y - 1;
            if ((status == 0) && (ix >= 0) && (iy >= 0) && (ix < naxes[0]) && (iy < naxes[1])){
                diffuse_data[ix + iy*naxes[0]] += hipstars[istars].flux_at_earth;
                tot_stars += hipstars[istars].flux_at_earth;
            }
        }

    }/*Read all the stars into memory*/
    
    /*Write FITS file*/
    remove(inp_par.out_file);
    FITS_WRITE_FILE(inp_par.out_file, wcs_in, naxes, diffuse_data, "g");
    printf("Bkgd: %f, Zodiacal: %f, Stars: %f\n", tot_bkgd, tot_zod, tot_stars);
    
/*Free variables*/
    free(data);
    free(err);
    free(diffuse_data);
    free(zod_data);
    free(zod_spec.wavelength);
    free(zod_spec.spectrum);
    free(bkgd_spec.wavelength);
    free(bkgd_spec.spectrum);
    free(filter_spec.wavelength);
    free(filter_spec.spectrum);
    for (i = 0; i < N_CASTELLI_MODELS; ++i){
        free(stellar_spectra[i].spectrum);
        free(stellar_spectra[i].wavelength);
    }
    free(stellar_spectra);
    free(hipstars);
}   /* End of MAIN */	
