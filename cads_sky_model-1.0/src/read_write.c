/*
 *  read_write.c
 *  diffuse
 *
 *  Created by Jayant Murthy on 21/04/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include "sky_model.h"


/*Function to create FITS header*/
void FITS_WRITE_FILE(char file_out[FLEN_FILENAME], struct WCS_HEAD wcs_out,
					long naxes[2], float *data, char gal[3])
{
    int status=0;
    char str_out[FLEN_KEYWORD];
    fitsfile *fptr;

    if (fits_create_file(&fptr, file_out, &status)){/*Create FITS file*/
		printf("Error in creating %s", file_out);
		exit(EXIT_FAILURE);
    }
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &wcs_out.xrefval,
				   "REF POINT VALUE IN DEGREES", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &wcs_out.xrefpix,
				   "REF POINT PIXEL LOCATION", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", &wcs_out.xinc,
				   "DEGREES PER PIXEL", &status);
    fits_write_key(fptr, TDOUBLE, "CROTA1", &wcs_out.rot,
				   "ROTATION FROM ACTUAL AXIS", &status);
    if (strcmp(gal, "g") == 0) strcpy(str_out,"GLON");
    else strcpy(str_out,"RA--");
    strcat(str_out, wcs_out.coordtype);
    fits_write_key(fptr, TSTRING, "CTYPE1", &str_out,
                   "COORDINATE TYPE", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &wcs_out.yrefval,
                   "REF POINT VALUE IN DEGREES", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &wcs_out.yrefpix,
                   "REF POINT PIXEL LOCATION", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", &wcs_out.yinc,
                   "DEGREES PER PIXEL", &status);
    fits_write_key(fptr, TDOUBLE, "CROTA2", &wcs_out.rot,
                   "ROTATION FROM ACTUAL AXIS", &status);
    if (strcmp(gal, "g") == 0) strcpy(str_out,"GLAT");
    else strcpy(str_out,"DEC-");
    strcat(str_out, wcs_out.coordtype);
    fits_write_key(fptr, TSTRING, "CTYPE2", &str_out,
                   "COORDINATE TYPE", &status);
    fits_write_2d_flt(fptr, 1, naxes[0], naxes[0], naxes[1], data, &status);
	fits_close_file(fptr, &status);
}/*End of FITS_WRITE_FILE*/
/*********************************************************************/

void GEN_DEFAULT_PARAMS() {

	FILE *fparam = NULL;

	fparam = fopen(PARAM_FILE, "w");
	if(fparam == NULL){
		fprintf(stderr, "ERROR  : Unable to create %s\n", PARAM_FILE);
		exit(EXIT_FAILURE);
	}

	fprintf(fparam,
			"#-----------------------------------------------------#\n"
			"# Parameter FILE for sky model                         \n"
			"#                                                      \n"
			"# Purpose: Input parameters for running the Diffuse    \n"
			"#          light modelling software will be read from  \n"
			"#          this file.                                  \n"
			);

	fprintf(fparam,
			"# NOTE:    Any line beginning with a hash (\"#\") or   \n"
			"#          those enclosed within \"/* -- */ \" are     \n"
			"#          comments and will be ignored by the software\n"
			"#                                                      \n"
			"# format of this file is: KEY = VALUE                  \n"
			"#                                                      \n"
			"#-----------------------------------------------------#\n\n"
			);

	fprintf(fparam,
			"# TIME & DATE Information:\n"
			"# -----------------------\n"
			"# 1. Time of the day (UT, hours. range: 0.0 - 24.0)\n"
			"# 2. Day of the month (unitless, range: 0 - 31)\n"
			"# 3. Month of the year (unitless, range: 1 - 12) \n"
			"# 4. Year (unitless, range: 1900 - 2100 for max. accuracy) \n\n"
			"TIME_UT      = 12.00\n"
			"DAY_OF_MONTH = 1.0\n"
			"MONTH        = 8\n"
			"YEAR         = 2007\n\n"
			);

	/*The files */
        fprintf(fparam,
                "# List of files:\n"
                "# -----------------------\n"
               );
        fprintf(fparam,
                "OUTPUT_FILE      = diffuse_output.fits\n"
               );
        fprintf(fparam,
                "ZODIAC_MODEL     = %szodiacal_data/leinert_dist.txt\n",
                CADSDATA);
        fprintf(fparam,
                "ZODIAC_SPECTRUM  = %szodiacal_data/zodiacal_spec.txt\n",
                CADSDATA);
        fprintf(fparam,
                "BKGD_MODEL       = %sdiffuse_data/nuv_allsky.txt\n",
                CADSDATA);
        fprintf(fparam,
                "BKGD_SPECTRUM    = %sdiffuse_data/bkgd_spectrum.txt\n",
                CADSDATA);
        fprintf(fparam,
                "STELLAR_SPECTRUM = %sstellar_data/castelli\n",
                CADSDATA);
        fprintf(fparam,
                "CROSS_SEC        = %sstellar_data/cross_sec.txt\n",
                CADSDATA);
        fprintf(fparam,
                "STAR_CAT         = %sstellar_data/hip_main.txt\n",
                CADSDATA);
        fprintf(fparam,
                "FILTER_SPECTRUM  = %ssample_filters/galex_FUV.txt\n\n"
               );

	fprintf(fparam,
			"# WCS Parameters for output FITS file\n"
            "#GAL_COORD = 1 => output is in galactic coordinates, othrwise equatorial"
			"# ------------------------\n"
			"XREFVAL = 0\n"
			"YREFVAL = 0\n"
			"XINC = -0.5\n"
			"YINC = 0.5\n"
            "ANGSIZE = 360\n"
			"COORDTYPE = -AIT\n"
            "GAL_COORD = 1\n\n"
			);

	fprintf(fparam,
			"# Background parameters\n"
			"# ---------------------\n"
			"ANG_LIMIT = 3\n"
			"CSC_LAW = 450.0\n\n"
            "DC_LEVEL = 0\n\n"
			);

	fprintf(fparam,
			"# Which elements to include\n"
			"# -----------------------\n"
			"ZOD_INC = 1\n"
			"GAL_INC = 1\n"
			"STAR_INC = 1\n\n"
			);

	printf("Default Parameter file has been created. Please check and modify, if necessary.\n");
	fclose(fparam);
	exit(EXIT_FAILURE);

}/* END GEN_DEFAULT_PARAMS */
/*********************************************************************/

int READ_PARAMS(struct INP_PAR *inp_par, struct WCS_HEAD *wcs_in)
/*
                float *hour, float *day, float *month, float *year,
				float *wave_ref, float *ang_limit, float *csc_law,
				int *zod_inc, int *gal_inc,
                char *zod_file, char *spec_file, char *bkgd_file, char *out_file,
				struct WCS_HEAD *wcs_in) {
*/
                {
	FILE    *fp;
    char    line[FLEN_FILENAME], name[FLEN_FILENAME];
    char    equal[FLEN_FILENAME], data[FLEN_FILENAME];

	if((fp = fopen(PARAM_FILE,"r")) == NULL){
		GEN_DEFAULT_PARAMS();
		if((fp = fopen(PARAM_FILE,"r")) == NULL)
			exit(EXIT_FAILURE);
	}
	/*Inititalize parameters*/
    inp_par->hour = 0;
    inp_par->day = 0;
                    inp_par->month = 0;
                    inp_par->csc_law=0;
                    inp_par->ang_limit = 5;
                    inp_par->dc_level= 0;
                    inp_par->gal_inc = 1;
                    inp_par->star_inc = 1;
                    inp_par->zod_inc = 1;

	/* Read input parameters */
	while (!feof(fp)) {
		fgets(line, FLEN_FILENAME, fp);
		if ((line[0] != '#') && (line[0] != '\n')) {
			sscanf(line, "%s %s %s ", name, equal, data);
			if(!strcmp("TIME_UT", name))        inp_par->hour = atof(data);
			if(!strcmp("DAY_OF_MONTH", name))   inp_par->day = atof(data);
			if(!strcmp("MONTH", name))          inp_par->month = atof(data);
			if(!strcmp("YEAR", name))           inp_par->year = atof(data);
			if(!strcmp("OUTPUT_FILE", name))    strcpy(inp_par->out_file, data);
			if(!strcmp("ZODIAC_MODEL", name))      strcpy(inp_par->zod_file, data);
			if(!strcmp("ZODIAC_SPECTRUM", name))   strcpy(inp_par->zod_spec, data);
			if(!strcmp("FILTER_SPECTRUM", name))   strcpy(inp_par->filter_file, data);
            if(!strcmp("BKGD_SPECTRUM", name))   strcpy(inp_par->bkgd_spec, data);
			if(!strcmp("BKGD_MODEL", name))		strcpy(inp_par->bkgd_model, data);
            if(!strcmp("CROSS_SEC", name))		strcpy(inp_par->sigma_file, data);
            if(!strcmp("STELLAR_SPECTRUM", name))
                                        strcpy(inp_par->castelli_file, data);
            if(!strcmp("STAR_CAT", name))
                strcpy(inp_par->hipparcos_file, data);
			if(!strcmp("ANG_LIMIT", name))		inp_par->ang_limit = atof(data);
			if(!strcmp("CSC_LAW", name))		inp_par->csc_law = atof(data);
            if(!strcmp("DC_LEVEL", name))		inp_par->dc_level = atof(data);
            if (!strcmp("ZOD_INC", name))       inp_par->zod_inc = atoi(data);
            if (!strcmp("GAL_INC", name))       inp_par->gal_inc = atoi(data);
            if (!strcmp("STAR_INC", name))      inp_par->star_inc = atoi(data);
			if(!strcmp("XREFVAL", name))        wcs_in->xrefval = atof(data);
			if(!strcmp("YREFVAL", name))		wcs_in->yrefval = atof(data);
			if(!strcmp("XINC", name))			wcs_in->xinc = atof(data);
			if(!strcmp("YINC", name))           wcs_in->yinc = atof(data);
			if(!strcmp("ANGSIZE", name))        wcs_in->angsize = atof(data);
            if(!strcmp("COORDTYPE", name))      strcpy(wcs_in->coordtype, data);

		}
	}
    wcs_in->rot = 0;
	strcpy(line, "End of File");
	/*This is to take care of the end of the file*/
	fclose(fp);
	return(EXIT_SUCCESS);
}/* END READ_PARAMS */

struct INP_PAR READ_INPUT_PARAMETERS(int argc, char *argv[])
{
    int i_arg = 2;
    FILE *para_file=NULL;
    struct INP_PAR inp_par;

    if (PARSE_PAR_CHAR(argc, argv, inp_par.par_file, &i_arg)){
        printf("Parameter file name is required on command line\n.");
        exit(EXIT_FAILURE);
    }/*File which contains names of other files*/

    if((para_file = fopen(inp_par.par_file,"r")) == NULL){
        printf("Could not open file list:\n");
        exit(EXIT_FAILURE);
    }
    strcpy(inp_par.out_dir, "");
    /*Read list of files*/
    fscanf(para_file, "%s",inp_par.castelli_file); /*Location of Castelli data*/
    fscanf(para_file, "%s",inp_par.sigma_file); /* Optical constants */
    fscanf(para_file, "%s",inp_par.hipparcos_file); /*Location of Hipparcos data*/
    fscanf(para_file, "%s",inp_par.location_file); /* Locations for calculations */
    fscanf(para_file, "%s", inp_par.out_dir);   /*Writes all output files in this directory*/
    fclose(para_file);

    /*Output file name is derived from parameters*/
    strcpy(inp_par.out_file,"total_scattered.dat");

    return (inp_par);
}

/**************** SPECT_READ**************************/
/*Reads the spectral (2 column) files*/
int SPECT_READ(float *wave, float *spec, char *spec_file)
{
	FILE *spec_file_ptr = NULL;
	int ispec;
	float f_tmp1, f_tmp2;
	char s_tmp[FLEN_FILENAME];

	if ((spec_file_ptr = fopen(spec_file, "r")) == NULL){
		fprintf(stderr, "ERROR  : Unable to open %s \n", spec_file);
		exit(EXIT_FAILURE);
	}/*Open spectrum if it is exists*/

    ispec =0;
    while (!feof(spec_file_ptr)) {
        if (ispec >  MAX_SPEC_ELEM){
            printf("%s %s\n", "Max. size exceeded for file:",spec_file);
            exit(EXIT_FAILURE);
        }
        /*Read line*/
        fgets(s_tmp, FLEN_FILENAME, spec_file_ptr);
		if (sscanf(s_tmp, "%f %f", &f_tmp1, &f_tmp2) == 2){
            *(wave + ispec) = f_tmp1;
            *(spec + ispec) = f_tmp2;
            ++ispec;
        }
    }
	fclose(spec_file_ptr);
    return (ispec);
}/*END SPECT_READ*/
/*********************************************************************/


