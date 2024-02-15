
#include "diffuse_model.h"
#include "parse_par.h"

/****************Begin Program READ_INPUT_PARAMETERS*****************************/
/*Program to read parameter file*/
struct INP_PAR READ_INPUT_PARAMETERS(int argc, char *argv[])
{
    FILE    *para_file=NULL;
    struct  INP_PAR inp_par;
    int     i_arg = 2;
    char line[FLEN_FILENAME], name[FLEN_FILENAME], data[FLEN_FILENAME];


    if (PARSE_PAR_CHAR(argc, argv, inp_par.par_file, &i_arg)){
        printf("Parameter file name is required on command line\n.");
        exit(EXIT_FAILURE);
    }//File which contains names of other files

    para_file = fopen(inp_par.par_file, "r");/*Open parameter file*/
    if (para_file == NULL){
        printf("No input file\n");
        exit(EXIT_FAILURE);
    }

    //Default debugging parameters
    strcpy(inp_par.print_debug, "no");
    inp_par.min_gl_debug = 0;
    inp_par.max_gl_debug = 360;
    inp_par.min_gb_debug = -90;
    inp_par.max_gb_debug = 90;
    /*Read parameter list*/
    while (!feof(para_file)) {
        fgets(line, FLEN_FILENAME, para_file);
        if ((line[0] != '#') && (line[0] != '\n')) {
            sscanf(line, "%s %s ", name, data);
            if (!strcmp("CASTELLI_FILE", name))
                strcpy(inp_par.castelli_file, data);
            if (!strcmp("SIGMA_FILE", name))
                strcpy(inp_par.sigma_file, data);
            if (!strcmp("HIPPARCOS_FILE", name))
                strcpy(inp_par.hipparcos_file, data);
            if (!strcmp("DUST_FILE", name))
                strcpy(inp_par.dust_file, data);
            if (!strcmp("COL_DENSITY_FILE", name))
                strcpy(inp_par.dust_col_density, data);
            if (!strcmp("WCS_FILE", name))
                strcpy(inp_par.wcs_file, data);
            if (!strcmp("EXCLUDE_FILE", name))
                strcpy(inp_par.star_file, data);
            if (!strcmp("EXTRA_STAR_FILE", name))
                strcpy(inp_par.extra_stars, data);
            if (!strcmp("NO_OF_PHOTONS", name))
                inp_par.num_photon = strtol(data, NULL, 10);
            if (!strcmp("NO_OF_SCATTER", name))
                inp_par.nscatter = atoi(data);
            if (!strcmp("WAVELENGTH", name))
                inp_par.wave = atof(data);
            if (!strcmp("ALBEDO", name))
                inp_par.albedo = atof(data);
            if (!strcmp("PHASE_FUNCTION", name))
                inp_par.g = atof(data);
            if (!strcmp("PRINT_DEBUG", name))
                strcpy(inp_par.print_debug, data);
            if (!strcmp("MIN_GL_DEBUG", name))
                inp_par.min_gl_debug = atof(data);
            if (!strcmp("MAX_GL_DEBUG", name))
                inp_par.max_gl_debug = atof(data);
            if (!strcmp("MIN_GB_DEBUG", name))
                inp_par.min_gb_debug = atof(data);
            if (!strcmp("MAX_GB_DEBUG", name))
                inp_par.max_gb_debug = atof(data);

        }//If the line is not a comment
        }//Read to the end of the file

    fclose(para_file);
    return (inp_par);

}
/*******************************************************************/

/**************Begin Program DUST_READ****************************/
/*The dust file and the column density file are both 3-d FITS files.*/
float *DUST_READ(struct INP_PAR *inp_par, char filename[MAX_FILE_LENGTH])
{
    fitsfile *fptr;
    float   *dust_arr;
    int     status = 0, anynull=0;
    float   nullval = 0;
    long    naxes[3], ltmp;
    char    comment[FLEN_COMMENT];
    double  dtmp;

    fits_open_file(&fptr, filename, READONLY, &status);
    if (status != 0){
        printf("Could not open %s \n", filename);
        exit(EXIT_FAILURE);
    }/*End if*/

/*Read parameters from dust file: xsize, ysize, zsize, binsize*/
    fits_read_key(fptr, TLONG, "NAXIS1", &ltmp, comment, &status);
    naxes[0] = ltmp;
    fits_read_key(fptr, TLONG, "NAXIS2", &ltmp, comment, &status);
    naxes[1] = ltmp;
    fits_read_key(fptr, TLONG, "NAXIS3", &ltmp, comment, &status);
    naxes[2] = ltmp;
    dust_arr = (float *) malloc(sizeof(float) * naxes[0]*naxes[1]*naxes[2]);
    fits_read_3d_flt(fptr, 1, nullval, naxes[0], naxes[1], naxes[0], naxes[1], naxes[2],
                     dust_arr, &anynull, &status);
    inp_par->dust_xsize = naxes[0];
    inp_par->dust_ysize = naxes[1];
    inp_par->dust_zsize = naxes[2];
    fits_read_key(fptr, TDOUBLE, "CRDELT1", &dtmp, comment, &status);
    inp_par->dust_binsize = dtmp;
    fits_read_key(fptr, TDOUBLE, "CRPIX1", &dtmp, comment, &status);
    inp_par->sun_x = dtmp;
    fits_read_key(fptr, TDOUBLE, "CRPIX2", &dtmp, comment, &status);
    inp_par->sun_y = dtmp;
    fits_read_key(fptr, TDOUBLE, "CRPIX3", &dtmp, comment, &status);
    inp_par->sun_z = dtmp;
    fits_close_file(fptr, &status);
    return(dust_arr);
}
/**************************************************************************/
double READ_CROSS_SEC(struct INP_PAR file_list)
{
    FILE    *fp;
    float   wave, cross_sec, t1, t2, t3;
    char    str[132], *read_res;
    int i = 0;
    struct SPECTRA temp;
    double sigma;

    /*  The following steps are to find out the Absorption CrossSection from
     the cross-section file corresponding to our wavelength of interest. This
     file was taken directly from Draine's page:
     ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all
     */

    if ((fp = fopen(file_list.sigma_file, "r")) == NULL){
        printf("Error opening Cross Section File:::\n");
        exit(EXIT_FAILURE);
    }
    temp.spectrum = malloc(sizeof(float) * N_CROSS_SEC);
    temp.wavelength = malloc(sizeof(float) * N_CROSS_SEC);

    //Skip all lines until we reach -
    strcpy(str, "          ");
    while ((strncmp("-", str, 1) != 0)) {
        read_res = fgets(str, 132, fp);
    }
    /*Find the cross section at the nearest wavelength*/
    i = 0;
    while (feof(fp) == 0){
        read_res = fgets(str, 132, fp);/*Read first bracket*/
        sscanf(str, "%f %f %f %e %e",
               &wave, &t1, &t2, &cross_sec, &t3);
        temp.wavelength[i] = wave*10000;//convert to Å
        temp.spectrum[i]   = cross_sec;
        i = i + 1;
    }

    fclose(fp);

/*We have to reverse the wavelengths because of Draine's format*/
    while (temp.wavelength[i] < file_list.wave)
        --i;
//    printf("Using a wavelength of %lf instead of %lf.\n",
//           temp.wavelength[i], file_list.wave);

/*
 Sigma is extinction cross section per H atom. We will scale everything in terms of atoms/cm^3
 but distances in terms of pc. Therefore we multiply sigma by a scale factor (number of cm in a parsec
 to make everything consistent.
*/
    sigma = temp.spectrum[i]*PCinCM(1);

    free(temp.spectrum);
    free(temp.wavelength);
    return(sigma);
}
/**************************************************************************/

/*****************  Begin WRITE_TO_GRID  *******************************/
/*Convert into tangential projection for later write into FITS file*/
void WRITE_TO_GRID(struct WCS_HEAD wcs, double ra, double dec, double flux,
                   float *grid)
{

    double xout, yout;
    int status = 0;
    long ipixel;

    if (fits_world_to_pix(ra, dec, wcs.xrefval, wcs.yrefval, wcs.xrefpix, wcs.yrefpix,
                          wcs.xinc, wcs.yinc, wcs.rot, wcs.coordtype, &xout, &yout, &status))
        return;
    else if ((xout >= 1) && (xout <= wcs.nx) && (yout >= 1) && (yout <= wcs.ny)){
        --xout;//Correct because FITS standard is to start at 1
        --yout;
        ipixel = (long) xout + ((long) yout) * wcs.nx;
        grid[ipixel] += (float) (flux / fabs(DEGRAD(wcs.xinc) * DEGRAD(wcs.yinc)));
    }
    return;
}/*End program*/
/******************************************************************/
/*****************  Begin WRITE_FITS_FILE  *******************************/
/*Store data into 3 dimensional array so we maintain distance information also*/
void WRITE_FITS_FILE(struct WCS_HEAD wcs_out, float *grid, long nphoton,
                     float tot_star, struct INP_PAR inp_par)
{

    int status = 0;
    long i, j, ipixel;
    fitsfile *fptr;
    char filename[FLEN_FILENAME], str_out[FLEN_KEYWORD], str_tmp[50];
    long axes[2];
    float datamin, datamax;
    float *dust_out;

    sprintf(filename, "%s", "scattered.fits");
    remove(filename);
    fits_create_file(&fptr, filename, &status);

    axes[0] = wcs_out.nx;
    axes[1] = wcs_out.ny;
    dust_out = (float *) malloc(wcs_out.nx*wcs_out.ny*sizeof(float));
    datamin = 1e6;
    datamax = 1e-6;
    for (i = 0; i < wcs_out.nx; ++i) {
        for (j = 0; j < wcs_out.ny; ++j) {
            ipixel = i + j*wcs_out.nx;
            dust_out[ipixel] = grid[ipixel]*tot_star/nphoton;

            if (dust_out[ipixel] > datamax)
                datamax = dust_out[ipixel];
            if (dust_out[ipixel] < datamin)
                datamin = dust_out[ipixel];
        }
    }
    /*HDU keywords*/
    fits_create_img(fptr, FLOAT_IMG, 2, axes, &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &wcs_out.xrefval, "REF POINT VALUE IN DEGREES", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &wcs_out.xrefpix, "REF POINT PIXEL LOCATION", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", &wcs_out.xinc, "DEGREES PER PIXEL", &status);
    fits_write_key(fptr, TDOUBLE, "CROTA1", &wcs_out.rot, "ROTATION FROM ACTUAL AXIS", &status);
    strcpy(str_out,"GLON");
    strcat(str_out, wcs_out.coordtype);
    fits_write_key(fptr, TSTRING, "CTYPE1", &str_out, "COORDINATE TYPE", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &wcs_out.yrefval, "REF POINT VALUE IN DEGREES", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &wcs_out.yrefpix, "REF POINT PIXEL LOCATION", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", &wcs_out.yinc, "DEGREES PER PIXEL", &status);
    fits_write_key(fptr, TDOUBLE, "CROTA2", &wcs_out.rot, "ROTATION FROM ACTUAL AXIS", &status);
    strcpy(str_out,"GLAT");
    strcat(str_out, wcs_out.coordtype);
    fits_write_key(fptr, TSTRING, "CTYPE2", &str_out, "COORDINATE TYPE", &status);
    fits_write_key(fptr, TFLOAT, "DATAMIN", &datamin, "Minimum Data Pixel", &status);
    fits_write_key(fptr, TFLOAT, "DATAMAX", &datamax, "Maximum Data Pixel", &status);
    fits_write_key(fptr, TLONG, "NPHOT", &nphoton, "No. of Photon", &status);
    fits_write_key(fptr, TDOUBLE, "ALBEDO", &inp_par.albedo, "Simulation Albedo", &status);
    fits_write_key(fptr, TDOUBLE, "G", &inp_par.g, "Simulation G", &status);
    fits_write_key(fptr, TDOUBLE, "WAVELENG", &inp_par.wave, "Simulation Wavelength", &status);
    strcpy(str_out, "Dust file: ");
    if (strlen(inp_par.dust_file) < 40)
        strcat(str_out, inp_par.dust_file);
    else {
        strncpy(str_tmp, inp_par.dust_file+strlen(inp_par.dust_file) - 30, 39);
        strcat(str_out, str_tmp);
    }
    fits_write_comment(fptr, str_out, &status);
    fits_write_comment(fptr, inp_par.version, &status);
    /*Write image data into first HDU */
    fits_write_2d_flt(fptr, 1, axes[0], axes[0], axes[1], dust_out, &status);
    /*Close file*/
    fits_close_file(fptr, &status);
    free(dust_out);
}/*End write FITS file*/

/******************************** WRITE_ENERGY_FILE *********************************/
void WRITE_ENERGY_FILE(float *energy_arr, struct INP_PAR inp_par, char filename[FLEN_FILENAME],
                       long nphoton, float tot_star)
{
    int status = 0;
    fitsfile *fptr;
    long axes[3], i, j, k, ipixel;
    float *energy_out;

    energy_out = (float *) calloc(inp_par.dust_xsize*inp_par.dust_ysize*inp_par.dust_zsize,
                                  sizeof(float));


    for (i = 0; i < inp_par.dust_xsize; ++i){
        for (j = 0; j < inp_par.dust_ysize; ++j) {
            for (k = 0; k < inp_par.dust_zsize; ++k){
                ipixel = GET_DUST_INDEX(i, j, k, inp_par);
                energy_out[ipixel]  = energy_arr[ipixel]/nphoton;
            }
        }
    }
    remove(filename);
    axes[0] = inp_par.dust_xsize;
    axes[1] = inp_par.dust_ysize;
    axes[2] = inp_par.dust_zsize;

    fits_create_file(&fptr, filename, &status);
    fits_create_img(fptr, FLOAT_IMG, 3, axes, &status);
    fits_write_3d_flt(fptr, 1, inp_par.dust_xsize, inp_par.dust_ysize,
                      inp_par.dust_xsize, inp_par.dust_ysize, inp_par.dust_zsize,
                      energy_out, &status);
    fits_close_file(fptr, &status);
    free(energy_out);
 }
/********************************* CHECKPOINT *********************************/
void  CHECKPOINT(float *dust_arr, struct INP_PAR inp_par, long nphoton,
                 float tot_star, struct WCS_HEAD wcs, struct STARS *hipstars,
                 int *starlog, int *misslog, float *totlog, float *distlog,
                 float *scatlog)
{
    //Time related
    time_t ltime;
    FILE *logfile;
    long i;

    time(&ltime);
    printf("Checkpoint of %li photons at %s", nphoton, ctime(&ltime));

    /*
     Add to cumulative grids. Scale by the number of photons from the star over
     the number of photons in the simulation. Write them out to FITS files.
     */

    WRITE_FITS_FILE(wcs, dust_arr, nphoton, tot_star, inp_par);

    logfile = fopen("datalogger.txt", "w");
    fprintf(logfile, "HIP_NO Dist.  Star_flux Star_phot Miss_phot Dist_sca scat_flux tot_flux\n");
    for (i = 0; i < NHIP_STARS; ++i)
        fprintf(logfile, "%i %f %e %i %i %e %e %e\n", hipstars[i].HIP_NO, hipstars[i].distance,
                hipstars[i].tot_photons, starlog[i], misslog[i], distlog[i], scatlog[i], totlog[i]);
    fclose(logfile);

}
