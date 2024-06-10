/*
  stars.c
  sky_model

  Created by Jayant Murthy on 07/10/12.
  Copyright (c) 2012 Jayant Murthy. All rights reserved.
*/

#include "sky_model.h"
/*Function to read single line from current position in stars file*/

int HIP_READ_LINE(FILE *star_file, struct STARS *line)
{
    int i, tst;
    char hip_line[HIP_MAIN_REC_LEN], *token, *read_res, s_tmp[13], s_tmp1[13];
    double  gl, gb;

    read_res = fgets(hip_line, HIP_MAIN_REC_LEN, star_file);/*Read first bracket*/
    if (read_res == NULL)
        return(EXIT_FAILURE);
    token = strtok(hip_line, HIP_DELIM); /*Set up tokens*/
    for (i = 1; i < 77; ++i) {
        token = strtok(NULL, HIP_DELIM); /*read token*/
        switch (i) {
            case 1:  line->HIP_NO   = atoi(token);
                break;
            case 8:  line->ra       = atof(token);
                break;
            case 9:  line->dec      = atof(token);
                break;
            case 11: if (atof(token) > 0)
                line->distance = 1000./atof(token);
            else line->distance = 1.e6;
                break;
            case 32: line->B_mag    = atof(token);
                break;
            case 34: line->V_mag    = atof(token);
                break;
            case 71: line->HD_NO    = atoi(token);
                break;
            case 76: strcpy(s_tmp, token);
                sscanf(s_tmp, "%s", s_tmp1);
                strcpy(line->sp_type, s_tmp1);
                break;
        }/*End switch*/
        tst = CONV_EQ_TO_GAL(line->ra, line->dec, &gl, &gb);
        line->gl = gl;
        line->gb = gb;

        line->x = cos(DEGRAD(line->gl))*cos(DEGRAD(line->gb))*line->distance;
        line->y = sin(DEGRAD(line->gl))*cos(DEGRAD(line->gb))*line->distance;
        line->z = sin(DEGRAD(line->gb))*line->distance;
    }/*End for*/
    return(0);
}/*End READ_LINE*/
/*********************************************************************/
int READ_CASTELLI_SPECTRA(char *spec_dir, struct SPECTRA
                          *stellar_spectra)
{
    int status = 0, anynull;
    fitsfile *fptr;
    char filename[MAX_FILE_LENGTH];
    int i;
    int temper[N_CASTELLI_MODELS], gindex[N_CASTELLI_MODELS];
    float nullval;
    char stemper[6];

    /*Lookup table for temperatures and g*/
    for (i = 0; i <= 37; ++i)
        temper[i] = 50000 - i*1000;
    for (i = 38; i < 76; ++i)
        temper[i] = 13000 - (i - 37)*250;
    for (i = 0; i <= 4; ++i)
        gindex[i] = 12;
    for (i = 5; i <= 10; ++i)
        gindex[i] = 11;
    for (i = 11; i <= 63; ++i)
        gindex[i] = 10;
    for (i = 64; i <= 75; ++i)
        gindex[i] = 11;

    /*Read stellar spectra*/
    for (i = 0; i < N_CASTELLI_MODELS; ++i) {
        strcpy(filename, spec_dir);
        strcat(filename, "/ckp00_");
        sprintf(stemper, "%i", temper[i]);
        strcat(filename, stemper);
        strcat(filename, ".fits");
        sprintf(stellar_spectra[i].filename, "%i", temper[i]);
        fits_open_file(&fptr, filename, READONLY, &status);
        fits_movabs_hdu(fptr, 2, NULL, &status);
        fits_read_col(fptr, TFLOAT, gindex[i], 1, 1, N_CASTELLI_SPECTRA,
                      &nullval, stellar_spectra[i].spectrum, &anynull,
                      &status);
        fits_read_col(fptr, TFLOAT, 1, 1, 1, N_CASTELLI_SPECTRA,
                      &nullval, stellar_spectra[i].wavelength, &anynull,
                      &status);

        fits_close_file(fptr, &status);
    }

    return(EXIT_SUCCESS);
}

int GET_STAR_TEMP(struct STARS *hipline)
{
    char sptype[10];
    char *str=sptype;

    if (strncmp(hipline->sp_type, "sd", 2) == 0) {
        strcpy(sptype, hipline->sp_type);
        str = str + 2;
        strcpy(hipline->sp_type, str);
    }
    if (strncmp(hipline->sp_type, "O3", 2) == 0) {
        hipline->temperature = 5;
    }
    else if (strncmp(hipline->sp_type, "O4", 2) == 0) {
        hipline->temperature = 7;
    }
    else if (strncmp(hipline->sp_type, "O5", 2) == 0) {
        hipline->temperature = 10;
    }
    else if (strncmp(hipline->sp_type, "O6", 2) == 0) {
        hipline->temperature = 11;
    }
    else if (strncmp(hipline->sp_type, "O7", 2) == 0) {
        hipline->temperature = 13;
    }
    else if (strncmp(hipline->sp_type, "O8", 2) == 0) {
        hipline->temperature = 15;
    }
    else if (strncmp(hipline->sp_type, "O9", 2) == 0) {
        hipline->temperature = 18;
    }
    else if (strncmp(hipline->sp_type, "B0", 2) == 0) {
        hipline->temperature = 20;
    }
    else if (strncmp(hipline->sp_type, "B1", 2) == 0) {
        hipline->temperature = 25;
    }
    else if (strncmp(hipline->sp_type, "B2", 2) == 0) {
        hipline->temperature = 28;
    }
    else if (strncmp(hipline->sp_type, "B3", 2) == 0) {
        hipline->temperature = 31;
    }
    else if (strncmp(hipline->sp_type, "B4", 2) == 0) {
        hipline->temperature = 33;
    }
    else if (strncmp(hipline->sp_type, "B5", 2) == 0) {
        hipline->temperature = 35;
    }
    else if (strncmp(hipline->sp_type, "B6", 2) == 0) {
        hipline->temperature = 36;
    }
    else if (strncmp(hipline->sp_type, "B7", 2) == 0) {
        hipline->temperature = 37;
    }
    else if (strncmp(hipline->sp_type, "B8", 2) == 0) {
        hipline->temperature = 41;
    }
    else if (strncmp(hipline->sp_type, "B9", 2) == 0) {
        hipline->temperature = 49;
    }
    else if (strncmp(hipline->sp_type, "A0", 2) == 0) {
        hipline->temperature = 51;
    }
    else if (strncmp(hipline->sp_type, "A1", 2) == 0) {
        hipline->temperature = 52;
    }
    else if (strncmp(hipline->sp_type, "A2", 2) == 0) {
        hipline->temperature = 53;
    }
    else if (strncmp(hipline->sp_type, "A3", 2) == 0) {
        hipline->temperature = 56;
    }
    else if (strncmp(hipline->sp_type, "A4", 2) == 0) {
        hipline->temperature = 56;
    }
    else if (strncmp(hipline->sp_type, "A5", 2) == 0) {
        hipline->temperature = 56;
    }
    else if (strncmp(hipline->sp_type, "A6", 2) == 0) {
        hipline->temperature = 57;
    }
    else if (strncmp(hipline->sp_type, "A7", 2) == 0) {
        hipline->temperature = 58;
    }
    else if (strncmp(hipline->sp_type, "A8", 2) == 0) {
        hipline->temperature = 59;
    }
    else if (strncmp(hipline->sp_type, "A9", 2) == 0) {
        hipline->temperature = 60;
    }
    else if (strncmp(hipline->sp_type, "F0", 2) == 0) {
        hipline->temperature = 60;
    }
    else if (strncmp(hipline->sp_type, "F1", 2) == 0) {
        hipline->temperature = 63;
    }
    else if (strncmp(hipline->sp_type, "F2", 2) == 0) {
        hipline->temperature = 61;
    }
    else if (strncmp(hipline->sp_type, "F3", 2) == 0) {
        hipline->temperature = 62;
    }
    else if (strncmp(hipline->sp_type, "F4", 2) == 0) {
        hipline->temperature = 62;
    }
    else if (strncmp(hipline->sp_type, "F5", 2) == 0) {
        hipline->temperature = 63;
    }
    else if (strncmp(hipline->sp_type, "F6", 2) == 0) {
        hipline->temperature = 63;
    }
    else if (strncmp(hipline->sp_type, "F7", 2) == 0) {
        hipline->temperature = 63;
    }
    else if (strncmp(hipline->sp_type, "F8", 2) == 0) {
        hipline->temperature = 64;
    }
    else if (strncmp(hipline->sp_type, "F9", 2) == 0) {
        hipline->temperature = 65;
    }
    else if (strncmp(hipline->sp_type, "G0", 2) == 0) {
        hipline->temperature = 65;
    }
    else if (strncmp(hipline->sp_type, "G1", 2) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "G2", 2) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "G3", 2) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "G4", 2) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "G5", 2) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "G6", 2) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "G7", 2) == 0) {
        hipline->temperature = 67;
    }
    else if (strncmp(hipline->sp_type, "G8", 2) == 0) {
        hipline->temperature = 67;
    }
    else if (strncmp(hipline->sp_type, "G9", 2) == 0) {
        hipline->temperature = 67;
    }
    else if (strncmp(hipline->sp_type, "K0", 2) == 0) {
        hipline->temperature = 68;
    }
    else if (strncmp(hipline->sp_type, "K1", 2) == 0) {
        hipline->temperature = 69;
    }
    else if (strncmp(hipline->sp_type, "K2", 2) == 0) {
        hipline->temperature = 70;
    }
    else if (strncmp(hipline->sp_type, "K3", 2) == 0) {
        hipline->temperature = 70;
    }
    else if (strncmp(hipline->sp_type, "K4", 2) == 0) {
        hipline->temperature = 71;
    }
    else if (strncmp(hipline->sp_type, "K5", 2) == 0) {
        hipline->temperature = 72;
    }
    else if (strncmp(hipline->sp_type, "K6", 2) == 0) {
        hipline->temperature = 72;
    }
    else if (strncmp(hipline->sp_type, "K7", 2) == 0) {
        hipline->temperature = 73;
    }
    else if (strncmp(hipline->sp_type, "K8", 2) == 0) {
        hipline->temperature = 73;
    }
    else if (strncmp(hipline->sp_type, "K9", 2) == 0) {
        hipline->temperature = 73;
    }
    else if (strncmp(hipline->sp_type, "M0", 2) == 0) {
        hipline->temperature = 74;
    }
    else if (strncmp(hipline->sp_type, "M1", 2) == 0) {
        hipline->temperature = 74;
    }
    else if (strncmp(hipline->sp_type, "M2", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M3", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M4", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M5", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M6", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M7", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M8", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "M9", 2) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "O", 1) == 0) {
        hipline->temperature = 13;
    }
    else if (strncmp(hipline->sp_type, "B", 1) == 0) {
        hipline->temperature = 35;
    }
    else if (strncmp(hipline->sp_type, "M", 1) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "C", 1) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "A", 1) == 0) {
        hipline->temperature = 56;
    }
    else if (strncmp(hipline->sp_type, "R", 1) == 0) {
        hipline->temperature = 75;
    }
    else if (strncmp(hipline->sp_type, "G", 1) == 0) {
        hipline->temperature = 66;
    }
    else if (strncmp(hipline->sp_type, "W", 1) == 0) {
        hipline->temperature = 35;
    }
    else if (strncmp(hipline->sp_type, "K", 1) == 0) {
        hipline->temperature = 72;
    }
    else if (strncmp(hipline->sp_type, "N", 1) == 0) {
        hipline->temperature = 72;
    }
    else if (strncmp(hipline->sp_type, "S", 1) == 0) {
        hipline->temperature = 72;
    }
    else if (strncmp(hipline->sp_type, "F", 1) == 0) {
        hipline->temperature = 63;
    }
    else if (strncmp(hipline->sp_type, "DA", 2) == 0) {
        hipline->temperature = 35;
    }
    /*Special Cases*/
    else {
        printf("Unimplemented spectral type: %s\n", hipline->sp_type);
        hipline->temperature = 66;
    }
    return (EXIT_SUCCESS);
}
/**************************************************************************/

int READ_CROSS_SEC(struct SPECTRA line, struct INP_PAR file_list)
{
    FILE    *fp;
    float   wave, cross_sec, t1, t2, t3;
    char    str[132], *read_res;
    int i = 0, j;
    struct SPECTRA temp;

    /*  The following steps are to find out the Absorption CrossSection from
     the cross-section file corresponding to our wavelength of interest. This
     file was taken directly from Draine's page:
     ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all
     */

    if ((fp = fopen(file_list.sigma_file, "r")) == NULL){
        printf("Error opening Cross Section File:::\n");
        exit(EXIT_FAILURE);
    }
    temp.spectrum = malloc(sizeof(float) * MAX_SPEC_ELEM);
    temp.wavelength = malloc(sizeof(float) * MAX_SPEC_ELEM);

    /*Skip all lines until we reach -*/
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
        temp.wavelength[i] = wave*10000;/*convert to Å*/
        temp.spectrum[i]   = cross_sec;
        i = i + 1;
    }
    /*We have to reverse the wavelengths because of Draine's format*/
    for (j=0; j < i; ++j) {
        line.wavelength[j] = temp.wavelength[i-j-1];
        line.spectrum[j] = temp.spectrum[i-j-1];
    }
    fclose(fp);
    free(temp.spectrum);
    free(temp.wavelength);
    return (i);
}
/**************************************************************************/

float CALC_FLUX(struct STARS *hipstars, struct INP_PAR files,
                struct SPECTRA cross_sec,
                struct SPECTRA *stellar_spectra,
                struct SPECTRA filters,
                float gl, float gb, float dist, float lambda,
                float ave_sigma)
{
    int istars, windex, sindex, cindex;
    float x, y, z;
    float dstarsq, tau, flux = 0;
    struct SPECTRA star_flux;

    /*Distance from star to point of interest*/
    x = cos(DEGRAD(gl))*cos(DEGRAD(gb))*dist;
    y = sin(DEGRAD(gl))*cos(DEGRAD(gb))*dist;
    z = sin(DEGRAD(gb))*dist;

    star_flux.wavelength = (float *) malloc(sizeof(float) * N_CASTELLI_SPECTRA);
    star_flux.spectrum   = (float *) malloc(sizeof(float) * N_CASTELLI_SPECTRA);
    star_flux.nelements  = N_CASTELLI_SPECTRA;

    istars = 0;
    dstarsq = pow((hipstars[istars].x - x),2) +
    pow((hipstars[istars].y - y),2) +
    pow((hipstars[istars].z - z),2);
    sindex = hipstars[istars].temperature;
    cindex = 0;
    for (windex = 0; windex < star_flux.nelements; ++windex) {
        star_flux.wavelength[windex] = stellar_spectra[sindex].wavelength[windex];
        while ((cross_sec.wavelength[cindex] < star_flux.wavelength[windex])
               && (cindex < cross_sec.nelements))
            ++cindex;
        tau = cross_sec.spectrum[cindex]*sqrt(dstarsq)*ave_sigma*PCinCM(1.);
        star_flux.spectrum[windex] =
        stellar_spectra[sindex].spectrum[windex] * expf(-tau) *
        hipstars[istars].scale/dstarsq*
        1e9*stellar_spectra[sindex].wavelength[windex]/6.626/3.;
    }
    flux = CALC_SCALE_FACTOR(filters, star_flux);
    return(flux);
}

int GET_SCALE_FACTOR(struct STARS *hipstars,
                     struct SPECTRA *stellar_spectra)
{
    float scale;
    int bindex, vindex, sindex, doprint;
    float ebv, b_mag, v_mag, bflux, vflux, bv;

    if (hipstars->V_mag == 0){
        hipstars->scale = 0;/* A few stars are bad*/
    } else {
        sindex = hipstars->temperature;
        /*We scale at 4400 and 5500 Å*/
        vindex = 0;
        while (stellar_spectra[sindex].wavelength[vindex] < 5500)
            ++vindex;
        bindex = 0;
        while (stellar_spectra[sindex].wavelength[bindex] < 4400)
            ++bindex;

        bflux = stellar_spectra[sindex].spectrum[bindex];
        vflux = stellar_spectra[sindex].spectrum[vindex];
        b_mag = -2.5 * log10(bflux/6.61);/*Convert flux into magnitudes*/
        v_mag = -2.5 * log10(vflux/3.64);
        bv = hipstars->B_mag - hipstars->V_mag;
        ebv = bv - (b_mag - v_mag);
        if (ebv < 0) ebv = 0;
        hipstars->ebv = ebv;
        /*Scale factor to convert model to Earth with no extinction*/
        scale = 3.64e-9*pow(10, -0.4 * (hipstars->V_mag - 3.2*ebv))/
        vflux;
        /*Scale to a distance of 1 pc*/
        hipstars->scale = scale*hipstars->distance*hipstars->distance;

        doprint = 0;
        if (doprint > 0)
            printf("%li %s %s %i %10.3e %10.3e\n", hipstars->HIP_NO, hipstars->sp_type,
                   stellar_spectra[sindex].filename, sindex, bflux, vflux);
    }

    return(EXIT_SUCCESS);
}

int READ_HIPPARCOS_CAT(struct INP_PAR inp_par, struct SPECTRA *stellar_spectra,
                       struct STARS *hipstars, struct SPECTRA filter)
{
    FILE *hipfile_ptr=NULL;
    int i, tst;
    struct STARS hipline;
    struct SPECTRA cross_sec;

    tst = READ_CASTELLI_SPECTRA(inp_par.castelli_file, stellar_spectra);
    cross_sec.spectrum = malloc(sizeof(float) * MAX_SPEC_ELEM);
    cross_sec.wavelength = malloc(sizeof(float) * MAX_SPEC_ELEM);
    cross_sec.nelements = READ_CROSS_SEC(cross_sec, inp_par);

    hipfile_ptr = fopen(inp_par.hipparcos_file, "r");
    if (hipfile_ptr == NULL){
        printf("%s %s\n", "Can't find stellar file:", inp_par.hipparcos_file);
        exit(EXIT_FAILURE);
    }
    i = 0;
    while (feof(hipfile_ptr) == 0) {
        tst   = HIP_READ_LINE(hipfile_ptr, &hipline);
        tst = GET_STAR_TEMP(&hipline);
        tst = GET_SCALE_FACTOR(&hipline, stellar_spectra);
        hipline.flux_at_earth = CALC_FLUX(&hipline, inp_par, cross_sec, stellar_spectra, filter,
                                          0, 0, 0, 2300, 1);
        hipstars[i] = hipline;
        ++i;
    }
    hipstars[0].nstars = i-1;
    fclose(hipfile_ptr);
    return (EXIT_SUCCESS);

}