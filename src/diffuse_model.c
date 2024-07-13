/********* Monte Carlo program for Milky Way Scattering***************/

/********************* Modification History:*******************************
 NOTE THAT I WILL ONLY CHANGE THE VERSION NUMBER IF THE RESULTS ARE CHANGED
 Dec. 11, 2014: Corrected fluxes for Almach. Hipparcos said it was a 2nd mag. B8.
 Jan. 10, 2015: Minor modification to remove extraneous calculations
 Jan. 10, 2015: Added a and g to the FITS header
 Jan. 17, 2015: Added options for saving files
 Mar. 8,  2015: Added more complete photon log
 Mar. 21, 2015: Corrected two major errors
 Mar. 22, 2015: Replaced squares by pow(...,2)
 Mar. 26, 2015: Checked a bunch of stuff with minor corrections.
 Apr. 29, 2015: Problem when tau was 0 or 1. Changed to exclude those.
 May   6, 2015: Minor changes in log files
 May   8, 2015: Fixed error in converting string to long
 May   9, 2015: Extra warning messages
 May  14, 2015: Added wavelength to FITS HEADER
 May  24, 2015: Added distance to datalogger
 May  24, 2015: Minor changes in printing output
 May  25, 2015: Corrected printing option
 Aug   3, 2015: Memory leak in printing FITS comment
 Aug. 26, 2015; Modified diagnostics
 Oct. 27, 2015; Possible error in distlog
 Jan. 28, 2016: Used too long keyword for wavelength
 Jul. 1,  2016: Memory leak in star spectral type
 Jan. 31, 2017: Changed random number call to include double precision.
 May  31, 2020: Reject photons if they don't fall within the solid angle of the star to the box.
 May 24, 2023: Corrected error where blank spectral types were incorrectly assigned
 May 24, 2023: Spectral type of HIP 105259 was assigned B1 instead of M1
********************End Modification History************************************/

#include "diffuse_model.h"
#define VERSION "May. 24, 2023 V3.2"
//Windows code
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)

/*****************   Begin Main Program  *******************************/
int main (int argc, char *argv[])
{
    //Input parameters
    struct INP_PAR inp_par;
    struct WCS_HEAD wcs;

    //Star related
    struct STARS    star, *hipstars;
    struct SPECTRA  stellar_spectra[N_CASTELLI_MODELS];
    FILE    *hipfile, *star_file = NULL, *ordered_list = NULL, *extra_star_file;
    double  *star_random, tot_star;
    long    *star_list;
    long     istar;
    int     *exclude_stars;
    int     NSTARS;

    //Dust related
    float   *dust_dist, *dust_arr, *col_density;
    long    dust_index;
    double  sigma;
    double  intens, extinct;
    FILE    *phot_log_file;
     float     x1[8], y1[8], z1[8], min_x, min_y, min_z, max_x, max_y, max_z,phibox[8];

    //SCATTERING
    double  *angle, *fscat, theta, phi, tau;
    long    nphoton, nscatter;
    double  xp, yp, zp, x_new, y_new, z_new;
    double  delta_x, delta_y, delta_z;
    double  cum_tau;
    float   flux, cum_flux;
    int     NEW_PHOTON;
    float   max_tau;
    double  ra, dec;
    double dst;

    //Random number
    uint32_t   init_seed = 1;
    double   random;
    dsfmt_t dsfmt;
    int ran_size;
    int nrandom = 10000000;
    double *ran_array;
    int ran_ctr;

    //Time related
    time_t time1, time2;

    //Other Variables
    int     i, j, tst;
    int     *starlog;
    float   *totlog;
    int     *misslog;
    float   *scatlog, *distlog;
    double  ftmp=0;
    long    ltmp;

    printf("diffuse_model: Version %s\n", VERSION);
    time(&time1);
    printf("Starting at %s", ctime(&time1));

    /***********  Read parameters       ************************/
    inp_par = READ_INPUT_PARAMETERS(argc, argv);
    strcpy(inp_par.version, VERSION);
    dust_dist = DUST_READ(&inp_par, inp_par.dust_file); /*Read dust file*/
    col_density = DUST_READ(&inp_par, inp_par.dust_col_density);//File containing column densities
    wcs = WCS_READ(inp_par.wcs_file);

    //Probabilities for scattering
    angle = (double *) malloc(sizeof(double) * NANGLE);
    fscat = (double *) malloc(sizeof(double) * NANGLE);
    SETUP_ANGLE(angle);//Angles from star
    SETUP_SCATTER(fscat, inp_par.g);

    /************************   Initializion ************************/

    x_new = y_new = z_new = theta = phi = xp = yp = zp = dust_index = tau = nphoton = 0;
    delta_x = delta_y = delta_z = 0;
//random numbers
    int no_random = 0;
    if (no_random == 0) {
//As it turns out, I won't need to take the mod for 100 years but let's plan for the future.
        time1 = pow(2,32);
        time2 = time(0);
        init_seed = (uint32_t) (time2 % time1);
    } else printf("%s\n", "************WARNING: NOT RANDOM");
//    setall(init_seed, init_seed + 1);// Sets seed
    dsfmt_init_gen_rand(&dsfmt, init_seed);
    init_seed = dsfmt_genrand_uint32(&dsfmt);
//A check to make sure we meet the minimum size criterion
    ran_size = dsfmt_get_min_array_size();
    if (ran_size < nrandom)
        ran_size = nrandom;
//Initialize random array
    ran_array = (double *) malloc(sizeof(double)*ran_size);
//    if (posix_memalign((void**) &ran_array, 16, sizeof(double)*nrandom) != 0)
//        printf("Couldn't allocate memory for random numbers\n");
    ran_ctr = nrandom;

    //The dust and energy arrays are continuously built up
    dust_arr = (float *) calloc(wcs.nx*wcs.ny, sizeof(float));

    //Get number of additional stars
    extra_star_file = fopen(inp_par.extra_stars, "r");
    if (extra_star_file != NULL) {
        fscanf(extra_star_file, "%i", &i);
        NSTARS = NHIP_STARS + i;
    } else NSTARS = NHIP_STARS;

    //Arrays related to the stars
    hipstars    = (struct STARS *) malloc(sizeof(struct STARS) * NSTARS);
    star_random = (double *) malloc(sizeof(double) * NSTARS);
    starlog     = (int *)    calloc(NSTARS, sizeof(int));
    distlog     = (float *)  calloc(NSTARS, sizeof(float));
    scatlog     = (float *)  calloc(NSTARS, sizeof(float));
    totlog      = (float *)  calloc(NSTARS, sizeof(float));
    misslog     = (int *)    calloc(NSTARS, sizeof(int));
    exclude_stars = (int *)  calloc(NSTARS, sizeof(int));
    star_list   = (long *)   calloc(NSTARS, sizeof(long));

/********************************* End Inititalization ******************************/

    /*Stellar spectra from Castelli*/
    for (i = 0; i < N_CASTELLI_MODELS; ++i){
        stellar_spectra[i].spectrum = malloc(sizeof(float) * N_CASTELLI_SPECTRA);
        stellar_spectra[i].wavelength = malloc(sizeof(float) * N_CASTELLI_SPECTRA);
    }
    tst = READ_CASTELLI_SPECTRA(inp_par.castelli_file, stellar_spectra);
    /*
     Read  Cross-sections from Draine's file. The Draine results are per H atom.
     The densities in the dust file are scaled to hydrogen atom per cm-3. Because
     we use distances in pc we have to scale the sigma to cm-2 pc-1. This helps
     keep all the numbers manageable.
     */
    sigma = READ_CROSS_SEC(inp_par);

    /******************** Read in the stars from the Hipparcos file *********************/
    hipfile = fopen(inp_par.hipparcos_file, "r");
    if (hipfile == NULL){
        printf("%s %s\n", "Can't find stellar file:", inp_par.hipparcos_file);
        exit(EXIT_FAILURE);
    }
    tot_star = 0;  
    for (i = 0; i < NHIP_STARS; ++i){
        tst = HIP_READ_LINE(hipfile, &star);
        tst = GET_STAR_TEMP(&star);
        tst = GET_SCALE_FACTOR(&star, stellar_spectra, inp_par);
        hipstars[i] = star;
        star_list[i] = i;
        int doprint=0;

        if (hipstars[i].HD_NO == 158408)doprint = 1;
        if (hipstars[i].HD_NO == 158926)doprint = 1;
        if (hipstars[i].HD_NO == 164794)doprint = 1;


        if (doprint == 1)
        printf("%i %li %s %i %f %f %f %f %f %f\n", hipstars[i].HIP_NO, hipstars[i].HD_NO, hipstars[i].sp_type, hipstars[i].temperature,
               hipstars[i].gl, hipstars[i].gb, hipstars[i].distance, hipstars[i].V_mag, hipstars[i].ebv,
               hipstars[i].tot_photons);
    }
    fclose(hipfile);

    //If there are any additional stars, we read them
    for (i = NHIP_STARS; i < NSTARS; ++i){
        tst = EXTRA_READ_LINE(extra_star_file, &star);
        tst = GET_STAR_TEMP(&star);
        tst = GET_SCALE_FACTOR(&star, stellar_spectra, inp_par);
        hipstars[i] = star;
        star_list[i] = i;
        int doprint=0;
        if (doprint == 1)
            printf("%i %li %s %f %f %f %f %f\n", hipstars[i].HIP_NO, hipstars[i].HD_NO, hipstars[i].sp_type,
                   hipstars[i].gl, hipstars[i].gb, hipstars[i].distance, hipstars[i].V_mag,
                   hipstars[i].tot_photons);
    }
    if (extra_star_file != NULL){
        printf("*************** No Additional Stars *************\n");
        fclose(extra_star_file);
    }

    /*
     Exclude stars from list. If the file exists and the star is in the exclusion list
     we set the number of photons from that star to 0. That removes it from consideration.
     */
    star_file = fopen(inp_par.star_file, "r");
    istar = 0;
    if (star_file != NULL) {
        while (feof(star_file) == 0) {
            fscanf(star_file, "%i", &tst);
            exclude_stars[istar] = tst;
            ++istar;
        }
        fclose(star_file);
    } else printf("********** All stars will be included **********\n");
    for (j = 0; j < istar; ++j) {
        for (i = 0; i < NSTARS; ++i) {
            if (hipstars[i].HIP_NO == exclude_stars[j]){
                hipstars[i].tot_photons = 0;
                break;
            }
        }
    }
/*How many photons from each star come close to our FOV? I have to find the limiting angle.
 First find out where the edges of the dust distribution are.
 */
    x1[0] = -inp_par.sun_x*inp_par.dust_binsize;
    y1[0] = -inp_par.sun_y*inp_par.dust_binsize;
    z1[0] = -inp_par.sun_z*inp_par.dust_binsize;
    x1[1] = -inp_par.sun_x*inp_par.dust_binsize;
    y1[1] = (inp_par.dust_ysize - inp_par.sun_y)*inp_par.dust_binsize;
    z1[1] = -inp_par.sun_z*inp_par.dust_binsize;
    x1[2] = (inp_par.dust_xsize -inp_par.sun_x)*inp_par.dust_binsize;
    y1[2] = -inp_par.sun_y*inp_par.dust_binsize;
    z1[2] = -inp_par.sun_z*inp_par.dust_binsize;
    x1[3] = (inp_par.dust_xsize - inp_par.sun_x)*inp_par.dust_binsize;
    y1[3] = (inp_par.dust_ysize - inp_par.sun_y)*inp_par.dust_binsize;
    z1[3] = -inp_par.sun_z*inp_par.dust_binsize;
    x1[4] = -inp_par.sun_x*inp_par.dust_binsize;
    y1[4] = -inp_par.sun_y*inp_par.dust_binsize;
    z1[4] = (inp_par.dust_zsize - inp_par.sun_z)*inp_par.dust_binsize;
    x1[5] = -inp_par.sun_x*inp_par.dust_binsize;
    y1[5] = (inp_par.dust_ysize - inp_par.sun_y)*inp_par.dust_binsize;
    z1[5] = (inp_par.dust_zsize - inp_par.sun_z)*inp_par.dust_binsize;
    x1[6] = (inp_par.dust_xsize -inp_par.sun_x)*inp_par.dust_binsize;
    y1[6] = -inp_par.sun_y*inp_par.dust_binsize;
    z1[6] = (inp_par.dust_zsize - inp_par.sun_z)*inp_par.dust_binsize;
    x1[7] = (inp_par.dust_xsize - inp_par.sun_x)*inp_par.dust_binsize;
    y1[7] = (inp_par.dust_ysize - inp_par.sun_y)*inp_par.dust_binsize;
    z1[7] = (inp_par.dust_zsize - inp_par.sun_z)*inp_par.dust_binsize;
    min_x = min_y = min_z = 1e6;
    max_x = max_y = max_z =-1e6;

    for (i = 0; i < 7; ++i){
        if (x1[i] < min_x) min_x = x1[i];
        if (y1[i] < min_y) min_y = y1[i];
        if (z1[i] < min_z) min_z = z1[i];
        if (x1[i] > max_x) max_x = x1[i];
        if (y1[i] > max_y) max_y = y1[i];
        if (z1[i] > max_z) max_z = z1[i];
    }

    for (i = 0; i < NSTARS; ++i) {
        hipstars[i].min_phi   = 360;
        hipstars[i].max_phi   = 0;
        hipstars[i].min_theta = 180;
        hipstars[i].max_theta = 0;
//First find the theta angle. This is just the ratio of the z over the total distance
        for (j = 0; j < 8; ++j){
            theta = (z1[j] - hipstars[i].z);
            ftmp = sqrt(pow(x1[j] - hipstars[i].x,2) + pow(y1[j] - hipstars[i].y,2) +
                    pow(z1[j] - hipstars[i].z,2));
            theta /= ftmp;
//Theta is the angle from the z axis.
            theta = 90 - RADDEG(asin(theta));
            if (theta > hipstars[i].max_theta)
                hipstars[i].max_theta = theta;
            if (theta < hipstars[i].min_theta)
                   hipstars[i].min_theta = theta;
        }
//If the star is in or over or under the box
        if ((hipstars[i].x > min_x) && (hipstars[i].x < max_x) &&
            (hipstars[i].y > min_y) && (hipstars[i].y < max_y)){
            if (hipstars[i].z > min_z) hipstars[i].max_theta  = 180.;
            if (hipstars[i].z < max_z) hipstars[i].min_theta  = 0.;

        }
//Now the phi angle
        for (j = 0; j < 8; ++j){
           phibox[j] = (x1[j] - hipstars[i].x);
           ftmp = sqrt(pow(x1[j] - hipstars[i].x, 2) + pow(y1[j] - hipstars[i].y, 2));
           phibox[j] /= ftmp;
           phibox[j] = RADDEG(acos(phibox[j]));
        }

//The different cases for where the star is
        if ((hipstars[i].x < min_x) && (hipstars[i].x < max_x) && (hipstars[i].y < min_y) && (hipstars[i].y < max_y)){
                hipstars[i].min_phi  = phibox[2];
                hipstars[i].max_phi  = phibox[1];
        }
        if ((hipstars[i].x < min_x) && (hipstars[i].x < max_x) && (hipstars[i].y > min_y) && (hipstars[i].y < max_y)){
                hipstars[i].min_phi  = 0;
                hipstars[i].max_phi  = 360;
        }
        if ((hipstars[i].x < min_x) && (hipstars[i].x < max_x) && (hipstars[i].y > min_y) && (hipstars[i].y > max_y)){
                hipstars[i].min_phi  = 360 - phibox[0];
                hipstars[i].max_phi  = 360 - phibox[3];
        }
        if ((hipstars[i].x > min_x) && (hipstars[i].x < max_x) && (hipstars[i].y < min_y) && (hipstars[i].y < max_y)){
                hipstars[i].min_phi  = phibox[2];
                hipstars[i].max_phi  = phibox[0];
}
        if ((hipstars[i].x > min_x) && (hipstars[i].x < max_x) && (hipstars[i].y > min_y) && (hipstars[i].y < max_y)){
                hipstars[i].min_phi  = 0;
                hipstars[i].max_phi  = 360;
}
        if ((hipstars[i].x > min_x) && (hipstars[i].x < max_x) && (hipstars[i].y > min_y) && (hipstars[i].y > max_y)){
                hipstars[i].min_phi  = 360 - phibox[1];
                hipstars[i].max_phi  = 360 - phibox[3];
        }
        if ((hipstars[i].x > min_x) && (hipstars[i].x > max_x) && (hipstars[i].y < min_y) && (hipstars[i].y < max_y)){
                hipstars[i].min_phi  = phibox[3];
                hipstars[i].max_phi  = phibox[0];
        }
        if ((hipstars[i].x > min_x) && (hipstars[i].x > max_x) && (hipstars[i].y > min_y) && (hipstars[i].y < max_y)){
                hipstars[i].min_phi  = phibox[3];
                hipstars[i].max_phi  = 360 - phibox[2];
        }
        if ((hipstars[i].x > min_x) && (hipstars[i].x > max_x) && (hipstars[i].y > min_y) && (hipstars[i].y > max_y)){
                hipstars[i].min_phi  = 360 - phibox[1];
                hipstars[i].max_phi  = 360. - phibox[2];
        }
    }



    for (i = 0; i < NSTARS; ++i){
        tot_star += hipstars[i].tot_photons;
    }

    /*
     Setup the array to pick a star. The stars are weighted by the number of photons so
     that brighter stars are more likely to be picked. The numbers go from 0 to the total number
     of photons from the star, which is quite a lot.
     */

    //First reorder the list so that the largest probabilities are at the beginning
    ordered_list = fopen("ordered_list_of_stars.list", "r");
    if (ordered_list == NULL){
        for (i = 0; i < NSTARS - 1; ++i) {
            printf("Sorting stars %i %i\r", i, NSTARS);
            for (j = i + 1; j < NSTARS; ++j) {
                if (hipstars[star_list[j]].tot_photons > hipstars[star_list[i]].tot_photons){
                    SWAP(star_list[i], star_list[j]);
                }
            }
        }
        ordered_list = fopen("ordered_list_of_stars.list", "w");
        for (i = 0; i < NSTARS; ++i){
            fprintf(ordered_list, "%li\n", star_list[i]);
        }
    fclose(ordered_list);
    } else {
        printf("***********  WARNING *************\n");
        printf("READING STAR LIST FROM FILE\n");
        for (i = 0; i < NSTARS; ++i) {
            fscanf(ordered_list, "%li", &ltmp);
            star_list[i] = ltmp;
        }
        fclose(ordered_list);
    }
    star_random[0] = hipstars[star_list[0]].tot_photons;
    for (i = 1; i < NSTARS; ++i)
        star_random[i] = star_random[i-1] + hipstars[star_list[i]].tot_photons;
    /**************************** End Hipparcos star read ******************************/
 
//Finally begin the actual scattering
    nphoton = 0;
    time(&time1);
    printf("Beginning Scattering\n");
    phot_log_file = fopen("every_photon.log", "w");
    int fix_rnd = 0;

    while (nphoton < inp_par.num_photon) {//Total number of photons in the simulation
//        getsd(&iseed1, &iseed2);
        if (fix_rnd == 1){
            printf("***********Warning: Fixed Random Numbers*************\n");
            init_seed = 1645310043;
            fix_rnd = -1;
//            setall(iseed1, iseed2);
        }
        //Save work every 1,000,000 photons
        if (((nphoton % 10000000) == 0) && (nphoton > 0)){
            fclose(phot_log_file); 
            if (strcmp(inp_par.print_debug, "yes") == 0)
                phot_log_file = fopen("every_photon.log", "a");
            else phot_log_file = fopen("every_photon.log", "w");
            CHECKPOINT(dust_arr, inp_par, nphoton, tot_star,
                       wcs, hipstars, starlog, misslog, totlog, distlog, scatlog);
            time(&time2);
            printf("Time taken for loop: %ld\n",time2 - time1);
            time1 = time2;
        }

        NEW_PHOTON  = TRUE;

        /*
         Select star. We have created an array and now just have to pick a random number
         to select the star
         */
        random =  GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, &ran_ctr, &init_seed)*tot_star;

        istar=0;
        while ((star_random[istar] < random) && (istar < NSTARS))
            ++istar;

        istar = star_list[istar];
        star = hipstars[istar];
        ++starlog[istar];
        /*
         We emit a photon from the star. If it is out of the box, we start from the point it
         reaches the box. Note that theta is measured from the z axis; otherwise we
         use Cartesian coordinates. The x and delta_x are in units of bins as
         determined by inp_par.binsize
         */
        tst = FIRST_PHOTON(&x_new, &y_new, &z_new, &delta_x, &delta_y, &delta_z,
                           inp_par, star, angle, dsfmt, ran_array, nrandom, &ran_ctr, &init_seed);

        /*
         From here we don't care which scattering it is (from star or from dust)
         */

        ++nphoton;
        xp = x_new;
        yp = y_new;
        zp = z_new;
        max_tau = 0;
        intens = 1.;
        nscatter = 0;
        cum_flux = 0;

        //If none of the photons enter the box, we have no scattering there.
        if (tst == -1) {
            intens = 0;
            nscatter = -1;
            ++misslog[istar];
        }
        while ((intens >= MIN_INTENS) && (nscatter <= inp_par.nscatter)){
            x_new = xp;
            y_new = yp;
            z_new = zp;

            random = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, &ran_ctr, &init_seed);
            if (random == 1) {
                tau = 1e-10;//Set to some tiny number greater than 0
            }else if (random == 0){
                tau = 100000;//Set to a large number less than infinity
            }else tau = -log(random);

            cum_tau = MATCH_TAU(&x_new, &y_new, &z_new, inp_par,
                                delta_x, delta_y, delta_z,
                                sigma, dust_dist, tau);
            /* We check to see if the photon is still inside the box. If it is
             we can continue.
             */
            if ((CHECK_LIMITS(x_new, y_new, z_new, inp_par) == TRUE) && (cum_tau >= tau)){

                /*
                 Here we scatter the photons. In order to save computing, we always send some
                 part back to the detector. The remaining energy is left in the cell and will
                 heat the ISM.
                 */
                dust_index = GET_DUST_INDEX(x_new, y_new, z_new, inp_par);
                intens *= inp_par.albedo;//We have to scale by the albedo
                /*Part of the photon is scattered back to the detector*/
                flux = DETECT(x_new, y_new, z_new, delta_x, delta_y, delta_z, inp_par);
                /*
                 Finally I have to divide the flux by the distance squared in parsecs so I have to
                 multiply the distance by the size of each bin.
                 */
                dst = CALC_DIST(x_new, y_new, z_new, inp_par.sun_x, inp_par.sun_y, inp_par.sun_z);
                dst *= pow(inp_par.dust_binsize, 2);
                flux *= intens/dst;
                intens -= flux;//We have to account for the small amount that was scattered
                extinct   = col_density[dust_index]*sigma;//Extinction back to the Sun
                flux *= exp(-extinct);
                cum_flux += flux;
                totlog[istar] += flux;//Flux from each star
                ltmp = (long) sqrt(dst);
                distlog[ltmp] += flux;//Keep track of flux as a function of distance
                scatlog[nscatter] += flux;//Flux for each scatter
                /*Convert cartesian back to galactic and project it onto the grid.*/
                CART_TO_CELEST(x_new, y_new, z_new, &ra, &dec, inp_par);
                WRITE_TO_GRID(wcs, ra, dec, flux, dust_arr);
//Only print out selected data
                if ((ra > inp_par.min_gl_debug) && (ra < inp_par.max_gl_debug) &&
                    (dec > inp_par.min_gb_debug) && (dec < inp_par.max_gb_debug))
                fprintf(phot_log_file, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %li %g\n",
                        star.HIP_NO, ra, dec, xp, yp, zp, x_new, y_new, z_new, extinct, nscatter, flux);
                /*
                 Now for the multiple scattering part. The photon is reemitted in
                 a random direction. After that the procedure continues as normal.
                 */
                random = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, &ran_ctr, &init_seed);
                theta = CALC_THETA(fscat, random);
                phi = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, &ran_ctr, &init_seed)*TWOPI;
                xp = x_new;// *New starting point for photon
                yp = y_new;
                zp = z_new;
                SCATTER(&delta_x, &delta_y, &delta_z, theta, phi);
                max_tau = 0;//Reset optical depth
                ++nscatter;

            } else {
/*
 I had originally reused part of the photon but this is actually double counting because
 I don't reduce the intensity of the photons that stay inside. To be consistent, I would
 have to reduce the intensity of both the photons that stay inside the box and those that
 exit. Therefore I just let them escape.
 max_tau = exp(-cum_tau);
 intens *= (1 - max_tau);
 */
                intens =0;
            }//CHECK_LIMITS
        }//Loop for each individual photon
    }//Continue until we reach the photon limit
    time(&time2);
    printf("Time taken for loop: %ld\n",time2 - time1);
    CHECKPOINT(dust_arr, inp_par, nphoton, tot_star,
               wcs, hipstars, starlog, misslog, totlog, distlog, scatlog);


    //Free Variables
    for (i = 0; i < N_CASTELLI_MODELS; ++i){
        free(stellar_spectra[i].spectrum);
        free(stellar_spectra[i].wavelength);
    }
    free(dust_dist);
    free(dust_arr);
    free(col_density);
    free(hipstars);
    free(star_random);
    free(star_list);
    free(misslog);
    free(starlog);
    free(totlog);
    fclose(phot_log_file);
}
/*End program*/
