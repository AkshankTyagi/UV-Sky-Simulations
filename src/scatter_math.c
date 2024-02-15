

#include "diffuse_model.h"

/*****************Begin PHASE_FUNCTION********************************/
double PHASE_FUNCTION(double cangle, double g)
{
    double pscat;
/* Henyey-Greenstein function*/
    
    pscat = 1./4./PI * (1 - g*g)/
    (pow(1. + g*g - 2.*g*cangle, 1.5));
    return(pscat);
//CHECKED ON OCT. 17, 2014
}
/*************** End PHASE_FUNCTION ****************************/
/*****************Begin SETUP_ANGLE********************************/
/*
  Note that the probability is proportional to sin(angle) where the
  angle is defined from the z axis (that is different from general
  astronomical use where the angle is defined from the equator.
 */
void SETUP_ANGLE(double *angle)
{
    
    long i;
    double theta;
    
    angle[0] = 0;
    for (i = 1; i < NANGLE; ++i){/*Cumulative probability*/
        theta = (double) (((double) i)/((double) NANGLE)*PI);
        angle[i] = angle[i-1] + sin(theta);
    }/*Endfor*/
    for (i = 1; i < NANGLE; ++i)
        angle[i] = angle[i]/angle[NANGLE - 1];
//CHECKED ON OCT. 17, 2014
}
/******************************************************************/
/*****************Begin SETUP_SCATTER********************************/
/*Use the Henyey-Greenstein function */
void SETUP_SCATTER(double *fscat, double g)
{
    
    long i;
    double theta;
    
    fscat[0] = 0;
    
    for (i = 1; i < NANGLE; ++i){/*Cumulative probability*/
        theta = (double) (((double) i)/((double) NANGLE)*PI);
        fscat[i] = fscat[i-1] + PHASE_FUNCTION(cos(theta), g);
    }/*Endfor*/
    for (i = 1; i < NANGLE; ++i)
        fscat[i] /= fscat[NANGLE - 1];
//CHECKED ON OCT. 17, 2014
}
/******************************************************************/
/*****************Begin CALC_THETA********************************/
/*Calculate theta doing interpolation between steps*/
double CALC_THETA(double angle[NANGLE], float random)
{
    long i = 0, max_i=NANGLE, min_i=0;
    double index, theta;
    
//Binary Search
    while ((max_i - min_i) > 3) {
        i = min_i + (max_i - min_i)/2;
        while ((angle[i] < random) && ((max_i - i) > 1)){
            min_i = i;
            i += (max_i - i)/2;
        }
        if (angle[i] > random)
            max_i = i;
        while ((angle[i] > random) && ((i - min_i) > 1)){
            max_i = i;
            i -= (i - min_i)/2;
        }
        if (angle[i] < random)
            min_i = i;
    }
    i = min_i;
    while ((angle[i] < random) && (i < NANGLE))
        ++i;/*Theta selection*/
    index = ((double) (i - 1)) + (angle[i] - random)/(angle[i] - angle[i - 1]);
    theta = index/((double) NANGLE) * PI;/*Theta will be between 0 and PI in radians*/
    return(theta);
//Checked on Oct. 17, 2014. Gives a proper sin distribution.
}
/******************************************************************/
/*****************Begin DETECT********************************/
/*Some fraction of every photon is scattered to the detector. We work in energy
 units. The fraction that hits the detector is I * PHASE_FUNCTION/ distance^2 * da
 where da is the area of the detector and the 4*pi is included in the phase function.
 Eventually we will have to divide by the solid angle of the detector but that comes
 at the end when we write out the data. Note that all distances are in parsecs as we
 cancel out all scale factors.
*/
double DETECT(double x1, double y1, double z1,
              double dx,   double dy,   double dz,
              struct INP_PAR dust)
{
    double pscat, d1, cangle, d2, x2, y2, z2;

/*
 All distances are in units of binsize.
 */
    d1 = sqrt(dx*dx + dy*dy + dz*dz);
    x2 = dust.sun_x - x1;
    y2 = dust.sun_y - y1;
    z2 = dust.sun_z - z1;
    d2 = sqrt(x2*x2 + y2*y2 + z2*z2);
/*
 The angle is between the direction of motion and the Sun where the angle will be
 0 if the Sun is in the same direction as dx, dy, dz. The two vectors are therefore
 dx, dy, dz and (Sun - point). I only need the cosine of the angle.
*/
    cangle = (dx*x2 + dy*y2 + dz*z2)/d1/d2;
    pscat = PHASE_FUNCTION(cangle, dust.g);
    return(pscat);
}
/******************************************************************/
/*****************   Begin SCATTER  *******************************/
void SCATTER(double *delta_x, double *delta_y, double *delta_z, double theta, double phi)
/*Change frame of reference so that motion is in the z direction and then rotate back*/
{
    double xnew, ynew, znew;
    double dist, dxy, cbeta, sbeta, calpha, salpha, rot_mat[3][3];
    
/*
 The scattering is easiest when the direction of motion is taken as the z axis. We calculate
 the motion in this frame. Note that theta is the angle from the z axis hence
 znew = cos(theta)
*/
    xnew = cos(phi) * sin(theta);
    ynew = sin(phi) * sin(theta);
    znew = cos(theta);
    
/*
 In order to convert back to the actual direction of motion (delta_x, delta_y and delta_z
 we have to define a rotation matrix by first rotating around the z axis by an angle alpha
 and then the y axis by an angle beta. beta is defined from the z axis and alpha from the 
 x axis.
*/
    dist = sqrt(pow(*delta_x,2) + pow(*delta_y, 2) + pow(*delta_z, 2));
    dxy  = sqrt(pow(*delta_x,2) + pow(*delta_y, 2));
    cbeta = (*delta_z)/dist;
    sbeta = dxy/dist;   //Because beta is the angle from the z axis, sbeta is always > 0
    calpha = (*delta_x)/dxy; //By definition alpha is the angle from the x axis
    salpha = (*delta_y)/dxy;
/* We can now set up thr rotation matrix. The elements are defined as [column, row]*/
    rot_mat[0][0] = calpha*cbeta;
    rot_mat[1][0] = -salpha;
    rot_mat[2][0] = calpha*sbeta;
    rot_mat[0][1] = cbeta*salpha;
    rot_mat[1][1] = calpha;
    rot_mat[2][1] = sbeta*salpha;
    rot_mat[0][2] = -sbeta;
    rot_mat[1][2] = 0;
    rot_mat[2][2] = cbeta;
    
    /*Multiply by rotation matrix to give new direction of travel*/
    *delta_x = rot_mat[0][0]*xnew + rot_mat[1][0]*ynew + rot_mat[2][0]*znew;
    *delta_y = rot_mat[0][1]*xnew + rot_mat[1][1]*ynew + rot_mat[2][1]*znew;
    *delta_z = rot_mat[0][2]*xnew + rot_mat[1][2]*ynew + rot_mat[2][2]*znew;
}

/*****************   Begin CHECK_LIMITS  *******************************/
int CHECK_LIMITS(double x, double y, double z, struct INP_PAR dust)
{
    
    if ((x > 0) && (x < (dust.dust_xsize - 1)) &&
        (y > 0) && (y < (dust.dust_ysize - 1)) &&
        (z > 0) && (z < (dust.dust_zsize - 1)))
        return(TRUE);
    else return(FALSE);
}
/******************************************************************/

/*****************  Begin INCREMENT  *******************************/
void INCREMENT(double *x_new,  double *y_new,  double *z_new,
               double delta_x, double delta_y, double delta_z)
{
    (*x_new) += delta_x;
    (*y_new) += delta_y;
    (*z_new) += delta_z;
}
/******************************************************************/
/*****************  Begin CALC_DELTA_X  *******************************/
void CALC_DELTA_X(double *delta_x, double *delta_y, double *delta_z,
                  double theta, double phi)
/*
 Only for the first case where the photon is ejected from the star.
 Note that theta is the colatitude (angle from the North pole, as in the 
 probability.
*/
{
    
    (*delta_x) = cos(phi) * sin(theta);
    (*delta_y) = sin(phi) * sin(theta);
    (*delta_z) = cos(theta);
}
/******************************************************************/
/*****************  Begin GET_DUST_INDEX  *******************************/
long GET_DUST_INDEX(double x_new, double y_new, double z_new, struct INP_PAR dust)
{
   long dust_index;

    dust_index = ((long) (x_new + 0.5)) +
    ((long) (y_new + 0.5)*dust.dust_xsize) +
    ((long) (z_new + 0.5)*dust.dust_xsize * dust.dust_ysize);
    return(dust_index);
}
/******************************************************************/

/**************************FIRST_PHOTON*****************************/
/*
 The first photon is emitted from the star. If the star is outside the box, we step
 along until it enters the box. We may have some extra just as the photon enters the box
 because I have chosen the step size to be 1 bin.
*/
int FIRST_PHOTON(double *x_new_ptr, double *y_new_ptr, double *z_new_ptr,
                  double *delta_x_ptr, double *delta_y_ptr, double *delta_z_ptr,
                  struct INP_PAR inp_par, struct STARS star, double *angle,
                 dsfmt_t dsfmt, double *ran_array, int nrandom, int *ran_ctr, uint32_t *init_seed)
{
    double  theta, phi;
    double  xp, yp, zp, x_new = 0, y_new = 0, z_new = 0;
    double  dp, d_new, d_old, dtheta, dphi;
    double  delta_x = 0, delta_y = 0, delta_z = 0;
    int     NEW_PHOTON=TRUE;
    int     iter = 0;
    double  step_delta;
    
    //Random number
    double   unif_rand;
    
    while ((NEW_PHOTON == TRUE) && (iter <= inp_par.num_photon)) {
        ++iter;
        /*
         The x, y, and z position of the Sun are tagged to the cell of the dust array.
         The x, y, and z position of the star are in parsecs based on the coordinates of the star.
         xp, yp, zp are the indices of the star in the dust array. As dust_binsize is the size of each
         bin in parsecs, we convert the coordinates by dividing by the binsize.
         dp is the distance in bins. It has to be multiplied by the binsize squared to convert into
         parsecs squared.
        */
        xp = star.x/inp_par.dust_binsize + inp_par.sun_x;
        yp = star.y/inp_par.dust_binsize + inp_par.sun_y;
        zp = star.z/inp_par.dust_binsize + inp_par.sun_z;
        x_new = xp;
        y_new = yp;
        z_new = zp;
        dp = CALC_DIST(inp_par.sun_x, inp_par.sun_y, inp_par.sun_z, xp, yp, zp);
        d_new = dp;
        d_old = dp;

        unif_rand =     GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, ran_ctr, init_seed);
        theta = CALC_THETA(angle, unif_rand);
        phi = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, ran_ctr, init_seed)*TWOPI; //Phi is between 0 and 360 degrees
        
        /*
         The deltas are calculated from the standard coordinate system where North is along the z axis,
         except that the angle is measured from the z axis unlike the declination.
         */
        CALC_DELTA_X(&delta_x, &delta_y, &delta_z, theta, phi);//Calculate the direction of motion
        
        /* If the photon is outside the defined volume, we continue along its path until it
         enters the volume if it is going in the right direction. If it goes in the wrong direction,
         we skip to another photon. This is wasteful of photons but hopefully does not take much time.
         */
 /*
I know the relevant angles so I can use them to get rid of photons I don't like. I know those photons won't
reach so I return -1.
*/
        dtheta = RADDEG(theta);
        dphi   = RADDEG(phi);
        if ((dtheta < star.min_theta) || (dtheta > star.max_theta) || (dphi < star.min_phi) || (dphi > star.max_phi))
            return(-1);
       
        if (CHECK_LIMITS(x_new, y_new, z_new, inp_par) == FALSE){
            step_delta = 1000;
            while (step_delta >= 1) {
                while ((CHECK_LIMITS(x_new, y_new, z_new, inp_par) == FALSE)
                       && (d_new <= d_old))
                {
                    delta_x *= step_delta;
                    delta_y *= step_delta;
                    delta_z *= step_delta;
                    INCREMENT(&x_new, &y_new, &z_new, delta_x, delta_y, delta_z);
                    d_old = d_new;
                    d_new = CALC_DIST(x_new, y_new, z_new, inp_par.sun_x, inp_par.sun_y, inp_par.sun_z);
                }
                INCREMENT(&x_new, &y_new, &z_new, -delta_x, -delta_y, -delta_z);
                dp = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
                delta_x = delta_x/dp;
                delta_y = delta_y/dp;
                delta_z = delta_z/dp;
                d_new = CALC_DIST(x_new, y_new, z_new, inp_par.sun_x, inp_par.sun_y, inp_par.sun_z);
                d_old = d_new;
                step_delta = step_delta/2;;
            }
            INCREMENT(&x_new, &y_new, &z_new, delta_x, delta_y, delta_z);
        }
        //If the photon is not going in the right direction, we have to pick a new photon.
        //I already know that either d_ne > d_old or CHECK_LIMITS is FALSE
        if (CHECK_LIMITS(x_new, y_new, z_new, inp_par) == FALSE)
            iter = -1;
        NEW_PHOTON = FALSE; //Because I'm using all the stars, I don't want to iterate.
    }
    *x_new_ptr = x_new;
    *y_new_ptr = y_new;
    *z_new_ptr = z_new;
    *delta_x_ptr = delta_x;
    *delta_y_ptr = delta_y;
    *delta_z_ptr = delta_z;
    return (iter);
}

float MATCH_TAU(double *x_new_ptr, double *y_new_ptr, double *z_new_ptr, struct INP_PAR inp_par,
                    double delta_x_ptr, double delta_y_ptr, double delta_z_ptr,
                    float sigma, float *dust_dist, float tau)
{
    
    double x_new, y_new, z_new;
    double delta_x, delta_y, delta_z;
    float cum_tau, step_size;
    long dust_index=0;
    int cont;
    
    
    x_new = *x_new_ptr;
    y_new = *y_new_ptr;
    z_new = *z_new_ptr;
    delta_x = delta_x_ptr;
    delta_y = delta_y_ptr;
    delta_z = delta_z_ptr;
    cum_tau = 0;
    step_size = 1.;//We travel around in fractions of a bin
    cont = 1;
    /*
     We keep going on in the same direction incrementing tau until we reach our limiting tau.
     We start by using a step size of 1 bin and then use a finer grid to be more precise.
     Note conversion factor in converting dust column into tau.
     */
    while (cont == 1) {
        while ((cum_tau <= tau) &&
               (CHECK_LIMITS(x_new + delta_x, y_new + delta_y, z_new + delta_z, inp_par) == TRUE)){
            INCREMENT(&x_new, &y_new, &z_new, delta_x, delta_y, delta_z);
            dust_index = GET_DUST_INDEX(x_new, y_new, z_new, inp_par);
            cum_tau += dust_dist[dust_index] * sigma * step_size * inp_par.dust_binsize;
        }
        /*
         There are two ways we can leave the previous loop. Either we reach the edge or we reach the maximum tau.
         if we leave because we are out of the array, we can just skip out to the next step.
         If we are still in the box, we decrease the step size to increase our resolution.
         */
        if (step_size < 0.1)
            cont = 0;
        if ((cont == 1) && (CHECK_LIMITS(x_new, y_new, z_new, inp_par) == TRUE)){
            cum_tau -= dust_dist[dust_index] * sigma * step_size *inp_par.dust_binsize;
            INCREMENT(&x_new, &y_new, &z_new, -delta_x, -delta_y, -delta_z);
            delta_x = delta_x/10.;
            delta_y = delta_y/10.;
            delta_z = delta_z/10.;
        }
        step_size = step_size/10.;
    }
    *x_new_ptr = x_new;
    *y_new_ptr = y_new;
    *z_new_ptr = z_new;
    
    return (cum_tau);
}

double GET_RANDOM_ARRAY(dsfmt_t dsfmt, double *ran_array, int nrandom, int *ran_ctr, uint32_t *init_seed)
{
    double ran_number;
    uint32_t seed;
    
    if (*ran_ctr < nrandom){
        ran_number = ran_array[*ran_ctr];
        *ran_ctr = *ran_ctr + 1;
    } else {
        dsfmt_init_gen_rand(&dsfmt, *init_seed);
        seed =  dsfmt_genrand_uint32(&dsfmt);
        dsfmt_fill_array_close_open(&dsfmt, ran_array, nrandom);
        *ran_ctr = 0;
        ran_number = ran_array[*ran_ctr];
        *ran_ctr = *ran_ctr + 1;
        *init_seed = seed;
//        printf("%u %lf\n", seed, ran_number);
    }
        
    return (ran_number);
}
