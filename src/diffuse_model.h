/*******************************************************************
Header file
 *********************************************************************/

#include <math.h>
#include <stdlib.h>
//#include "parse_par.h"
//#include "fitsio.h"
#include <string.h>
#include <time.h> /*Note that this is only needed for informational purposes*/
//#include "randlib.h"
#include "dSFMT.h"


#define PRGM_ID      "Diffuse Light"
#define PI            3.141592653589793238462643
#define TWOPI         6.283185307179586476925286
#define NODATA        -1
#define DEGRAD(x)    ((x)*PI/180.)
#define RADDEG(x)    ((x)*180./PI)
#define NANGLE       100000
#define MIN_INTENS   1e-20
#define PCinAU       206265
#define PCinCM(x)    (3.08568024696E+18 * (x)) /* Converts from Parsec to cm */
#define eRStar       0.0       /* Average Stellar radius (in cm ) */
#define rd           0        /* distance to the reference plane from earth(in PC)*/
#define ERG_TO_PHOT  50306871.92
#define NUMPHOT      5000000
#define WEIGHT       0.05
#define NHIP_STARS      118218  //Number of stars in Hipparcos catalog
static double dsqrarg;
#define SWAP(x, y) ftmp = (x); (x) = (y); (y) = ftmp
#define DSFMT_MEXP 19937    //For the Mersenne Twister algorithm
//For Windows
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ?0 :errno)


/*Common variables*/
#define HIP_DELIM "|"
#define HIP_MAIN_REC_LEN 1000

//Variables
#define MAX_FILE_LENGTH 100
#define N_CASTELLI_SPECTRA 1121
#define N_CASTELLI_MODELS 76
#define N_CROSS_SEC 10000 /*Maximum number of lines in cross-section file*/

struct INP_PAR {
    char par_file[MAX_FILE_LENGTH];
    char dust_file[MAX_FILE_LENGTH];
    char dust_col_density[MAX_FILE_LENGTH];
    char castelli_file[MAX_FILE_LENGTH];
    char sigma_file[MAX_FILE_LENGTH];
    char hipparcos_file[MAX_FILE_LENGTH];
    char wcs_file[MAX_FILE_LENGTH];
    char out_dir[MAX_FILE_LENGTH];
    char star_file[MAX_FILE_LENGTH];
    char extra_stars[MAX_FILE_LENGTH];
    int  hipparcos_star_no;
    long  dust_xsize;
    long  dust_ysize;
    long  dust_zsize;
    double sun_x;
    double sun_y;
    double sun_z;
    double dust_binsize;
    long num_photon;
    long nscatter;
    double wave;
    double albedo;
    double g;
    char   version[FLEN_COMMENT];
    char   print_debug[FLEN_COMMENT];
    double min_gl_debug;
    double max_gl_debug;
    double min_gb_debug;
    double max_gb_debug;
};

/*Definition of the data structure for Hipparcos stars*/
struct STARS
{
    int    HIP_NO;             /* Hipparcos Catalog No             */
    long    HD_NO;              /* HD Number                        */
    char    sp_type[13];        /* Spectral Type.                   */
    float   distance;           /* Distance to the Star (in pc)     */
    float   ra;                 /* Right Ascension (2000)           */
    float   dec;                /* Declination                      */
    double  x;                  /* X Cartesian Coordinate.          */
    double  y;                  /* Y Cartesian Coordinate.          */
    double  z;                  /* Z Cartesian Coordinate.          */
    float   B_mag;              /* B mag of the star                */
    float   V_mag;              /* V mag of the star                */
    float   gl;                 /* Galactic longitude l             */
    float   gb;                 /* Galactic latitude  b             */
    int     temperature;        /* index of temperature of the star */
    int     nstars;             /*Number of elements*/
    float   scale;              /*Scale factor for stellar bightness*/
    float   ebv;                /*E(B-V)*/
    double  tot_photons;        /*Total number of photons*/
    float   min_theta;          /*Minimum angle for scattered photons*/
    float   max_theta;          /*Maximum angle for scattered photons*/
    float   min_phi;            /*Minimum angle for scattered photons*/
    float   max_phi;            /*Maximum angle for scattered photons*/
};

struct SPECTRA
{
    char filename[MAX_FILE_LENGTH];
    float *spectrum;
    float *wavelength;
};

struct WCS_HEAD {
                  double xrefval;
                  double yrefval;
                  double xrefpix;
                  double yrefpix;
                  double xinc;
                  double yinc;
                  double rot;
                  char coordtype[5];
                  long nx;
                  long ny;
                 };

/***************** FUNCTION Definitions ********************************/
//read_files.c
struct  INP_PAR READ_INPUT_PARAMETERS(int argc, char *argv[]);
float   *DUST_READ(struct INP_PAR *inp_par, char filename[MAX_FILE_LENGTH]);
double  READ_CROSS_SEC(struct INP_PAR file_list);
void    WRITE_FITS_FILE(struct WCS_HEAD wcs_out, float *grid, long nphoton,
                     float tot_star, struct INP_PAR inp_par);
void    WRITE_TO_GRID(struct WCS_HEAD wcs, double ra, double dec, double flux,
                   float *grid);
void   CHECKPOINT(float *dust_arr,
                 struct INP_PAR inp_par, long nphoton,
                 float tot_star, struct WCS_HEAD wcs, struct STARS *hipstars,
                 int *starlog, int *misslog, float *totlog, float *distlog,
                  float *scatlog);


//star_programs.c
int HIP_READ_LINE(FILE *hipfile, struct STARS *line);
int READ_CASTELLI_SPECTRA(char *, struct SPECTRA *);
int GET_SCALE_FACTOR(struct STARS *hipstars,
                     struct SPECTRA *stellar_spectra,
                     struct INP_PAR inp_par);
int GET_STAR_TEMP(struct STARS *hipline);
int  EXTRA_READ_LINE(FILE *extra_star_file, struct STARS *line);


//scatter_math.c
double  PHASE_FUNCTION(double cangle, double g);
void    SETUP_ANGLE(double *angle);
void    SETUP_SCATTER(double *fscat, double g);
double  CALC_THETA(double angle[NANGLE], float random);
double  DETECT(double xstart, double ystart, double zstart,
              double xend,   double yend,   double zend, struct INP_PAR dust);
void    SCATTER(double *delta_x, double *delta_y, double *delta_z, double theta, double phi);
int     CHECK_LIMITS(double x, double y, double z, struct INP_PAR dust);
void    INCREMENT(double *x_new,  double *y_new,  double *z_new,
               double delta_x, double delta_y, double delta_z);
void    CALC_DELTA_X(double *delta_x, double *delta_y, double *delta_z,
                  double theta, double phi);
long    GET_DUST_INDEX(double x_new, double y_new, double z_new, struct INP_PAR dust);
int    FIRST_PHOTON(double *x_new, double *y_new, double *z_new,
                  double *delta_x, double *delta_y, double *delta_z,
                  struct INP_PAR inp_par, struct STARS star, double *angle,
                    dsfmt_t dsfmt, double *ran_array, int nrandom, int *ran_ctr, uint32_t *init_seed);
double  EXTINCT(double xnew, double ynew, double znew, float *dust_arr,
               struct INP_PAR inp_par, double sigma);
float MATCH_TAU(double *x_new, double *y_new, double *z_new, struct INP_PAR inp_par,
                double delta_x, double delta_y, double delta_z,
                float sigma, float *dust_dist, float tau);
double GET_RANDOM_ARRAY(dsfmt_t dsfmt, double *ran_array, int nrandom, int *ran_ctr, uint32_t *init_seed);


//coordinates
int     CONV_EQ_TO_GAL(double ra, double dec, double *gl, double *gb);
void    CART_TO_CELEST(double x, double y, double z, double *ra, double *dec,
                       struct INP_PAR inp_par);
double  CALC_DIST(double a, double b, double c, double x, double y, double z);
struct  WCS_HEAD WCS_READ(char wcs_file[MAX_FILE_LENGTH]);
void WRITE_ENERGY_FILE(float *energy_arr, struct INP_PAR inp_par,
                       char filename[FLEN_FILENAME], long nphoton, float tot_star);
