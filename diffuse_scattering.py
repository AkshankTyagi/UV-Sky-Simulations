from astropy.io import fits
import numpy as np
import ast
import time
from configparser import ConfigParser
from view_orbit import get_folder_loc, get_cords_from_ra_dec
from star_spectrum import *
from star_data import *
from dust_scatter import *
from read_dust_files import *
import csv

# albedo = 0.36
# g = 0.5

NHIP_STARS = 118218  #Number of stars in Hipparcos catalog
nrandom = 100
NA = 100000  # NANGLE

folder_loc = get_folder_loc()
params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global dust_file, dust_col_file, sigma_file, hipp_file, castelli_file
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    dust_file = config.get(param_set, 'dust_file')
    dust_col_file = config.get(param_set, 'dust_col_file')
    sigma_file = config.get(param_set, 'sigma_file')


read_parameter_file()
dust_par = read_input_parameters()
# global dust_file, dust_col_file, sigma_file, hipp_file, castelli_file
VERSION = "May 2024"
dust_par.version = VERSION  # Assuming VERSION is a predefined constant

# dust_dist = dust_read(dust_par, dust_file)  # Read dust file
# print('dust_dist\n',dust_dist[0])
# col_density = dust_read(dust_par, dust_col_file)  # File containing column densities
# print('col_density\n',col_density[0])
# sigma = read_cross_sec(sigma_file, dust_par)
# print(f"sigma \n{sigma}")  

# # Probabilities for scattering
# angle = np.zeros(NA)
# fscat = np.zeros(NA)
# SETUP_ANGLE(angle)
# SETUP_SCATTER(fscat, dust_par.g)
# print('angle\n',angle)
# print('fscat\n',fscat)


# # Initialize variables
# x_new = y_new = z_new = theta = phi = xp = yp = zp = dust_index = tau = nphoton = 0
# delta_x = delta_y = delta_z = 0

# # Set the seed for random number generation
# time1 = 2**30
# time2 = int(time.time())
# init_seed = time2 % time1
# print(time2, time1, init_seed)
# # init_seed = np.random.seed(init_seed)
# # print(init_seed)

# # Initialize random array using numpy
# ran_array = np.random.rand(nrandom)
# print(ran_array)
# ran_ctr = nrandom
# print(ran_ctr)
# # The dust and energy arrays are continuously built up
# dust_arr = np.zeros((3600, 1800))

# Open and Read Stellar Files

STARS = read_hipparcos_data(hipp_file)
all_spectra = READ_CASTELLI_SPECTRA(castelli_file) #, mag_threshold=False
hipstars = []
star_list = []
print(STARS[0:3])
print('working1')
for j in range(len(STARS['hip'])):
    star = stars()
    star.hip_no = STARS['hip'].iloc[j]
    star.mag = STARS['mag'].iloc[j]
    star.sp_type = STARS['Spectral_type'].iloc[j]
    star.parallax = STARS['trig_parallax'].iloc[j]
    star.ra = STARS['ra_deg'].iloc[j]
    star.dec = STARS['de_deg'].iloc[j]
    star.ebv = STARS['B-V'].iloc[j]
    star.temperature= GET_STAR_TEMP(star.sp_type)
    star.distance= GET_DISTANCE(star.parallax)
    star.x, star.y, star.z = get_cords_from_ra_dec(star.ra, star.dec, star.distance )
    data = all_spectra[star.temperature]
    star.wavelengths = data['wavelength']
    star.Iflux = data['spectrum']
    star.scale, star.tot_photons = GET_SCALE_FACTOR(star.mag, star.parallax, star.ebv, star.wavelengths, star.Iflux)
    # print(star.hip_no, list(zip(star.wavelengths,star.tot_photons)))
    # star.min_phi = 0
    # star.max_phi = 0
    # star.min_theta = 0
    # star.max_theta = 0
    star_list.append([j]*len(star.wavelengths))
    hipstars.append(star)
del STARS
NSTARS = len(hipstars)
# print(hipstars)


# Arrays related to the stars
starlog = np.zeros(NSTARS, dtype=int)
distlog = np.zeros(NSTARS, dtype=float)
scatlog = np.zeros(NSTARS, dtype=float)
totlog = np.zeros(NSTARS, dtype=float)
misslog = np.zeros(NSTARS, dtype=int)
exclude_stars = np.zeros(NSTARS, dtype=int)

# Calculate coordinates of the dust distribution edges
x1 = [-dust_par.sun_x * dust_par.dust_binsize, -dust_par.sun_x*dust_par.dust_binsize,
      (dust_par.dust_xsize -dust_par.sun_x)*dust_par.dust_binsize, (dust_par.dust_xsize - dust_par.sun_x)*dust_par.dust_binsize,
      -dust_par.sun_x*dust_par.dust_binsize, -dust_par.sun_x*dust_par.dust_binsize,
      (dust_par.dust_xsize -dust_par.sun_x)*dust_par.dust_binsize, (dust_par.dust_xsize - dust_par.sun_x)*dust_par.dust_binsize]
y1 = [-dust_par.sun_y * dust_par.dust_binsize, (dust_par.dust_ysize - dust_par.sun_y)*dust_par.dust_binsize,
      -dust_par.sun_y*dust_par.dust_binsize, (dust_par.dust_ysize - dust_par.sun_y)*dust_par.dust_binsize,
      -dust_par.sun_y*dust_par.dust_binsize, (dust_par.dust_ysize - dust_par.sun_y)*dust_par.dust_binsize,
      -dust_par.sun_y*dust_par.dust_binsize, (dust_par.dust_ysize - dust_par.sun_y)*dust_par.dust_binsize]
z1 = [-dust_par.sun_z * dust_par.dust_binsize, -dust_par.sun_z*dust_par.dust_binsize,
      -dust_par.sun_z*dust_par.dust_binsize, -dust_par.sun_z*dust_par.dust_binsize,
      (dust_par.dust_zsize - dust_par.sun_z)*dust_par.dust_binsize, (dust_par.dust_zsize - dust_par.sun_z)*dust_par.dust_binsize,
      (dust_par.dust_zsize - dust_par.sun_z)*dust_par.dust_binsize, (dust_par.dust_zsize - dust_par.sun_z)*dust_par.dust_binsize]


print (x1, y1, z1)
print('working2')

# Initialize variables for min and max coordinates
min_x = min(x1)
max_x = max(x1)
min_y = min(y1)
max_y = max(y1)
min_z = min(z1)
max_z = max(z1)

import math


# Calculate theta angle
for star in hipstars:
    # Initialize hipstars attributes
    star.min_phi = 360
    star.max_phi = 0
    star.min_theta = 180
    star.max_theta = 0
    for j in range(8):
        theta = (z1[j] -  star.z) / math.sqrt((x1[j] - star.x)**2 + (y1[j] - star.y)**2 + (z1[j] - star.z)**2)
        theta = 90 - math.degrees(math.asin(theta))
        if theta > star.max_theta :
            star.max_theta  = theta
        if theta < star.min_theta:
            star.min_theta = theta

    # Adjust theta angle based on position
    if min_x < star.x < max_x and min_y < star.y < max_y:
        if star.z > min_z:
            star.max_theta = 180
        if star.z < max_z:
            star.min_theta = 0

# Calculate phi angle
    phibox = [0] * 8  # Initialize the phibox list

    for j in range(8):
        phibox[j] = (x1[j] - star.x)
        ftmp = math.sqrt(pow(x1[j] - star.x, 2) + pow(y1[j] - star.y, 2))
        phibox[j] /= ftmp
        phibox[j] = math.degrees(math.acos(phibox[j]))

    if (star.x < min_x) and (star.x < max_x) and (star.y < min_y) and (star.y < max_y):
        star.min_phi = phibox[2]
        star.max_phi = phibox[1]
    if (star.x < min_x) and (star.x < max_x) and (star.y > min_y) and (star.y < max_y):
        star.min_phi = 0
        star.max_phi = 360
    if (star.x < min_x) and (star.x < max_x) and (star.y > min_y) and (star.y > max_y):
        star.min_phi = 360 - phibox[0]
        star.max_phi = 360 - phibox[3]
    if (star.x > min_x) and (star.x < max_x) and (star.y < min_y) and (star.y < max_y):
        star.min_phi = phibox[2]
        star.max_phi = phibox[0]
    if (star.x > min_x) and (star.x < max_x) and (star.y > min_y) and (star.y < max_y):
        star.min_phi = 0
        star.max_phi = 360
    if (star.x > min_x) and (star.x < max_x) and (star.y > min_y) and (star.y > max_y):
        star.min_phi = 360 - phibox[1]
        star.max_phi = 360 - phibox[3]
    if (star.x > min_x) and (star.x > max_x) and (star.y < min_y) and (star.y < max_y):
        star.min_phi = phibox[3]
        star.max_phi = phibox[0]
    if (star.x > min_x) and (star.x > max_x) and (star.y > min_y) and (star.y < max_y):
        star.min_phi = phibox[3]
        star.max_phi = 360 - phibox[2]
    if (star.x > min_x) and (star.x > max_x) and (star.y > min_y) and (star.y > max_y):
        star.min_phi = 360 - phibox[1]
        star.max_phi = 360. - phibox[2]


# Assuming NSTARS, hipstars, and star_list are defined somewhere in your code
print('working3')
# Define a function to sort stars based on tot_photons
def sort_stars_by_photons(star_list, hipstars ):
    print(f"Sorting stars {NSTARS}")
    for w in range(len(hipstars[0].wavelengths)):
        for i in range(NSTARS):
            for j in range(i + 1, NSTARS):
                # print (star_list[j][w], hipstars[star_list[j][w]].tot_photons[w], hipstars[j].tot_photons[w] )
                if hipstars[star_list[j][w]].tot_photons[w] > hipstars[star_list[i][w]].tot_photons[w]:
                    star_list[i][w], star_list[j][w] = star_list[j][w], star_list[i][w]

# Open the file "ordered_list_of_stars.list" for reading or writing
try:
    ordered_list = open(f'{folder_loc}ordered_list_of_stars_{NSTARS}.list', "r")
    print("READING STAR LIST FROM FILE")
    for i in range(NSTARS):
        line = ordered_list.readline()
        ltmp = ast.literal_eval(line.strip())
        star_list[i] = ltmp
    ordered_list.close()
except FileNotFoundError:
    # Sort stars if the file is not found
    sort_stars_by_photons(star_list, hipstars)
    # Write sorted star list to the file
    ordered_list = open(f'{folder_loc}ordered_list_of_stars_{NSTARS}.list', "w")
    for i in range(NSTARS):
        ordered_list.write(f"{star_list[i]}\n")
    ordered_list.close()
# print(star_list[0:2])

# Calculate weighted probabilities for selecting stars
def calc_weigthed_probab():
    star_wgt = [[0]*len(hipstars[0].wavelengths)]* NSTARS
    # print(star_wgt)
    for w in range(len(hipstars[0].wavelengths)):
        star_wgt[0][w] = hipstars[star_list[0][w]].tot_photons[w]
        for i in range(1, NSTARS):
            if hipstars[star_list[i][w]].tot_photons[w] < 0:
                print(f'{hipstars[star_list[i][w]].hip}--- photons @ {hipstars[star_list[i][w]].wavelengths[w]}:  {hipstars[star_list[0][w]].tot_photons[w]} ')
            star_wgt[i][w] = star_wgt[i-1][w] + hipstars[star_list[i][w]].tot_photons[w]
    return star_wgt

star_wgt = calc_weigthed_probab()
# print(star_wgt[0])
print('working4')
def CHECKPOINT(dust_arr, inp_par, nphoton, tot_star, wcs, hipstars, starlog, misslog, totlog, distlog, scatlog):
    # Time related
    current_time = time()
    print(f"Checkpoint of {nphoton} photons at {ctime(current_time)}")

    # Add to cumulative grids. Scale by the number of photons from the star over
    # the number of photons in the simulation. Write them out to FITS files.

    write_fits_file(wcs, dust_arr, nphoton, tot_star, inp_par)

    with open("datalogger.txt", "w") as logfile:
        logfile.write("HIP_NO Dist.  Star_flux Star_phot Miss_phot Dist_sca scat_flux tot_flux\n")
        for i in range(len(hipstars)):
            logfile.write(f"{hipstars[i]['HIP_NO']} {hipstars[i]['distance']} {hipstars[i]['tot_photons']} "
                          f"{starlog[i]} {misslog[i]} {distlog[i]} {scatlog[i]} {totlog[i]}\n")

def write_fits_file(wcs_out, grid, nphoton, tot_star, inp_par):
    filename = "scattered.fits"
    axes = [wcs_out.nx, wcs_out.ny]
    dust_out = np.zeros((wcs_out.ny, wcs_out.nx), dtype=np.float32)

    for i in range(wcs_out.nx):
        for j in range(wcs_out.ny):
            ipixel = i + j * wcs_out.nx
            dust_out[j, i] = grid[ipixel] * tot_star / nphoton

    hdu = fits.PrimaryHDU(dust_out)
    hdu.header["CRVAL1"] = wcs_out.xrefval
    hdu.header["CRPIX1"] = wcs_out.xrefpix
    hdu.header["CDELT1"] = wcs_out.xinc
    hdu.header["CROTA1"] = wcs_out.rot
    hdu.header["CTYPE1"] = f"GLON{wcs_out.coordtype}"
    hdu.header["CRVAL2"] = wcs_out.yrefval
    hdu.header["CRPIX2"] = wcs_out.yrefpix
    hdu.header["CDELT2"] = wcs_out.yinc
    hdu.header["CROTA2"] = wcs_out.rot
    hdu.header["CTYPE2"] = f"GLAT{wcs_out.coordtype}"
    hdu.header["DATAMIN"] = np.min(dust_out)
    hdu.header["DATAMAX"] = np.max(dust_out)
    hdu.header["NPHOT"] = nphoton
    hdu.header["ALBEDO"] = inp_par.albedo
    hdu.header["G"] = inp_par.g
    hdu.header["WAVELENG"] = inp_par.wave
    hdu.header["COMMENT"] = f"Dust file: {inp_par.dust_file[:30]}"
    hdu.header["COMMENT"] = inp_par.version

    hdu.writeto(filename, overwrite=True)

############

nphoton = 0
time1 = time.time()
print("Beginning Scattering")
phot_log_file = open("every_photon.log", "w")
fix_rnd = 0

while nphoton < dust_par.num_photon:
    if fix_rnd == 1:
        print("***********Warning: Fixed Random Numbers*************")
        init_seed = 1645310043
        fix_rnd = -1

    if (nphoton % 10000000 == 0) and (nphoton > 0):
        phot_log_file.close()
        if dust_par.print_debug == "yes":
            phot_log_file = open("every_photon.log", "a")
        else:
            phot_log_file = open("every_photon.log", "w")
        CHECKPOINT(dust_arr, dust_par, nphoton, tot_star, wcs, hipstars, starlog, misslog, totlog, distlog, scatlog)
        time2 = time.time()
        print("Time taken for loop:", time2 - time1)
        time1 = time2

    NEW_PHOTON = True

    random_val = GET_RANDOM_ARRAY( ran_array, nrandom, ran_ctr, init_seed) * tot_star
    istar = 0
    while star_wgt[istar] < random_val and istar < NSTARS:
        istar += 1

    istar = star_list[istar]
    star = hipstars[istar]
    starlog[istar] += 1

    tst = FIRST_PHOTON(x_new, y_new, z_new, delta_x, delta_y, delta_z, dust_par, star, angle, dsfmt, ran_array, nrandom, ran_ctr, init_seed)

    nphoton += 1
    xp = x_new
    yp = y_new
    zp = z_new
    max_tau = 0
    intens = 1
    nscatter = 0
    cum_flux = 0

    if tst == -1:
        intens = 0
        nscatter = -1
        misslog[istar] += 1

    while intens >= MIN_INTENS and nscatter <= dust_par.nscatter:
        x_new = xp
        y_new = yp
        z_new = zp

        random_val = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, ran_ctr, init_seed)
        if random_val == 1:
            tau = 1e-10
        elif random_val == 0:
            tau = 100000
        else:
            tau = -log(random_val)

        cum_tau = MATCH_TAU(x_new, y_new, z_new, dust_par, delta_x, delta_y, delta_z, sigma, dust_dist, tau)

        if CHECK_LIMITS(x_new, y_new, z_new, dust_par) and cum_tau >= tau:
            dust_index = GET_DUST_INDEX(x_new, y_new, z_new, dust_par)
            intens *= dust_par.albedo
            flux = DETECT(x_new, y_new, z_new, delta_x, delta_y, delta_z, dust_par)
            dst = CALC_DIST(x_new, y_new, z_new, dust_par.sun_x, dust_par.sun_y, dust_par.sun_z)
            dst *= dust_par.dust_binsize**2
            flux *= intens / dst
            intens -= flux
            extinct = col_density[dust_index] * sigma
            flux *= exp(-extinct)
            cum_flux += flux
            totlog[istar] += flux
            ltmp = int(sqrt(dst))
            distlog[ltmp] += flux
            scatlog[nscatter] += flux

            ra, dec = CART_TO_CELEST(x_new, y_new, z_new, dust_par)
            WRITE_TO_GRID(wcs, ra, dec, flux, dust_arr)

            if dust_par.min_gl_debug < ra < dust_par.max_gl_debug and dust_par.min_gb_debug < dec < dust_par.max_gb_debug:
                phot_log_file.write(f"{star.HIP_NO} {ra} {dec} {xp} {yp} {zp} {x_new} {y_new} {z_new} {extinct} {nscatter} {flux}\n")

            random_val = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, ran_ctr, init_seed)
            theta = CALC_THETA(fscat, random_val)
            phi = GET_RANDOM_ARRAY(dsfmt, ran_array, nrandom, ran_ctr, init_seed) * TWOPI
            xp = x_new
            yp = y_new
            zp = z_new
            SCATTER(delta_x, delta_y, delta_z, theta, phi)
            max_tau = 0
            nscatter += 1

        else:
            intens = 0

time2 = time.time()
print("Time taken for loop:", time2 - time1)
CHECKPOINT(dust_arr, dust_par, nphoton, tot_star, wcs, hipstars, starlog, misslog, totlog, distlog, scatlog)

# Free variables
for i in range(N_CASTELLI_MODELS):
    stellar_spectra[i].spectrum = None
    stellar_spectra[i].wavelength = None

dust_dist = None
dust_arr = None
col_density = None
hipstars = None
star_random = None
star_list = None
misslog = None
starlog = None
totlog = None
phot_log_file.close()


