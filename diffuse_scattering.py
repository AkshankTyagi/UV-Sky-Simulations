from astropy.io import fits
import numpy as np
import ast
import time
import math
from astropy.wcs import WCS
import cProfile
import pstats
import io
import multiprocessing as mp
from configparser import ConfigParser


from Params_configparser import get_folder_loc
from Coordinates import *
from star_spectrum import *
from star_data import *
from dust_scatter import *
from read_dust_files import *
from create_csv_files import calc_weigthed_probab, Get_hipstars, CHECKPOINT
import csv

# albedo = 0.36
# g = 0.5

# NHIP_STARS = 118218  #Number of stars in Hipparcos catalog
nrandom = 100
NANGLE = 100000  # NANGLE
MIN_INTENS =1e-20

folder_loc , params_file = get_folder_loc()
# params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global dust_file, dust_col_file, sigma_file, hipp_file, castelli_file
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    dust_file = config.get(param_set, 'dust_file')
    dust_col_file = config.get(param_set, 'dust_col_file')
    sigma_file = config.get(param_set, 'sigma_file')
    min_lim = float(config.get(param_set, 'limit_min'))
    max_lim = float(config.get(param_set, 'limit_max'))

    return min_lim, max_lim


print('working1------ Reading dust_files and parameter files ')

wave_min, wave_max = read_parameter_file()
dust_par = read_scatter_parameters()


VERSION = "June 2024"
dust_par.version = VERSION  

dust_dist = dust_read(dust_par, dust_file)  # Read dust file
# print('dust_dist\n',dust_dist[0,0,2])
col_density = dust_read(dust_par, dust_col_file)  # File containing column densities
# print('col_density\n',col_density[0])
sigma = read_cross_sec(sigma_file, dust_par)

# wcs data
wcs_param = read_wcs_parameters()
wcs = WCS(naxis=2)
wcs.wcs.crval = [wcs_param[0], wcs_param[1]]
wcs.wcs.crpix = [wcs_param[2], wcs_param[3]]
wcs.wcs.cdelt = [wcs_param[4], wcs_param[5]]
wcs.wcs.crota = [wcs_param[6], wcs_param[6]]
wcs.wcs.ctype = [f"GLON{wcs_param[7]}", f"GLAT{wcs_param[7]}"]
wcs.array_shape = [ wcs_param[9], wcs_param[8],] 

print(f"sigma \n{sigma}")
print(f"dust_par.dust_xsize:{dust_par.dust_xsize}, dust_par.dust_ysize:{dust_par.dust_ysize}, dust_par.dust_zsize:{dust_par.dust_zsize}\ndust_par.dust_binsize :{dust_par.dust_binsize }")

# Probabilities for scattering
angle = SETUP_ANGLE()
fscat = SETUP_SCATTER( dust_par.g)
# print('angle\n',angle)
# print('fscat\n',fscat)


# Initialize variables
x_new = y_new = z_new = theta = phi = xp = yp = zp = tau = nphoton = 0
delta_x = delta_y = delta_z = 0

# Set the seed for random number generation
# time1 = 2**30
# time2 = int(time.time())
# init_seed = time2 % time1
# print(time2, time1, init_seed)
# init_seed = np.random.seed(init_seed)
# print(init_seed)

# Initialize random array using numpy
# ran_array = np.random.rand(nrandom)
# print(ran_array)
# ran_ctr = nrandom
# print(ran_ctr)

print('\nworking2------ Reading hipstars and initialise log files')


# Read Stellar Files
hipstars = Get_hipstars()
NSTARS = int(len(hipstars))
wavelengths_list = hipstars.loc[0, 'wavelengths']
# wavelengths_str = wavelengths_str.replace('  ', ',')
wavelengths_list = ast.literal_eval(wavelengths_list)
print(f"NSTARS: {NSTARS},\nwavelengths_list:{wavelengths_list}, {wavelengths_list[0]} ")

# np.savetxt(f'excluded_stars_{len(excluded_stars)}.txt', np.array(excluded_stars), fmt='%d')
# print(f'working1:----- Nstars:{NSTARS}, excluded_stars_{len(excluded_stars)}.txt saved, {hipstars[-1].index, hipstars[-1].hip_no}' )


# Arrays related to the stars
# starlog = np.zeros((len(wavelengths_list), NSTARS), dtype=int)   # Counts which star chosen how many times
# distlog = np.zeros((len(wavelengths_list), NSTARS), dtype=float) # flux as a function of distance
# scatlog = np.zeros((len(wavelengths_list), NSTARS), dtype=float) # flux for each scatter
# totlog = np.zeros((len(wavelengths_list), NSTARS), dtype=float)  #keeps count flux from each star
# # print(totlog.shape, totlog[0][2])
# misslog = np.zeros((len(wavelengths_list), NSTARS), dtype=int)
# exclude_stars = np.zeros(NSTARS, dtype=int)

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


# print (x1, y1, z1)
print('\nworking3------ reading sorted csv files')

#sorting and weighting star photons data
weighted_list, sorted_stars = calc_weigthed_probab(wavelengths_list, hipstars)


# Calculate theta angle for each star
def calculate_theta(star_x, star_y, star_z, x1, y1,z1):
    star_min_theta = 180
    star_max_theta = 0

    for j in range(8):
        theta = (z1[j] -  star_z) / math.sqrt((x1[j] - star_x)**2 + (y1[j] - star_y)**2 + (z1[j] - star_z)**2)
        theta = 90 - math.degrees(math.asin(theta))
        if theta > star_max_theta :
            star_max_theta  = theta
        if theta < star_min_theta:
            star_min_theta = theta

    # Initialize variables for min and max coordinates
    min_x = min(x1)
    max_x = max(x1)
    min_y = min(y1)
    max_y = max(y1)
    min_z = min(z1)
    max_z = max(z1)

    # Adjust theta angle based on position
    if min_x < star_x < max_x and min_y < star_y < max_y:
        if star_z > min_z:
            # print('a')
            star_max_theta = 180
        if star_z < max_z:
            # print('b')
            star_min_theta = 0
    # else:
    #     print(f"{min_x, star_x , max_x , min_y , star_y , max_y}")
    # print(f"star_min_theta, star_max_theta:{star_min_theta, star_max_theta}")
    return star_min_theta, star_max_theta

# Calculate phi angle for each star
def calculate_phi(star_x, star_y, star_z, x1, y1, z1):
    phibox = [0] * 8  # Initialize the phibox list
    star_min_phi = 360
    star_max_phi = 0

    for j in range(8):
        phibox[j] = (x1[j] - star_x)
        ftmp = math.sqrt(pow(x1[j] - star_x, 2) + pow(y1[j] - star_y, 2))
        phibox[j] /= ftmp
        phibox[j] = math.degrees(math.acos(phibox[j]))

    # Initialize variables for min and max coordinates
    min_x = min(x1)
    max_x = max(x1)
    min_y = min(y1)
    max_y = max(y1)
    min_z = min(z1)
    max_z = max(z1)

    if (star_x < min_x) and (star_x < max_x) and (star_y < min_y) and (star_y < max_y):
        # print(1)
        star_min_phi = phibox[2]
        star_max_phi = phibox[1]
    if (star_x < min_x) and (star_x < max_x) and (star_y > min_y) and (star_y < max_y):
        # print(2)
        star_min_phi = 0
        star_max_phi = 360
    if (star_x < min_x) and (star_x < max_x) and (star_y > min_y) and (star_y > max_y):
        # print(3)
        star_min_phi = 360 - phibox[0]
        star_max_phi = 360 - phibox[3]
    if (star_x > min_x) and (star_x < max_x) and (star_y < min_y) and (star_y < max_y):
        # print(4)
        star_min_phi = phibox[2]
        star_max_phi = phibox[0]
    if (star_x > min_x) and (star_x < max_x) and (star_y > min_y) and (star_y < max_y):
        # print(5)
        star_min_phi = 0
        star_max_phi = 360
    if (star_x > min_x) and (star_x < max_x) and (star_y > min_y) and (star_y > max_y):
        # print(6)
        star_min_phi = 360 - phibox[1]
        star_max_phi = 360 - phibox[3]
    if (star_x > min_x) and (star_x > max_x) and (star_y < min_y) and (star_y < max_y):
        # print(7)
        star_min_phi = phibox[3]
        star_max_phi = phibox[0]
    if (star_x > min_x) and (star_x > max_x) and (star_y > min_y) and (star_y < max_y):
        # print(8)
        star_min_phi = phibox[3]
        star_max_phi = 360 - phibox[2]
    if (star_x > min_x) and (star_x > max_x) and (star_y > min_y) and (star_y > max_y):
        # print(9)
        star_min_phi = 360 - phibox[1]
        star_max_phi = 360. - phibox[2]
    # print(f"star_min_phi, star_max_phi: {star_min_phi, star_max_phi}")
    return star_min_phi, star_max_phi

if hipstars.loc[0,'max_phi'] == 0:
    print('\nworking4------- Calculate theta & phi angle for each star')
    for index, star in hipstars.iterrows():
        hipstars.at[index, 'min_theta'], hipstars.at[index, 'max_theta'] = calculate_theta(star['x'], star['y'], star['z'], x1, y1,z1)
        hipstars.at[index, 'min_phi'], hipstars.at[index, 'max_phi'] = calculate_phi(star['x'], star['y'], star['z'], x1, y1,z1)
    hipstars.to_csv(f'hipstars_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv', index=False)
else:
    print('\nworking4------- not required to Calculate theta & phi angle again !!!! ')

print(f"min_theta={hipstars.loc[0,'min_theta']}, max_theta={hipstars.loc[0,'max_theta']}, min_phi = {hipstars.loc[0,'min_phi']}, max_phi={hipstars.loc[0,'max_phi']}\n")


################################################################################################

print('Begining scattering')

# Define the scatter function
def scattered_light(data):
    i, w = data
    nphoton = 0
    time1 = time.time()
    # phot_log_file = open("every_photon.log", "w")
    # fix_rnd = 0
    tot_star = weighted_list.iloc[int(i)][-1]
    print(f'---{i+1}--- wavelength={w}, N_stars = {NSTARS}, tot_star:{tot_star}')

    # The dust and energy arrays are continuously built up
    dust_arr = np.zeros(( 1800, 3600))
    weighted_list_i = weighted_list.iloc[i].values
    sorted_stars_i = sorted_stars.iloc[i].values
    # dust_arr2 = dust_arr
    
    while nphoton < dust_par.num_photon:
        # start_loop_time = time.time()
        # if fix_rnd == 1:
            # print("***********Warning: Fixed Random Numbers*************")
            # init_seed = 1645310043
            # fix_rnd = -1

        if (nphoton % 1000000 == 0) and nphoton!= 0: # and (nphoton > 0):
            timez = time.time()
            print(f"wavelength: {w}, nphoton: {nphoton}, time for loop: {timez - time1},")
            if (nphoton % 10000000 == 0):
                CHECKPOINT(dust_arr, dust_par, nphoton, tot_star, wcs_param, hipstars, w)
                # plot_diffused_bg(dust_arr * tot_star / nphoton, w, dust_par.num_photon)
            #, starlog, misslog, totlog, distlog, scatlog)
            # plot_diffused_bg(dust_arr*tot_star / nphoton, 1105)
            # phot_log_file.close()
            # if dust_par.print_debug == "yes":
            #     phot_log_file = open("every_photon.log", "a")
            # else:
            #     phot_log_file = open("every_photon.log", "w")
            # time2 = time.time()
            # print("Time taken for loop:", time2 - time1)
            # time1 = time2


        NEW_PHOTON = True
        random_value = GET_RANDOM_ARRAY() * tot_star
        # istar = 1
        # while weighted_list.iloc[i][istar] < random_val and istar < NSTARS:
        #     istar += 1
        # lstar = int(sorted_stars.iloc[i][istar])
        istar = np.searchsorted(weighted_list_i[1:NSTARS+1], random_value) + 1
        lstar = int(sorted_stars_i[istar])
        star = hipstars.loc[lstar]

        # starlog[i][lstar] += 1
    # print(f"nphoton: {nphoton}, istar:{istar}, lstar:{lstar}, star.hip_no:{star['hip_no']}, random_value used:{random_val}")#,   random_val:{random_val}, check value:{weighted_list.iloc[i][istar]}
        # print(f'star.hip_no:{star.hip_no}, star.sp_type:{star.sp_type}, star.distance:{star.distance}, star.scale:{star.distance}')
 
        x_new, y_new, z_new, delta_x, delta_y, delta_z, tst = FIRST_PHOTON(dust_par, star, angle) #,  ran_array, nrandom, ran_ctr, init_seed)


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
            nphoton += 1
            # misslog[i][lstar] += 1
            continue
        else:
            nphoton += 1
        # print(f"{nphoton}, istar:{istar}, hip:{star['hip_no']}, tst:{tst}")#,   random_val:{random_val}, check value:{weighted_list.iloc[i][istar]}

        # print(f"tst:{tst}, intens:{intens}, nscatter:{nscatter}")
        while intens >= MIN_INTENS and nscatter <= dust_par.nscatter:

            x_new = xp
            y_new = yp
            z_new = zp

            random_val = GET_RANDOM_ARRAY( )#ran_array, nrandom, ran_ctr, init_seed)
            # if random_val == 1:
            #     tau = 1e-10
            # elif random_val == 0:
            #     tau = 100000
            # else:
            tau = -np.log(random_val)

            cum_tau = MATCH_TAU(x_new, y_new, z_new, dust_par, delta_x, delta_y, delta_z, sigma, dust_dist, tau)


            if CHECK_LIMITS(x_new, y_new, z_new, dust_par) and cum_tau >= tau:
                dust_index = GET_DUST_INDEX(x_new, y_new, z_new)
                intens *= dust_par.albedo
                # print(x_new, y_new, z_new, delta_x, delta_y, delta_z)
                flux = DETECT(x_new, y_new, z_new, delta_x, delta_y, delta_z, dust_par)
                # print(f"flux1: {flux}")
                dst = CALC_DIST(x_new, y_new, z_new, dust_par.sun_x, dust_par.sun_y, dust_par.sun_z)
                dst *= dust_par.dust_binsize**2
                flux *= intens / dst
                # print(f"flux2: {flux}, factor:{intens / dst}")
                intens -= flux
                extinct = col_density[dust_index[0], dust_index[1], dust_index[2]] * sigma[0]
                flux *= np.exp(-extinct)
                # print(f"flux3: {flux}, factor:{np.exp(-extinct)}")
                cum_flux += flux 
                # totlog[i][lstar] += flux
                ltmp = int(np.sqrt(dst))
                # distlog[i][ltmp] += flux
                # scatlog[i][nscatter] += flux

                ra, dec = conv_cart_to_eq(x_new, y_new, z_new, dust_par)
                # print(f'{nphoton}) istar:{istar}, rand_val:{random_value:.5f}, checkval:{weighted_list.iloc[i][istar]:.5f}, ra: {ra:.3f}, dec: {dec:.3f}, flux: {flux}')#;  extinct:{extinct:.3f}, dist:{ltmp}')
                # print(f"")#,   random_val:{random_val}, check value:{weighted_list.iloc[i][istar]}
                dust_arr = write_to_grid(wcs, ra, dec, flux, dust_arr)

                # dust_arr2 = write_to_gridgg(wcs, ra, dec, flux, dust_arr2)

                # if dust_par.min_gl_debug < ra < dust_par.max_gl_debug and dust_par.min_gb_debug < dec < dust_par.max_gb_debug:
                #     phot_log_file.write(f"{star.hip_no} {ra} {dec} {xp} {yp} {zp} {x_new} {y_new} {z_new} {extinct} {nscatter} {flux}\n")

                random_val = GET_RANDOM_ARRAY()
                theta = CALC_THETA(fscat, random_val)
                phi = GET_RANDOM_ARRAY() * 2*np.pi
                xp = x_new
                yp = y_new
                zp = z_new
                SCATTER(delta_x, delta_y, delta_z, theta, phi)
                max_tau = 0
                nscatter += 1
            else:
                intens = 0
    time2 = time.time()
    print("Time taken for this wavelength_process:", time2 - time1)
    CHECKPOINT(dust_arr, dust_par, nphoton, tot_star, wcs_param, hipstars, w)#, starlog, misslog, totlog, distlog, scatlog)
    plot_diffused_bg(dust_arr*tot_star/nphoton, w, dust_par.num_photon)
    # plot_diffused_bg(dust_arr2*tot_star / nphoton, 110500)
    # print(dust_arr)

    # # Free variables
    # for i in range(N_CASTELLI_MODELS):
    #     stellar_spectra[i].spectrum = None
    #     stellar_spectra[i].wavelength = None

    # dust_dist = None
    # dust_arr = None
    # col_density = None
    # hipstars = None
    # star_random = None
    # sorted_stars = None
    # misslog = None
    # starlog = None
    # totlog = None 
    # phot_log_file.close()

    # return dust_arr * tot_star / dust_par.num_photon



# main
if __name__ == '__main__':
    scatter_wavelengths = []
    for w in dust_par.wave:
        i = wavelengths_list.index(w)
        scatter_wavelengths.append([i,w])


    NProcessor = 10
    start_time = time.time()

    if len(scatter_wavelengths) < NProcessor:
        NProcessor = len(scatter_wavelengths)

    print(f"Input:{scatter_wavelengths}, Nprocessor:{NProcessor}")

    with mp.Pool(processes = NProcessor) as pool:
        print("working!!!!!!")
    # Use the pool to apply the function to each pair in the input list
        pool.map(scattered_light, scatter_wavelengths )
        print("working2!!!!!!")

    print(f"time taken: {time.time() - start_time}")

    
    # pr = cProfile.Profile()
    
    # # Enable profiling
    # pr.enable()
    
    # # Run the main function
    # main(0,1105)
    
    # # Disable profiling
    # pr.disable()
    
    # # Create a StringIO object to store the profiling results
    # s = io.StringIO()
    
    # # Create a Stats object
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
    # # Strip directories and print the stats
    # # ps.strip_dirs().print_stats('diffuse_scattering.py', 'dust_scatter.py', 'coordinates.py')
    # # ps.strip_dirs().print_stats('dust_scatter.py')
    # # print(s.getvalue())
    # ps.strip_dirs().print_stats('dust_scatter.py')
    # print(s.getvalue())


