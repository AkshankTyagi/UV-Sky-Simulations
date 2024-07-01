
# from star_spectrum import *
# from star_data import *
# from dust_scatter import *
from configparser import ConfigParser
from astropy.io import fits
import numpy as np
from Params_configparser import get_folder_loc
# from diffuse_scattering import pc_to_cm

class Dust_params:
    def __init__(self):
        self.dust_xsize = 0
        self.dust_ysize = 0
        self.dust_zsize = 0
        self.sun_x = 0.0
        self.sun_y = 0.0
        self.sun_z = 0.0
        self.dust_binsize = 0.0
        self.num_photon = 0
        self.nscatter = 0
        self.wave = 0.0
        self.albedo = 0.0
        self.g = 0.0
        self.version = ""
        self.print_debug = ""
        self.min_gl_debug = 0.0
        self.max_gl_debug = 0.0
        self.min_gb_debug = 0.0
        self.max_gb_debug = 0.0

class stars:
    def __init__(self):
        self.hip_no = 0
        self.index =0
        self.mag = 0
        self.sp_type = 0
        self.parallax = 0
        self.ra = 0
        self.dec = 0
        self.ebv = 0

        self.x = 0
        self.y = 0
        self.z = 0
        self.distance = 0
        self.temperature = 0
        self.wavelengths = []
        self.Iflux = []
        # self.nstars = nstars
        self.scale = 0
        self.tot_photons = []
        self.min_theta = 0
        self.max_theta = 0
        self.min_phi = 0
        self.max_phi = 0

folder_loc = get_folder_loc()
params_file = f'{folder_loc}init_parameter.txt'

def pc_to_cm(x):
    return 3.08568024696E+18 * x

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)

    global dust_file, dust_col_file, sigma_file, hipp_file, castelli_file
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    dust_file = config.get(param_set, 'dust_file')
    dust_col_file = config.get(param_set, 'dust_col_file')
    sigma_file = config.get(param_set, 'sigma_file')
    return 

def read_scatter_parameters(filename = params_file, param_set = 'Scatter_params'):
    dust_par = Dust_params()  # Create a Dust_params class
    config = ConfigParser()
    config.read(filename)

    # Dust parameters: num_photons, num_scatter, wavelength_list, albedo, phase_func
    dust_par.num_photon = float(config.get(param_set, 'No_photons'))
    dust_par.nscatter = float(config.get(param_set, 'No_scatter'))
    wavelength = config.get(param_set, 'wavelength').strip('[]')
    dust_par.wave = [float(value) for value in wavelength.split(',')]
    dust_par.albedo = float(config.get(param_set, 'Albedo'))
    dust_par.g = float(config.get(param_set, 'Phase_func'))

    # Default debugging parameters
    dust_par.print_debug = float(config.get(param_set, 'print_debug'))
    dust_par.min_gl_debug = float(config.get(param_set, 'min_gl_debug'))
    dust_par.max_gl_debug = float(config.get(param_set, 'max_gl_debug'))
    dust_par.min_gb_debug = float(config.get(param_set, 'min_gb_debug'))
    dust_par.max_gb_debug = float(config.get(param_set, 'max_gb_debug'))

    return dust_par

def dust_read(dust_par, filename):
    with fits.open(filename) as hdul:
        data = hdul[0].data
        naxes = data.shape
        dust_arr = np.array(data, dtype=float)
        dust_par.dust_xsize = float(naxes[0])
        dust_par.dust_ysize = float(naxes[1])
        dust_par.dust_zsize = float(naxes[2])
        dust_par.dust_binsize = float(hdul[0].header['CRDELT1'])
        dust_par.sun_x = float(hdul[0].header['CRPIX1'])
        dust_par.sun_y = float(hdul[0].header['CRPIX2'])
        dust_par.sun_z = float(hdul[0].header['CRPIX3'])

        print(f"dust_par.dust_xsize:{dust_par.dust_xsize}, dust_par.dust_binsize :{dust_par.dust_binsize }")
        print(f"dust_par.sun_x:{dust_par.sun_x}\ndust_par.sun_y: {dust_par.sun_y}, \ndust_par.sun_z : {dust_par.sun_z}")

    return dust_arr

def read_cross_sec(sigma_file, dust_par):
    wavelengths = []
    crossX = []
    start_line_found = False

    with open(sigma_file, "r") as fp:
        for line in fp:
            if not start_line_found:
                if line.startswith("-"):
                    start_line_found = True
                continue
            if line.startswith("-"):
                break
            parts = line.split()
            if len(parts) >= 4:
                wave = float(parts[0])
                cross_sec = float(parts[3])
                wavelengths.append(wave * 10000)  # Convert to Ångstroms
                crossX.append(cross_sec)
    

    wavelengths = np.array(wavelengths)
    crossX = np.array(crossX)
    wavel = np.array(dust_par.wave, dtype=float)
    sigma = []

    for wave in wavel:
        index = np.argmin(np.abs(wavelengths - wave))
        sigma_value = crossX[index] * pc_to_cm(1)  # Convert to the desired units
        sigma.append(sigma_value)

    print(f'wavel: {wavel}')
    return sigma

