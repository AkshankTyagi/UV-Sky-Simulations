
# from star_spectrum import *
# from star_data import *
# from dust_scatter import *
from configparser import ConfigParser
from astropy.io import fits
import numpy as np
from view_orbit import get_folder_loc
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
    # Dust parameters:
    num_photons = float(config.get(param_set, 'No_photons'))
    num_scatter = float(config.get(param_set, 'No_scatter'))
    albedo = float(config.get(param_set, 'Albedo'))
    wavelength = config.get(param_set, 'wavelength').strip('[]')
    wavelength_list = [float(value) for value in wavelength.split(',')]
    # print (wavelength_list[[0,1]])
    phase_func = float(config.get(param_set, 'Phase_func'))

    return num_photons, num_scatter, wavelength_list, albedo, phase_func

def read_input_parameters():
    dust_par = Dust_params()  # Create an instance of the Dust_params class
    
    num_photons, num_scatter, wavelength, albedo, phase_func = read_parameter_file()
    print(num_photons)
    dust_par.num_photon = int(num_photons)
    dust_par.nscatter = int(num_scatter)
    dust_par.wave = list(wavelength)
    dust_par.albedo = float(albedo)
    dust_par.g = float(phase_func)

    # Default debugging parameters
    dust_par.print_debug = "no"
    dust_par.min_gl_debug = 0
    dust_par.max_gl_debug = 360
    dust_par.min_gb_debug = -90
    dust_par.max_gb_debug = 90
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
    # print (sigma_file, hipp_file)
    with open(sigma_file, "r") as fp:
        lines = fp.readlines()
        # print(lines[0])

        # Find the line with the start of the cross-section data
        start_line_index = 0
        for i, line in enumerate(lines):
            if line.startswith("-"):
                start_line_index = i + 1
                print(start_line_index)
                break

        # Read cross-section data
        wavelengths = []
        crossX = []
        for line in lines[start_line_index:]:
            if line.startswith("-"):
                print(line)
                break
            parts = line.split()
            if len(parts) >= 4:
                wave = float(parts[0])
                cross_sec = float(parts[3])
                wavelengths.append(wave * 10000)  # Convert to Ångstroms
                crossX.append(cross_sec)

    # Find the nearest wavelength
    wavel = dust_par.wave
    sigma = []
    for wave in wavel:
        wave = float(wave)
        index = min(range(len(wavelengths)), key=lambda i: abs(wavelengths[i] - wave))
        # Calculate sigma (extinction cross-section per H atom scaled to atoms/cm^3)
        sigma.append(crossX[index]* pc_to_cm(1))
    print (f'wavel:{wavel}')

    return sigma

