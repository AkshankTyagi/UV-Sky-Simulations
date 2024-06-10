# Create and save photons file + sorted file for all stars and save it

from configparser import ConfigParser
from view_orbit import get_folder_loc
from star_spectrum import *
from star_data import *
from read_dust_files import *
import pandas as pd
import csv

folder_loc = get_folder_loc()
params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global hipp_file, castelli_file
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    min_lim = float(config.get(param_set, 'limit_min'))
    max_lim = float(config.get(param_set, 'limit_max'))
    return min_lim, max_lim

wave_min, wave_max = read_parameter_file()

def Get_hipstars():
    STARS = read_hipparcos_data(hipp_file)
    all_spectra = READ_CASTELLI_SPECTRA(castelli_file)
    hipstars = []
    print('read files')
    # star_list = []
    for j in range(len(STARS['hip'])):
        star = stars()
        star.hip_no = STARS['hip'].iloc[j]
        print(star.hip_no)
        star.mag = STARS['mag'].iloc[j]
        star.sp_type = STARS['Spectral_type'].iloc[j]
        star.parallax = STARS['trig_parallax'].iloc[j]
        star.ebv = STARS['B-V'].iloc[j]
        star.temperature= GET_STAR_TEMP(star.sp_type, star.hip_no)
        data = all_spectra[star.temperature]
        wavelengths = data['wavelength']
        Iflux = data['spectrum']
        star.scale, photons = GET_SCALE_FACTOR(star.mag, star.parallax, star.ebv, wavelengths, Iflux )
        star.wavelengths, star.Iflux, star.tot_photons = Trim_Spectral_data(data, photons)
        # print(star.scale, list(zip(star.wavelengths, star.Iflux, star.tot_photons)))
        # print(star.hip_no, list(zip(star.wavelengths,star.tot_photons)))
        # star_list.append([j]*len(star.wavelengths))
        hipstars.append(star)
    del STARS, data
    print("hipstars saved")

def Create_Allstars_flux(hipstars):
    data_dict = {'Wavelengths': hipstars[0].wavelengths}
    # Populate the dictionary with HIP numbers as keys and total photons as values
    for star in hipstars:
        data_dict[star.hip_no] = star.tot_photons
    # Convert data to ensure it is little-endian
    for key in data_dict:
        if isinstance(data_dict[key], np.ndarray) and data_dict[key].dtype.byteorder == '>':
            data_dict[key] = data_dict[key].byteswap().newbyteorder()

    # Create a DataFrame from the dictionary and save it as csv file
    df = pd.DataFrame(data_dict)
    df.to_csv(f'diffused_data\\Allstars_flux_data[{wave_min},{wave_max}].csv', index=False)
    print('csv created')


Create_Allstars_flux(hipstars)

