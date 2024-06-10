# Create and save photons file + sorted file for all stars and save it

from configparser import ConfigParser
from view_orbit import get_folder_loc, get_cords_from_ra_dec
from star_spectrum import *
from star_data import *
from read_dust_files import *
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

# Saving The photon Values for all stars to a csv file
photon_data = []
header = ['Wavelengths'] + [f'{hipstars[i].hip_no}' for i in range(len(hipstars))]
photon_data.append(header) 
num_wavelengths = len(hipstars[0].wavelengths)
print('header created')
# print(num_wavelengths,len(hipstars), photon_data)
for i in range(num_wavelengths):
    row = [hipstars[0].wavelengths[i]]  # First column: Wavelength
    row.extend([hipstars[j].tot_photons[i] for j in range(len(hipstars))])  
    photon_data.append(row)
print('photon_data created')

# Write the data to a CSV file
with open(f'allstars_flux_data[{wave_min},{wave_max}].csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(photon_data)

print('csv created')