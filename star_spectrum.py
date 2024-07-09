from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from configparser import ConfigParser
import math
import pandas as pd
from Params_configparser import get_folder_loc

folder_loc = get_folder_loc()

class StellarSpectrum:
    def __init__(self, v1 =0, v2= []):
        self.temperature = v1
        self.spectrum = v2
        # self.wavelength = v3
        # self.photons = v4
        # self.scale = v5


class Spectral_FOV:
    def __init__(self):
        self.frame = []
        self.wavelength = []
        self.spectra_per_star= []
        self.ra = []
        self.dec = []
        self.scale = []
        self.photons = []
        self.frame_size = []

ERG_TO_PHOT = 50306871.92
N_CASTELLI_MODELS = 76

# Castelli_data--
# spec_dir= r"C:\Users\Akshank Tyagi\Documents\GitHub\spg-iiap-UV-Sky-Simulation\Castelli\ckp00"
# params_file = r'C:\Users\Akshank Tyagi\Documents\GitHub\spg-iiap-UV-Sky-Simulation\init_parameter.txt'
params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global Castelli_data 
    Castelli_data = config.get(param_set, 'Castelli_data')
    # Wavelength range:
    min_lim = float(config.get(param_set, 'limit_min'))
    max_lim = float(config.get(param_set, 'limit_max'))
    return min_lim, max_lim

_, _ = read_parameter_file()

def index_greater_than(lst, value):
    for i, elem in enumerate(lst):
        if elem >= value:
            return i
    return None  

def GET_STAR_TEMP(sptype, hip_no = 0):
    sptype = str(sptype)
    temperature = 0 

    if sptype.startswith("sd"):
        sptype = sptype[2:]
    

    if sptype.startswith("O3"):
        temperature = 5
    elif sptype.startswith("O4"):
        temperature = 7
    elif sptype.startswith("O5"):
        temperature = 10
    elif sptype.startswith("O6"):
        temperature = 11
    elif sptype.startswith("O7"):
        temperature = 13
    elif sptype.startswith("O8"):
        temperature = 15
    elif sptype.startswith("O9"):
        temperature = 18
    elif sptype.startswith("B0"):
        temperature = 20
    elif sptype.startswith("B1"):
        temperature = 25
    elif sptype.startswith("B2"):
        temperature = 28
    elif sptype.startswith("B3"):
        temperature = 31
    elif sptype.startswith("B4"):
        temperature = 33
    elif sptype.startswith("B5"):
        temperature = 35
    elif sptype.startswith("B6"):
        temperature = 36
    elif sptype.startswith("B7"):
        temperature = 37
    elif sptype.startswith("B8"):
        temperature = 41
    elif sptype.startswith("B9"):
        temperature = 49
    elif sptype.startswith("A0"):
        temperature = 51
    elif sptype.startswith("A1"):
        temperature = 52
    elif sptype.startswith("A2"):
        temperature = 53
    elif sptype.startswith("A3"):
        temperature = 56
    elif sptype.startswith("A4"):
        temperature = 56
    elif sptype.startswith("A5"):
        temperature = 56
    elif sptype.startswith("A6"):
        temperature = 57
    elif sptype.startswith("A7"):
        temperature = 58
    elif sptype.startswith("A8"):
        temperature = 59
    elif sptype.startswith("A9"):
        temperature = 60
    elif sptype.startswith("F0"):
        temperature = 60
    elif sptype.startswith("F1"):
        temperature = 63
    elif sptype.startswith("F2"):
        temperature = 61
    elif sptype.startswith("F3"):
        temperature = 62
    elif sptype.startswith("F4"):
        temperature = 62
    elif sptype.startswith("F5"):
        temperature = 63
    elif sptype.startswith("F6"):
        temperature = 63
    elif sptype.startswith("F7"):
        temperature = 63
    elif sptype.startswith("F8"):
        temperature = 64
    elif sptype.startswith("F9"):
        temperature = 65
    elif sptype.startswith("G0"):
        temperature = 65
    elif sptype.startswith("G1"):
        temperature = 66
    elif sptype.startswith("G2"):
        temperature = 66
    elif sptype.startswith("G3"):
        temperature = 66
    elif sptype.startswith("G4"):
        temperature = 66
    elif sptype.startswith("G5"):
        temperature = 66
    elif sptype.startswith("G6"):
        temperature = 66
    elif sptype.startswith("G7"):
        temperature = 67
    elif sptype.startswith("G8"):
        temperature = 67
    elif sptype.startswith("G9"):
        temperature = 67
    elif sptype.startswith("K0"):
        temperature = 68
    elif sptype.startswith("K1"):
        temperature = 69
    elif sptype.startswith("K2"):
        temperature = 70
    elif sptype.startswith("K3"):
        temperature = 70
    elif sptype.startswith("K4"):
        temperature = 71
    elif sptype.startswith("K5"):
        temperature = 72
    elif sptype.startswith("K6"):
        temperature = 72
    elif sptype.startswith("K7"):
        temperature = 73
    elif sptype.startswith("K8"):
        temperature = 73
    elif sptype.startswith("K9"):
        temperature = 73
    elif sptype.startswith("M0"):
        temperature = 74
    elif sptype.startswith("M1"):
        temperature = 74
    elif sptype.startswith("M2"):
        temperature = 75
    elif sptype.startswith("M3"):
        temperature = 75
    elif sptype.startswith("M4"):
        temperature = 75
    elif sptype.startswith("M5"):
        temperature = 75
    elif sptype.startswith("M6"):
        temperature = 75
    elif sptype.startswith("M7"):
        temperature = 75
    elif sptype.startswith("M8"):
        temperature = 75
    elif sptype.startswith("M9"):
        temperature = 75
    elif sptype.startswith("O"):
        temperature = 13
    elif sptype.startswith("B"):
        temperature = 35
    elif sptype.startswith("M"):
        temperature = 75
    elif sptype.startswith("C"):
        temperature = 75
    elif sptype.startswith("A"):
        temperature = 56
    elif sptype.startswith("R"):
        temperature = 75
    elif sptype.startswith("G"):
        temperature = 66
    elif sptype.startswith("W"):
        temperature = 35
    elif sptype.startswith("K"):
        temperature = 72
    elif sptype.startswith("N"):
        temperature = 72
    elif sptype.startswith("S"):
        temperature = 72
    elif sptype.startswith("F"):
        temperature = 63
    elif sptype.startswith("DA"):
        temperature = 35
    else:
        # print(" Star_sptype_error", sptype,'----HIP No: ', hip_no)
        # Handle other spectral types or unknown cases
        temperature = 66  # Default temperature

    return temperature

def GET_DISTANCE(parallax):
    if parallax > 0:
        distance = 1000/parallax  #distance in parsec
    else:
        distance = 1e6  
    return distance

def READ_CASTELLI_SPECTRA(spec_dir = Castelli_data):
    stellar_spectra = [{"temp": None, "spectrum": None, "wavelength": None} for i in range(N_CASTELLI_MODELS)]

    temper = [50000 - i*1000 for i in range(38)] + [13000 - (i - 37)*250 for i in range(38, 76)]
    gindex = [12]*5 + [11]*6 + [10]*53 + [11]*12  # g_effective = (g_index[]-2) *5)))
    # print(list(zip(temper, gindex)))

    for i in range(len(temper)):
        filename = f"{spec_dir}/ckp00_{temper[i]}.fits"
        stellar_spectra[i]["temp"] = temper[i]
        with fits.open(filename) as hdul:
            data = hdul[1].data
            stellar_spectra[i]["spectrum"] = data.field(gindex[i] - 1)
            stellar_spectra[i]["wavelength"] = data.field(0)
    
    return stellar_spectra 

def GET_SPECTRA(spec_dir, data):
    all_spectra = READ_CASTELLI_SPECTRA(spec_dir)
    spectral_FOV = Spectral_FOV()

    wave_min , wave_max = read_parameter_file()
    print(wave_min, wave_max)
    low_lim = index_greater_than(all_spectra[0]['wavelength'], wave_min)
    high_lim = index_greater_than(all_spectra[0]['wavelength'], wave_max)

    for i in range(len(data)): # repeats over all frames i
        spectral_FOV.frame.append(str(i+1))
        spectral_flux = []
        spectral_wavelength = []
        scale_per_star = []
        photons_per_star = []

        _, d, frame_size = zip(data[i])
        if d[0]:
            c = list(zip(*d[0]))
            ra, dec = c[0], c[1]
            spectral_FOV.ra.append(ra)
            spectral_FOV.dec.append(dec)
            spectral_FOV.frame_size = frame_size
            spectral_type = c[7]

            for j in range(len(spectral_type)):     # repeats over all stars in the FOV in frame i
                t_index = GET_STAR_TEMP(spectral_type[j])
                data1 = all_spectra[t_index]
                wavelengths = data1['wavelength']
                flux = data1['spectrum']
                scale, photons = GET_SCALE_FACTOR(c[4][j], c[5][j], c[6][j], wavelengths, flux)

                spectral_wavelength = wavelengths[low_lim:high_lim]
                spectral_flux.append(flux[low_lim:high_lim])
                scale_per_star.append(scale)
                photons_per_star.append(photons[low_lim:high_lim])
            spectral_FOV.spectra_per_star.append(spectral_flux)
            spectral_FOV.wavelength.append(spectral_wavelength)
            spectral_FOV.scale.append(scale_per_star)
            spectral_FOV.photons.append(photons_per_star)


        else:
            # print('Frame',i +1,'is EMPTY', end="\n \n")
            spectral_FOV.wavelength.append([0])
            spectral_FOV.spectra_per_star.append([[0]])
            spectral_FOV.ra.append([0])
            spectral_FOV.dec.append([0])
            spectral_FOV.scale.append([0])
            spectral_FOV.photons.append([[0]])
    
    return spectral_FOV

def Trim_Spectral_data(data, photons):
    wave_min , wave_max = read_parameter_file()
    # print(wave_min, wave_max)
    low_lim = index_greater_than(data['wavelength'], wave_min)
    high_lim = index_greater_than(data['wavelength'], wave_max)
    return data['wavelength'][low_lim: high_lim], data['spectrum'][low_lim: high_lim], photons[low_lim: high_lim]

def GET_SCALE_FACTOR( V_mag, parallax, B_V, waveL_range, stellar_spectra):

    # V_mag, parallax, B_V= c[4][j], c[5][j], c[6][j]
    scale = 0
    tot_photons = []

    if V_mag == 0:
        scale = 0
        tot_photons.append(0)
    else:
        vindex = index_greater_than(waveL_range, 5450)
        bindex = index_greater_than(waveL_range, 4360)
        # windex = next(i for i, w in enumerate(stellar_spectra[sindex]['wavelength']) if w >= inp_par['wave'])
        # print(vindex, bindex, '\n' ,stellar_spectra)
        # print(stellar_spectra[vindex])
        bflux = stellar_spectra[bindex]
        vflux = stellar_spectra[vindex]
        # print (f'bflux:{bflux}, vflux:{vflux}')

        b_mag = -2.5 * math.log10(bflux / 6.61)
        v_mag = -2.5 * math.log10(vflux / 3.64)
        # B_V = hipstar['B_mag'] - hipstar['V_mag']

        ebv = B_V - (b_mag - v_mag)
        ebv = max(ebv, 0)
        distance = GET_DISTANCE(parallax)

        scale = 3.64e-9 * pow(10, -0.4 * (V_mag - 3.1 * ebv)) / vflux
        scale = scale * distance**2

        for w in  range (len(waveL_range)):
            photon_number = stellar_spectra[w] * scale * 4 * math.pi * ERG_TO_PHOT * waveL_range[w]
            if (photon_number <= 1e-10):
                photon_number = 0

            tot_photons.append(photon_number)

        # if hipstar['HD_NO'] in [158926, 160578]:
        #     print(f"{hipstar['HIP_NO']} {hipstar['sp_type']} {stellar_spectra[sindex]['filename']} {sindex} {hipstar['V_mag']} {bflux} {vflux} {hipstar['scale']} {stellar_spectra[sindex]['spectrum'][windex]}")
    # print(scale, tot_photons)
    return scale, tot_photons

def GET_All_STAR_TEMP(sptype_list, hip_no_list):
    t_index = []
    for i ,sptype in enumerate(sptype_list):
        t_index.append(GET_STAR_TEMP(sptype, hip_no_list[i]))
    return t_index

def GET_ALL_STAR_DISTANCE(parallax_list):
    distance = []
    for parallax in parallax_list:
        distance.append(GET_DISTANCE(parallax))
    return distance

def GET_SCALE_FACTOR_ALL_STARS(all_spectra):
    return 0



# # Example usage
# hipline = {"sp_type": "sdG2V", "temperature": 0}
# temperature = GET_STAR_TEMP(hipline['sp_type'])
# print(f"Temperature: {temperature}")

# all_spectra = READ_CASTELLI_SPECTRA(Castelli_data)
# stellar_spectra = StellarSpectrum()
# data1 = all_spectra[temperature]
# # print(data1)
# stellar_spectra.temperature = float(data1['temp'])
# stellar_spectra.wavelength = np.array(data1['wavelength'])
# stellar_spectra.spectrum = np.array(data1['spectrum'])
# print(stellar_spectra.temperature , list(zip(stellar_spectra.wavelength, stellar_spectra.spectrum)))
# # for i in range(len(stellar_spectra.wavelength)):
# #     stellar_spectra.spectrum[i] = 3.336e-19 * (4*np.pi)**(-1) * (stellar_spectra.wavelength[i])**2
# # print(list(zip(stellar_spectra.wavelength, stellar_spectra.spectrum)))
# # print (stellar_spectra.temperature,'\n',
# #         stellar_spectra.wavelength,'\n',
# #         stellar_spectra.spectrum)

# low_UV = index_greater_than(stellar_spectra.wavelength, 100)
# high_UV = index_greater_than(stellar_spectra.wavelength, 3800)
# low_vis = index_greater_than(stellar_spectra.wavelength, 3800)
# high_vis = index_greater_than(stellar_spectra.wavelength, 7500)

# fig, ax = plt.subplots()
# # ax.plot(wavelength[low_UV:high_UV], Surface_Flux[low_UV:high_UV], color='blue', label = r'ckp00_3500')
# ax.plot(stellar_spectra.wavelength[low_vis:high_vis], stellar_spectra.spectrum[low_vis:high_vis], color='blue', label = stellar_spectra.temperature)
# # ax.plot(stellar_spectra.wavelength, stellar_spectra.spectrum, color='blue', label = stellar_spectra.temperature)

# ax.set_xlabel(r'Wavelength- A')
# ax.set_ylabel(r'Surface Flux')
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# plt.savefig( 'C:\\Users\\Akshank Tyagi\\Desktop\\without_conversion.png' )
# plt.show()

# contents = []
# import pickle

# # Load data from a pickle file
# with open('star_data.pkl', 'rb') as file:
#     contents = pickle.load(file)
# print(len(contents))
# spectral_fov = GET_SPECTRA(Castelli_data, contents)
# # print(spectral_fov.photons)

# for i in range(len(spectral_fov.frame)):
#     fig, ax = plt.subplots()
#     for j in range(len(spectral_fov.photons[i])):
#         ax.plot(spectral_fov.wavelength[i], spectral_fov.photons[i][j], label = f'frame:{spectral_fov.frame[i]} -- star{j+1}')
#     ax.set_ylim(1e-3, 1e10)
#     ax.set_xlabel(r'Wavelength- A')
#     ax.set_ylabel(r'# of Photons')
#     # ax.set_xscale('log')
#     ax.set_yscale('log')
#     ax.legend()
#     plt.show()

