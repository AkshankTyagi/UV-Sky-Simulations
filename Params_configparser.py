from configparser import ConfigParser
from view_orbit import get_folder_loc
folder_loc = get_folder_loc()

config = ConfigParser()

config['Params_1'] = {
    'hipparcos_catalogue' : f'{folder_loc}hip_main.dat', #path to the Hipparcos file
    'Castelli_data' : f'{folder_loc}Castelli\ckp00', #path to the ckp00 file of the Castelli Kurucz Atlas
    'dust_file' : f'{folder_loc}green_dust_60\green_dust_60.fits', 
    'dust_col_file' : f'{folder_loc}green_dust_60_col\green_cum_tau_60.fits',
    'sigma_file' : f'{folder_loc}kext_albedo_WD_MW_3.1_60_D03.all.txt',    

    'sat_name' : 'Astrosat',
    'roll' : False,
    'roll_rate_hrs' : False,
    # TBA Directional Cosines of Detector from the velocity of Satellite 
    'number of Revolutions' : 1,
    # Specify either Number of frames or period in sec after which the next Frame is given
    'N_frames' : False,
    't_slice' : 400, # Seconds,

    # Camera Field of View in Deg default 9.3 X 7
    'width': 0.5, #RA width
    'height': 7, #Dec height
    'star_mag_threshold' : 15, #threshold for what apaarent magnitude stars we want to look at
    
    # Spectrum Parameters (UV Band Wavelengths in Angstroms)
    'limit_min': 100,
    'limit_max': 3800,

    # Dust Scatter Parameters
    "No_Photons":  10000000, 
    'wavelength': [1500,2000,3000] ,
    'No_scatter': 5,
    'Albedo': 0.36,
    'Phase_func': 0.5,

    #Animation parameters
    # set view
    'azm': 40,
    'ele': 25,
    'longitudinal_spectral_width' : 0.8, #Declination width of spectral spread to fall on detector in degrees
    'interval_bw_Frames' : 2000 # milliSec
}


with open(f'{folder_loc}init_parameter.txt',"w") as f:
    config.write(f)