from configparser import ConfigParser
import os


def get_folder_loc():
    # print(os.getcwd())
    # folder_loc = fr'/home/akshanktyagi/Documents/GitHUB/UV-Sky-Simulations-Backup'  
    folder_loc = fr'/home/akshanktyagi/Documents/GitHUB/UV-Sky-Simulations-Backup' 
    folder_loc = os.getcwd()
    if not folder_loc.endswith(os.sep):
        folder_loc += os.sep
    params_file = fr'{folder_loc}init_parameter.txt'

    return folder_loc, params_file 

folder_loc,  params_file = get_folder_loc()


config = ConfigParser()

config['Params_1'] = {
    'hipparcos_catalogue' : f'{folder_loc}hip_main.dat', #path to the Hipparcos file
    'Castelli_data' : f'{folder_loc}Castelli{os.sep}ckp00', #path to the ckp00 file of the Castelli Kurucz Atlas
    'dust_file' : f'{folder_loc}green_dust_60{os.sep}green_filled_3d.fits', 
    'dust_col_file' : f'{folder_loc}green_dust_60_col{os.sep}green_filled_cum_3d.fits',
    'sigma_file' : f'{folder_loc}kext_albedo_WD_MW_3.1_60_D03.all.txt',    
    'exclude_stars' : f'{folder_loc}exclude_stars_final.txt',

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
    'height': 2, #Dec height
    'star_mag_threshold' : 20, #threshold for what apaarent magnitude stars we want to look at
    
    # Spectrum Parameters (UV Band Wavelengths in Angstroms)
    'limit_min': 100,
    'limit_max': 4000,

    #Animation parameters
    # set view
    'azm': 40,
    'ele': 25,
    'longitudinal_spectral_width' : 0.8, #Declination width of spectral spread to fall on detector in degrees
    'interval_bw_Frames' : 2000 # milliSec
}
config['Scatter_params'] = {
    # Dust Scatter Parameters
    "No_Photons":  10000000, 
    'wavelength': [1105],
    'No_scatter': 5,
    'Albedo': 0.36, #0.36, 0.40
    'Phase_func': 0.5, #0.5, 0.6
    'print_debug': "no",
    'min_gl_debug':0,
    'max_gl_debug': 360,
    'min_gb_debug': -90,
    'max_gb_debug': 90,
}

config['WCS'] = {
    'CRVAL1':0,
    'CRVAL2':0,
    'CRPIX1':1800,
    'CRPIX2':900,
    'CDELT1':-.1,
    'CDELT2':-.1,
    'CROTA':0,
    'CTYPE': '-AIT',
    'NAXIS1':3600,
    'NAXIS2':1800
}

with open(params_file,"w") as f:
    config.write(f)
    print(f'{folder_loc}init_parameter.txt Saved')