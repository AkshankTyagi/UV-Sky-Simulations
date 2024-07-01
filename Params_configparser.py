from configparser import ConfigParser

# from view_orbit import get_folder_loc

def get_folder_loc():
    folder_loc = r'C:\Users\Akshank Tyagi\Documents\GitHub\UV-Sky-Simulations\\'
    print('1')
    return folder_loc

folder_loc = get_folder_loc()
config = ConfigParser()

config['Params_1'] = {
    'hipparcos_catalogue' : f'{folder_loc}hip_main.dat', #path to the Hipparcos file
    'Castelli_data' : f'{folder_loc}Castelli\ckp00', #path to the ckp00 file of the Castelli Kurucz Atlas
    'dust_file' : f'{folder_loc}green_dust_60\green_dust_3d.fits', 
    'dust_col_file' : f'{folder_loc}green_dust_60_col\green_cum_tau_3d.fits',
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
    'star_mag_threshold' : 7, #threshold for what apaarent magnitude stars we want to look at
    
    # Spectrum Parameters (UV Band Wavelengths in Angstroms)
    'limit_min': 100,
    'limit_max': 1110,

    #Animation parameters
    # set view
    'azm': 40,
    'ele': 25,
    'longitudinal_spectral_width' : 0.8, #Declination width of spectral spread to fall on detector in degrees
    'interval_bw_Frames' : 2000 # milliSec
}
config['Scatter_params'] = {
    # Dust Scatter Parameters
    "No_Photons":  1000, 
    'wavelength': [1500],
    'No_scatter': 5,
    'Albedo': 0.36,
    'Phase_func': 0.5,
    'print_debug': "no",
    'min_gl_debug':0,
    'max_gl_debug': 360,
    'min_gb_debug': -90,
    'max_gb_debug': 90,
}

config['WCS'] = {

}

with open(f'{folder_loc}init_parameter.txt',"w") as f:
    config.write(f)
    print(f'{folder_loc}init_parameter.txt Saved')