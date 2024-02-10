from configparser import ConfigParser

config = ConfigParser()

config['Params_1'] = {
    'Hipparcos_catalogue' : 'hip_main.dat',
    'sat_name' : 'Astrosat',
    'roll' : False,
    'roll_rate_hrs' : False,
    # Directional Cosines of Detector from the velocity of Satellite 
    'number of Revolutions' : 1,
    # Specify either Number of frames or period in sec after which the next Frame is given
    'N_frames' : False,
    't_slice' : 300, # Seconds,

    # Camera Field of View in Deg default 9.3 X 7
    'width': 0.5, #RA width
    'height': 7, #Dec height
    'Star_mag_threshold' : 8.5, #threshold for what apaarent magnitude stars we want to look at

    #Animation parameters
    # set view
    'azm': 40,
    'ele': 25,
    'interval_bw_Frames' : 500 # milliSec
}

config['Params_2'] ={
    'Hipparcos_catalogue' : 'hip_main.dat',
    'sat_name' : 'RISAT-2B',
    'roll' : False,
    'roll_rate_hrs' : False,
    # Directional Cosines of Detector from the velocity of Satellite 
    'number of Revolutions' : 1,
    # Specify either Number of frames or period in sec after which the next Frame is given
    'N_frames' : False,
    't_slice' : 100, # Seconds

    # Camera Field of View in Deg
    'width': 9.31, #RA width
    'height': 7, #Dec height
    'Star_mag_threshold' : 8.5,

    #Animation parameters
    # set view
    'azm': 60,
    'ele': 55,
    'interval_bw_Frames' : 50 # milliSec
}

with open('init_parameter.txt',"w") as f:
    config.write(f)