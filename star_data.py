# Functions to read and filter sky catalogue data
# Author: Ravi Ram

import datetime
import numpy as np
import pandas as pd
from configparser import ConfigParser
from star_spectrum import *#,GET_SPECTRA
from Params_configparser import get_folder_loc
folder_loc = get_folder_loc()

# Hipparcos Catalogue [hip_main.dat]
# http://cdsarc.u-strasbg.fr/ftp/cats/I/239/ 
# FILENAME = r'C:\Users\Akshank Tyagi\Documents\GitHub\spg-iiap-UV-Sky-Simulation\hip_main.dat'
params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global hipp_file , star_mag_threshold
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    # Field of View size:
    width = float(config.get(param_set, 'width'))
    height = float(config.get(param_set, 'height'))
    star_mag_threshold = float(config.get(param_set, 'star_mag_threshold'))
    return  width, height 

_, _ = read_parameter_file()

def get_abc(hip):
    return 3*hip

# read hipparcos catalogue 'hip_main.dat'
def read_hipparcos_data(FILENAME = hipp_file, mag_threshold = True):
    # Field H1: Hipparcos Catalogue (HIP) identifier
    # Field H5: V magnitude
    # Fields H8–9:  The right ascension, α , and declination, δ (in degrees)    
    threshold = star_mag_threshold
    try:
        df = pd.read_csv(FILENAME, header=None,
                         sep = '|', skipinitialspace=True).iloc[:, [1, 5, 8, 9, 11, 37, 76, 32, 34]]
        df.columns = ['hip', 'mag', 'ra_deg', 'de_deg', 'trig_parallax', 'B-V', 'Spectral_type', 'B_mag', 'V_mag']

        df['mar_size'] = 2*(threshold - df['mag'])
        df['distance'] = GET_ALL_STAR_DISTANCE(df['trig_parallax'])
        df['t_index'] = GET_All_STAR_TEMP(df['Spectral_type'], df['hip'])
        
        # filter data above
        if mag_threshold == True:
            print (f'Stars apparent magnitude Threshold= {threshold}')
            q = 'mag <= @threshold'
            df = df.query(q) 
        return df  
    
    except FileNotFoundError:
        print("df is empty. File not found.")

# camera fov : 9.31◦ × 7◦
def filter_by_fov(mdf, ra, de): 
    # frame field of view
    # get valid boundaries  
    width, height = read_parameter_file()
    xmin, ymin, xmax, ymax = get_frame_boundaries( width, height, ra, de)
    frame_boundaries = [xmin, ymin, xmax, ymax]
    # print(frame_boundaries)
    # if mdf[0]:
    # extract useful columns
    mdf = mdf[['ra_deg', 'de_deg', 'mar_size','hip','mag', 'trig_parallax', 'B-V', 'Spectral_type']]
    # filter data within the boundaries    
    q = 'ra_deg >= @xmin & ra_deg <= @xmax & de_deg >= @ymin & de_deg <= @ymax' 
    mdf = mdf.query(q)
    # print(mdf)
    # return filtered data
    return mdf, frame_boundaries


# get valid frame boundaries
def get_frame_boundaries(w, h, x, y):
    # set x boundaries
    xmin = float(x) - float(w) / 2.0
    xmin = 0 if xmin<0 else xmin   
    if xmin==0: xmax = w
    else: xmax = float(x) + float(w) / 2.0
    xmax = 360 if xmax>360 else xmax
    if xmax==360: xmin = 360-w
    
    # set y boundaries
    ymin = float(y) - float(h) / 2.0
    ymin = -90 if ymin<-90 else ymin
    if ymin==-90: ymax = (-90 + float(h))
    else: ymax = float(y) + float(h) / 2.0
    if ymax == 90: ymin = (90 - float(h))
    
    # return limits 
    return  xmin, ymin, xmax, ymax
