# Functions to calculate state vectors, ra, dec and plot
# Author: Ravi Ram
# Jayant Murthy Sirs code ascl:1512.012

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime, timedelta
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72

from configparser import ConfigParser
config = ConfigParser()
import pickle

def get_folder_loc():
    folder_loc = r'C:\Users\Akshank Tyagi\Documents\GitHub\UV-Sky-Simulations\\'
    return folder_loc

from Params_configparser import *
from Satellite_configparser import *
from star_data import *#, filter_by_fov, read_hipparcos_data
from plot import *#,animate
from star_spectrum import *#,GET_SPECTRA

# include the parameter file and sattelite TLE file
folder_loc = get_folder_loc()
params_file = f'{folder_loc}init_parameter.txt'
sat_file = f'{folder_loc}Satellite_TLE.txt'
print(params_file)

def read_parameter_file(filename=params_file, param_set = 'Params_1'):
    config.read(filename)

    global roll_rate_hrs

    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    sat_name = config.get(param_set, 'sat_name')
    roll = config.get(param_set, 'roll')
    if (roll == True):
        roll_rate_hrs = float(config.get(param_set, 'roll_rate_hrs'))
    else:
        roll_rate_hrs = False
    N_revolutions = config.get(param_set, 'number of Revolutions')
    N_frames = config.get(param_set, 'N_frames')
    T_slice = config.get(param_set, 't_slice')
    # width = float(config.get(param_set, 'width'))
    # height = float(config.get(param_set, 'height'))
    # Threshold = float(config.get(param_set,'star_mag_threshold'))
    
    print('sat_name:', sat_name, ', roll:',roll,',  roll_rate_hrs:',roll_rate_hrs, ',  N_revolutions:',N_revolutions, ',  N_frames:', N_frames, ',  T_slice:', T_slice)
    return hipp_file, castelli_file, sat_name, float(T_slice), N_frames, float(N_revolutions), roll #, Threshold

def read_satellite_TLE(filename= sat_file, sat_name = 'ISS'):
    config.read(filename)
    line1 = config.get(sat_name, 'line1')
    line2 = config.get(sat_name, 'line2')

    return line1, line2

# get satellite object from TLE (2 lines data)
def get_satellite(line1, line2):
    global mu, r, a
    
    # create satellite object from TLE
    satellite = twoline2rv(line1, line2, wgs72)
    # constants
    mu = satellite.mu           # Earth’s gravitational parameter (km³/s²)
    r = satellite.radiusearthkm # Radius of the earth (km).
    # orbital parameters
    a = satellite.a * r
    #e = satellite.ecco
    apo, peri = satellite.alta * r, satellite.altp * r
    print('Perigee : %5.2f km, Apogee : %5.2f km' % (peri, apo))
    print('mu =',mu, 'km^3/s^2, Earth Radius =',r, 'km \nDistances from Center of Earth: Perigee =',r+peri,'km, Semi major =', a, 'km, Apogee =',r + apo,'km' )
    # perigee and apogee
    # return
    return satellite

# get ra and dec from state vectors
def get_ra_dec_from_sv(r, v):
    # norm
    v_n = np.linalg.norm(v)
    # print(v_n)
    # direction cosines
    l = v[0]/v_n; m = v[1]/v_n; n = v[2]/v_n; 
    # print(l,m,n)
    # declination
    delta = np.arcsin(n)*180/np.pi                    
    # right ascension
    np.cosd = lambda x : np.cos( np.deg2rad(x) )
    if m >0:  alfa = np.arccos(l/np.cosd(delta))*180/np.pi
    else: alfa = 360 - np.arccos(l/np.cosd(delta))*180/np.pi
    # return
    return alfa, delta,

def get_cords_from_ra_dec (ra, dec, distance):
    x = np.cos(np.deg2rad(ra))*np.cos(np.deg2rad(dec))*distance
    y = np.sin(np.deg2rad(ra))*np.cos(np.deg2rad(dec))*distance
    z = np.sin(np.deg2rad(dec))*distance
    return x, y, z

# returns a list of state vectors, ra, dec for a
# given sgp4 satellite object
def propagate(sat, time_start, time_end, dt):
    # time
    end = np.arange(0.0, time_end, dt)
    # print(end)
    # list of datetime
    time_arr = time_start + end.astype('timedelta64[s]')    
    # state vectors, ra, dec for each time step
    position = []; velocity = []
    right_ascension = []; declination = []
    for j in time_arr.tolist():
        p, v = sat.propagate(j.year, j.month, j.day, j.hour, j.minute, j.second)
        # print(p,v)
        ra, dec = get_ra_dec_from_sv(p, v)
        # list packing
        position.append(p); velocity.append(v)
        right_ascension.append(ra); declination.append(dec)
        

    # slice into columns
    pos, vel   = list(zip(*position)), list(zip(*velocity))
    # print(pos, vel)
    # print (position)
    X, Y, Z    = np.array(pos[0]), np.array(pos[1]), np.array(pos[2])
    VX, VY, VZ = np.array(vel[0]), np.array(vel[1]), np.array(vel[2])
    state_vectors = [X, Y, Z, VX, VY, VZ]
    celestial_coordinates = [np.array(right_ascension), np.array(declination)]
    # print(celestial_coordinates)
    # return
    return time_arr, state_vectors, celestial_coordinates

# get list of star data in view along with satellite state_vectors
def get_simulation_data(sat, df, start_time, sim_secs, time_step, roll=False):
    # state_vectors, celestial_coordinates
    tr, sc, cc = propagate(sat, start_time, sim_secs, time_step)
    # parse celestial_coordinates
    ra, dec = cc

    # [TESTING] Roll about velocity direction 
    if roll:
        if roll_rate_hrs:         # deg per hr
            roll_rate_sec = roll_rate_hrs/3600.0 # deg per sec
        # modify ra with roll rate
            ra = [ri + roll_rate_sec * i for i, ri in enumerate (ra)]

    # Create an empty all frames data
    frame_row_list = []

    for frame, (r, d) in enumerate(zip(ra, dec)):
        # print(frame, (r, d))
        tdf_values, frame_boundary = filter_by_fov(df, r, d)
        tdf_values = tdf_values.values.tolist()
        # print (tdf_values)
        frame_row_list.append([frame, tdf_values, frame_boundary ])
    return tr, sc, frame_row_list,


def main():
    global data
    hipp_file, castelli_dir, sat_name, t_slice, n_frames, N_revolutions, roll  = read_parameter_file(params_file,'Params_1')
    line1, line2 = read_satellite_TLE(sat_file, sat_name)
    
    # create satellite object
    satellite = get_satellite(line1, line2)
    # read star data
    df = read_hipparcos_data(hipp_file)
    print (df[:20])
    # time period for one revolution
    t_period = N_revolutions* 2 * np.pi * (a**3/mu)**0.5
    print("Time period of Satellite =",t_period,'sec')

    # each time slice
    if t_slice:
        t_step = int(t_period / t_slice) + 1
        print("N_Frames:", t_step,"t_slice =",t_slice)
    else:
        # set Number of frames
        if n_frames:
            t_slice = t_period/int(n_frames)
            print("N_Frames:", n_frames ,"t_slice =",t_slice)
        else:
            print('T_slice not found')

    # simulation starts from current time to one full orbit
    start = np.datetime64(datetime.datetime.now())
    # times, state_vectors, celestial_coordinates  
    time_arr, state_vectors, celestial_coordinates = get_simulation_data(satellite, df, start, t_period, t_slice, roll)
    Spectra = GET_SPECTRA(castelli_dir, celestial_coordinates)

    # animate
    animate(time_arr, state_vectors, celestial_coordinates, Spectra, r)
    return

# main
if __name__ == '__main__':
    main()
    
    