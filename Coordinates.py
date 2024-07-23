import math
import numpy as np
from configparser import ConfigParser
from astropy.io import fits
from astropy.wcs import WCS
import csv
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mc

from Params_configparser import get_folder_loc

folder_loc, _ = get_folder_loc()


def conv_cart_to_eq(xp, yp, zp, dust_par):
    x = (xp - dust_par.sun_x) * dust_par.dust_binsize
    y = (yp - dust_par.sun_y) * dust_par.dust_binsize
    z = (zp - dust_par.sun_z) * dust_par.dust_binsize
    
    dist = np.sqrt(x**2 + y**2 + z**2)
    dxy = np.sqrt(x**2 + y**2)
    sdec = z / dist
    cra = x / dxy
    sra = y / dxy
    
    if sra >= 0:
        ra = np.degrees(np.arccos(cra))
    else:
        ra = 360 - np.degrees(np.arccos(cra))
    dec = np.degrees(np.arcsin(sdec))
    
    return ra, dec

def conv_eq_to_cart(ra, dec, distance):
    x = np.cos(np.deg2rad(ra))*np.cos(np.deg2rad(dec))*distance
    y = np.sin(np.deg2rad(ra))*np.cos(np.deg2rad(dec))*distance
    z = np.sin(np.deg2rad(dec))*distance
    return x, y, z

def conv_eq_to_gal(ra, dec):
    # The transformation matrix
    A = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762]
    ])

    # Convert input angles from degrees to radians
    r = np.radians(ra)
    d = np.radians(dec)

    # Calculate x vector
    cos_d = np.cos(d)
    x = np.array([np.cos(r) * cos_d, np.sin(r) * cos_d, np.sin(d)])

    # Matrix multiplication to get xp vector
    xp = A @ x

    # Calculate latitude and longitude in radians, ensuring xp[2] is within [-1, 1]
    lat = np.arcsin(np.clip(xp[2], -1, 1))
    lon = np.arctan2(xp[1], xp[0])
    
    # Adjust longitude to be in the range [0, 2*pi]
    lon = np.mod(lon, 2.0 * np.pi)

    # Convert latitude and longitude to degrees
    gl = np.degrees(lon)
    gb = np.degrees(lat)

    return gl, gb


def conv_gal_to_eq(gl, gb):
    # Convert Galactic coordinates from degrees to radians
    l = np.radians(gl)
    b = np.radians(gb)

    # Calculate Cartesian coordinates in Galactic system
    xyz_galactic = np.array([np.cos(l) * np.cos(b), np.sin(l) * np.cos(b), np.sin(b)])

    # Apply the inverse of the transformation matrix to get Equatorial coordinates
    A_inv = np.linalg.inv(A)
    xyz_equatorial = np.dot(A_inv, xyz_galactic)

    # Convert Equatorial coordinates from Cartesian to spherical (RA, Dec)
    ra_rad = np.arctan2(xyz_equatorial[1], xyz_equatorial[0])
    dec_rad = np.arcsin(xyz_equatorial[2])

    # Convert Equatorial coordinates from radians to degrees
    ra = np.degrees(ra_rad)
    dec = np.degrees(dec_rad)
    for i, ra_i in enumerate(ra):
        if ra_i<0:
            ra[i] =360 + ra_i
    
    return ra, dec



def get_world_coordinates(x, y, fits_file):
    with fits.open(fits_file) as hdul:
        wcs = WCS(hdul[0].header)  # Create WCS object from FITS header
        print(wcs)
        # Convert pixel coordinates to world coordinates (RA, Dec)
        ra, dec = wcs.all_pix2world(x, y, 1)  # Assumes pixel indices start from 1
        # print(ra,dec)
        return ra, dec

def plot_diffused_bg(data, wavelength, nphoton, a, g):

    data= data
    plt.figure(figsize=(5, 3))
    n = np.random.rand()
    colors = [(0, 0, 0), (0, 0, 1)]  # Black to blue
    cmap_name = 'black_to_blue'
    BtoB_cmap = mc.LinearSegmentedColormap.from_list(cmap_name, colors)
    # print(wavelength)
    plt.imshow(data, cmap= BtoB_cmap, vmin=0, vmax= 1000)
    plt.colorbar()
    plt.title(f'diffused_UV_background@{wavelength} for {nphoton}')
    # plt.savefig(fr'{folder_loc}diffused_data{os.sep}diffused_UV_BG_{wavelength}.jpg', dpi=300)
    plt.savefig(fr'{folder_loc}diffused_output{os.sep}trialN{int(nphoton)}_{int(wavelength)}.jpg', dpi=1000)

    plt.show(block=False)  # Show the plot non-blocking
    plt.pause(2)           # Pause for 2 seconds
    plt.close()            
