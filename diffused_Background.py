import math
import numpy as np
from configparser import ConfigParser
from astropy.io import fits
from astropy.wcs import WCS
import csv
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mc

from Params_configparser import get_folder_loc
# from Coordinates import *

# def DEGRAD(deg):
#     return math.radians(deg)

# def RADDEG(rad):
#     return math.degrees(rad)

folder_loc , params_file = get_folder_loc()
# params_file = f'{folder_loc}init_parameter.txt'


def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global sat_name, Interval, spectra_width
    sat_name = config.get(param_set, 'sat_name')
    diffused_wavelength = config.get(param_set, 'BG_wavelength')
    return diffused_wavelength

def invert_gal_to_eq(gl, gb):

    # Convert Galactic coordinates from degrees to radians
    l = np.radians(gl)
    b = np.radians(gb)

    # Transformation matrix from Galactic to Equatorial coordinates, // J.C.Liu et al(A&A 536, A102 (2011))//wrong
    A = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762]
    ])

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
        # print(hdul[0].header)
        # print(hdul[1].header)
        # print(hdul[2].header)
        # print(ra,dec)
        return ra, dec

def plot_diffused_bg(data, wavelength, nphoton):

    data= data/3000

    n = np.random.rand()
    colors = [(0, 0, 0), (0, 0, 1)]  # Black to blue
    cmap_name = 'black_to_blue'
    BtoB_cmap = mc.LinearSegmentedColormap.from_list(cmap_name, colors)
    # print(wavelength)

    plt.figure(figsize=(10, 6))
    plt.imshow(data, cmap= BtoB_cmap, vmin=0, vmax= 1)
    plt.gca().invert_yaxis()
    # plt.colorbar()
    plt.title(f'diffused_UV_background@{wavelength} for {nphoton}')
    # plt.savefig(fr'{folder_loc}diffused_data{os.sep}diffused_UV_BG_{wavelength}.jpg', dpi=300)
    plt.savefig(fr'{folder_loc}diffused_data{os.sep}trial1000N{nphoton}_{wavelength}.jpg', dpi=300)

    plt.show(block=False)  # Show the plot non-blocking
    plt.pause(115)           # Pause for 2 seconds
    plt.close() 


# Example usage
# wavelength = 1105
# # file_path = fr"C:\Users\Akshank Tyagi\Documents\GitHub\UV-Sky-Simulations\diffused_data\RA_sorted_flux_{wavelength}.csv"

# print ('working')
for x in [1100]:
    nphoton = 000000
    # fits_filename = f"{folder_loc}diffused_data{os.sep}scattered_1e10_{x}_a40_g6{os.sep}scattered.fits"
    # fits_filename = r"C:\Users\Akshank Tyagi\Documents\GitHub\UV-Sky-Simulations\diffused_data\scattered_1e10_1100_a40_g6\scattered.fits"
    fits_filename = fr'{folder_loc}diffused_data{os.sep}scattered_{nphoton}[(1100, 1110)]_mag20.fits'
    ra, dec = get_world_coordinates( 1800, 900, fits_filename)
    print(f"ra,dec:",ra,dec)
    with fits.open(fits_filename) as hdul:
        data = hdul[0].data
        # df = pd.DataFrame(data)
        # print(df)
        # df.to_csv(f'my_scatter_output{nphoton}.csv', index=False)
        plot_diffused_bg(data, 1100, nphoton)



# if os.path.exists(file_path):
#     print (f'{file_path} file exists.')
# else:
#     print('Running diffused_Background.py to create sorted files of diffused UV BG from Jayant Murthy (2016) data')
#     # wavelength = read_parameter_file()
#     for wavelength in [1100, 1500, 2300]:
#         fits_filename = f"{folder_loc}diffused_UV_data\scattered_1e10_{wavelength}_a40_g6\scattered.fits"
#         print(fits_filename)

#         plot_diffused_bg(fits_filename, wavelength)

#         gl= [0] #longitude
#         gb= [0] #latitude
#         ra, dec = invert_gal_to_eq(gl,gb)
#         print (gl,gb,'--->',ra,dec)

#         x_pixel = [1800]
#         y_pixel = [900]

#         glon, glat = get_world_coordinates(x_pixel, y_pixel, fits_filename)
#         ra, dec = invert_gal_to_eq(glon,glat)
#         print("RA:", ra)
#         print("Dec:", dec)

#         x_range= 3600
#         x_array = np.arange(1, 3601)
#         y_range = 1800

#         # Obtain the Flux values
#         with fits.open(fits_filename) as hdul:
#             global values
#             # Print the header of the PrimaryHDU (HDU 0)
#             fits.info(fits_filename)
#             values = hdul[0].data

#         print("values_obtained")

#         # from the Pixels converting to galctic to Equatorial coordinates + flux grid
#         grid = []
#         for i in range( y_range+1):
#             if i%20== 0:
#                 print("line-",i)
#             glon, glat = get_world_coordinates(x_array, [i]*3600, fits_filename)
#             # print(glon, glat)
#             ra_line, dec_line = invert_gal_to_eq(glon,glat)
            
#             line = zip(values[i-1],ra_line, dec_line)
#             grid.append(line)

#         print("cordinate_transformation_done")

#         # write the csv file with flux and respective coordinate
#         with open(fr"{folder_loc}diffused_UV_data\flux_{wavelength}.csv", 'w', newline='') as csvfile:
#             csv_writer = csv.writer(csvfile)
#             # Write the data array to the CSV file row by row
#             for row in grid:
#                 csv_writer.writerow(row)

#         print(f"flux_{wavelength}.csv saved")
#         del values
#         del grid
#         del x_array

#         list = []
#         with open(fr"{folder_loc}diffused_UV_data\flux_{wavelength}.csv", 'r', newline='') as csvfile:
#             csv_reader = csv.reader(csvfile)
#             for r,row in enumerate(csv_reader):
#                 print(r)
#                 for flux_location in row:
#                     if flux_location == '(o.o, nan, nan)':
#                         continue

#                     flux, ra, dec = flux_location.split(', ')

#                     if ra == 'nan':
#                         continue

#                     flux = float(flux[1:])
#                     ra = float(ra)
#                     dec = float(dec[:-1])
#                     list.append([ra, dec, flux])

#         sorted_entries = sorted(list, key= lambda x: x[0])

#         with open(fr"{folder_loc}diffused_UV_data\RA_sorted_flux_{wavelength}.csv", 'w', newline='') as csvfile:
#             csv_writer = csv.writer(csvfile)
#             # Write the data array to the CSV file row by row
#             for row in sorted_entries:
#                 csv_writer.writerow(row)

#         print(f"RA_sorted_flux_{wavelength}.csv saved")


# sorted_list=[]
# with open(fr"{folder_loc}diffused_UV_data/flux_{wavelength}.csv", 'r', newline='') as csvfile:
#     csv_reader = csv.reader(csvfile)
#     already_sorted = 0
#     for r,row in enumerate(csv_reader):
#         print(r, len(sorted_list),'---',len(sorted_list)-already_sorted)
#         already_sorted = len(sorted_list)

#         for flux_location in row:

#             if flux_location == '(o.o, nan, nan)':
#                 continue

#             flux, ra, dec = flux_location.split(', ')

#             if ra == 'nan':
#                 continue

#             flux = float(flux[1:])
#             ra = float(ra)
#             dec = float(dec[:-1])
#             if len(sorted_list) == 0:
#                 sorted_list.append([ra, dec, flux])
#             elif len(sorted_list) == 1:
#                 if ra <= sorted_list[0][0]:
#                     sorted_list.insert(0,[ra, dec, flux])
#                 else:
#                     sorted_list.append([ra,dec,flux])
#             else:
#                 for i in range(len(sorted_list)-1):
#                     if ra<= sorted_list[i][0]:
#                         sorted_list.insert(i, [ra,dec,flux])
#                         break
#                     elif (ra> sorted_list[i][0] and ra<= sorted_list[-1][0]): 
#                         continue
#                     else:
#                         sorted_list.append([ra,dec,flux])
#                         break


# with open(fr"{folder_loc}diffused_UV_data\RA_sorted_flux_{wavelength}.csv", 'w', newline='') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     # Write the data array to the CSV file row by row
#     for row in sorted_list:
#         csv_writer.writerow(row)

# print(f"RA_sorted_flux_{wavelength}.csv saved")


# with open(fr"{folder_loc}diffused_UV_data\ra_{wavelength}.csv", 'w', newline='') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     # Write the data array to the CSV file row by row
#     for row in ra_grid:
#         csv_writer.writerow(row)
# print(f"ra_{wavelength}.csv saved")
# with open(fr"{folder_loc}diffused_UV_data\dec_{wavelength}.csv", 'w', newline='') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     # Write the data array to the CSV file row by row
#     for row in dec_grid:
#         csv_writer.writerow(row)
# print(f"dec_{wavelength}.csv saved")

# # center is at BQF900

