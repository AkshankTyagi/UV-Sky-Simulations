from configparser import ConfigParser
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import time
from Params_configparser import get_folder_loc
from Coordinates import *
# from diffuse_scattering import pc_to_cm

class Dust_params:
    def __init__(self):
        self.dust_xsize = 0
        self.dust_ysize = 0
        self.dust_zsize = 0
        self.sun_x = 0.0
        self.sun_y = 0.0
        self.sun_z = 0.0
        self.dust_binsize = 0.0
        self.num_photon = 0
        self.nscatter = 0
        self.wave = 0.0
        self.albedo = 0.0
        self.g = 0.0
        self.version = ""
        self.print_debug = ""
        self.min_gl_debug = 0.0
        self.max_gl_debug = 0.0
        self.min_gb_debug = 0.0
        self.max_gb_debug = 0.0

class stars:
    def __init__(self):
        self.hip_no = 0
        self.index =0
        self.mag = 0
        self.sp_type = 0
        self.parallax = 0
        self.ebv = 0
        self.ra = 0
        self.dec = 0
        self.gl = 0
        self.gb = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.distance = 0
        self.temperature = 0
        self.wavelengths = []
        self.Iflux = []
        # self.nstars = nstars
        self.scale = 0
        self.tot_photons = []
        self.min_theta = 0
        self.max_theta = 0
        self.min_phi = 0
        self.max_phi = 0

folder_loc, params_file = get_folder_loc()


def pc_to_cm(x):
    return 3.08568024696E+18 * x

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)

    global dust_file, dust_col_file, sigma_file, hipp_file, castelli_file, star_mag_threshold
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    dust_file = config.get(param_set, 'dust_file')
    dust_col_file = config.get(param_set, 'dust_col_file')
    sigma_file = config.get(param_set, 'sigma_file')
    min_lim = float(config.get(param_set, 'limit_min'))
    max_lim = float(config.get(param_set, 'limit_max'))
    star_mag_threshold = float(config.get(param_set, 'star_mag_threshold'))

    return min_lim, max_lim

def read_scatter_parameters(filename = params_file, param_set = 'Scatter_params'):
    dust_par = Dust_params()  # Create a Dust_params class
    config = ConfigParser()
    config.read(filename)

    # Dust parameters: num_photons, num_scatter, wavelength_list, albedo, phase_func
    dust_par.num_photon = float(config.get(param_set, 'No_photons'))
    dust_par.nscatter = float(config.get(param_set, 'No_scatter'))
    wavelength = config.get(param_set, 'wavelength').strip('[]')
    dust_par.wave = [float(value) for value in wavelength.split(',')]
    dust_par.albedo = float(config.get(param_set, 'Albedo'))
    dust_par.g = float(config.get(param_set, 'Phase_func'))

    # Default debugging parameters
    dust_par.print_debug = config.get(param_set, 'print_debug')
    dust_par.min_gl_debug = float(config.get(param_set, 'min_gl_debug'))
    dust_par.max_gl_debug = float(config.get(param_set, 'max_gl_debug'))
    dust_par.min_gb_debug = float(config.get(param_set, 'min_gb_debug'))
    dust_par.max_gb_debug = float(config.get(param_set, 'max_gb_debug'))

    return dust_par

def read_wcs_parameters(filename = params_file, param_set = 'WCS'):
    config = ConfigParser()
    config.read(filename)

    crval1 = float(config[param_set]['CRVAL1'])
    crval2 = float(config[param_set]['CRVAL2'])
    crpix1 = float(config[param_set]['CRPIX1'])
    crpix2 = float(config[param_set]['CRPIX2'])
    cdelt1 = float(config[param_set]['CDELT1'])
    cdelt2 = float(config[param_set]['CDELT2'])
    crot = float(config[param_set]['CROTA'])
    ctype = config[param_set]['CTYPE']
    naxis1 = int(config[param_set]['NAXIS1'])
    naxis2 = int(config[param_set]['NAXIS2'])

    wcs_params = [
        crval1, crval2,
        crpix1, crpix2,
        cdelt1, cdelt2,
        crot, ctype,
        naxis1, naxis2
    ]
    return wcs_params

def dust_read(dust_par, filename):
    with fits.open(filename) as hdul:
        data = hdul[0].data
        naxes = data.shape
        dust_arr = np.array(data, dtype=float)
        dust_par.dust_xsize = float(naxes[0])
        dust_par.dust_ysize = float(naxes[1])
        dust_par.dust_zsize = float(naxes[2])
        dust_par.dust_binsize = float(hdul[0].header['CRDELT1'])
        dust_par.sun_x = float(hdul[0].header['CRPIX1'])
        dust_par.sun_y = float(hdul[0].header['CRPIX2'])
        dust_par.sun_z = float(hdul[0].header['CRPIX3'])

        # print(f"dust_par.dust_xsize:{dust_par.dust_xsize}, dust_par.dust_binsize :{dust_par.dust_binsize }")
        # print(f"dust_par.sun_x:{dust_par.sun_x}\ndust_par.sun_y: {dust_par.sun_y}, \ndust_par.sun_z : {dust_par.sun_z}")

    return dust_arr

def read_cross_sec(sigma_file, dust_par):
    wavelengths = []
    crossX = []
    start_line_found = False

    with open(sigma_file, "r") as fp:
        for line in fp:
            if not start_line_found:
                if line.startswith("-"):
                    start_line_found = True
                continue
            if line.startswith("-"):
                break
            parts = line.split()
            if len(parts) >= 4:
                wave = float(parts[0])
                cross_sec = float(parts[3])
                wavelengths.append(wave * 10000)  # Convert to Ångstroms
                crossX.append(cross_sec)
    

    wavelengths = np.array(wavelengths)
    crossX = np.array(crossX)
    wavel = np.array(dust_par.wave, dtype=float)
    sigma = []

    for wave in wavel:
        index = np.argmin(np.abs(wavelengths - wave))
        sigma_value = crossX[index] * pc_to_cm(1)  # Convert to the desired units
        sigma.append(sigma_value)

    # print(f'wavel: {wavel}')
    return sigma



def write_to_grid(wcs, ra, dec, flux, grid):
 # shape of the grid

    # Perform world to pixel transformation
    gl, gb = conv_eq_to_gal(ra, dec)
    xout, yout = wcs.world_to_pixel_values(gl, gb)
    # xout, yout = wcs.world_to_pixel_values(ra, dec)
    # print(f"ra:{ra}, dec:{dec}, gl:{gl}, gb:{gb}, xout:{xout}, yout:{yout}, xout2:{xout2}, yout:{yout2}" )

    # Check if the pixel coordinates are within the bounds of the grid
    if 1 <= xout <= wcs.array_shape[1] and 1 <= yout <= wcs.array_shape[0]:
        # Correct because FITS standard is to start at 1
        xout -= 1
        yout -= 1

        # print(xout, yout, flux / abs(np.deg2rad(wcs.wcs.cdelt[0]) * np.deg2rad(wcs.wcs.cdelt[1])))
        # Add the flux to the grid
        grid[int(yout), int(xout)] = grid[int(yout), int(xout)] + flux / abs(np.deg2rad(wcs.wcs.cdelt[0]) * np.deg2rad(wcs.wcs.cdelt[1]))

    return grid

# def write_to_gridgg(wcs_param, ra, dec, flux, grid2):
#     wcs = WCS(naxis=2)
#     wcs.wcs.crval = [wcs_param[0], wcs_param[1]]
#     wcs.wcs.crpix = [wcs_param[2], wcs_param[3]]
#     wcs.wcs.cdelt = [wcs_param[4], wcs_param[5]]
#     wcs.wcs.crota = [wcs_param[6], wcs_param[6]]
#     wcs.wcs.ctype = [f"RA--{wcs_param[7]}", f"DEC-{wcs_param[7]}"]
#     wcs.array_shape = [ wcs_param[9], wcs_param[8],]  # shape of the grid

#     # Perform world to pixel transformation
#     # gl, gb = conv_eq_to_gal(ra, dec)
#     xout, yout = wcs.world_to_pixel_values(ra, dec)
#     # print(f"ra{ra}, dec{dec}, gl{gl}, gb{gb}, xout{xout}, yout{yout}" )

#     # Check if the pixel coordinates are within the bounds of the grid
#     if 1 <= xout <= wcs.array_shape[1] and 1 <= yout <= wcs.array_shape[0]:
#         # Correct because FITS standard is to start at 1
#         xout -= 1
#         yout -= 1

#         # print(xout, yout, ipixel)
#         # Add the flux to the grid
#         grid2[int(yout), int(xout)] += flux / abs(np.deg2rad(wcs.wcs.cdelt[0]) * np.deg2rad(wcs.wcs.cdelt[1]))

#     return grid2

def CHECKPOINT(dust_arr, inp_par, nphoton, tot_star, wcs, hipstars, wavelength):#, starlog, misslog, totlog, distlog, scatlog):
    # Time related
    current_time = time.time()
    print(f"Checkpoint of {nphoton} for {wavelength} photons at {current_time}")

    # Add to cumulative grids. Scale by the number of photons from the star over
    # the number of photons in the simulation. Write them out to FITS files.

    write_fits_file(wcs, dust_arr, nphoton, tot_star, inp_par, wavelength)

    # with open("datalogger.txt", "w") as logfile:
        # logfile.write("HIP_NO Dist.  Star_flux Star_phot Miss_phot Dist_sca scat_flux tot_flux\n")
        # for i in range(len(hipstars)):
        #     logfile.write(f"{hipstars[i].hip_no} {hipstars[i].distance} {hipstars[i].tot_photons} "
        #                   f"{starlog[0][i]} {misslog[0][i]} {distlog[0][i]} {scatlog[0][i]} {totlog[0][i]}\n")


def write_fits_file(wcs_param, grid, nphoton, tot_star, inp_par, wavelength):#i, wavelengths_list):
    min_wave, max_wave = read_parameter_file()
    filename = f"diffused_output{os.sep}scattered_{int(inp_par.num_photon)}[{int(wavelength)}]_mag{int(star_mag_threshold)}.fits"

    dust_out = (np.round(grid * tot_star / inp_par.num_photon, decimals=3)).astype(np.float32)
    # dust_out = (grid * tot_star / inp_par.num_photon).astype(np.int32)

    # dust_out[1500, 500] = 1000
    # plot_diffused_bg(dust_out, wavelengths_list[i])

    # Create or open the FITS file
    # if i == 0:
        # Create a new FITS file with the primary HDU
    primary_hdu = fits.PrimaryHDU(dust_out)
    primary_hdu.header["CRVAL1"] = wcs_param[0]
    primary_hdu.header["CRPIX1"] = wcs_param[2]
    primary_hdu.header["CDELT1"] = wcs_param[4]
    primary_hdu.header["CROTA1"] = wcs_param[6]
    primary_hdu.header["CTYPE1"] = f"GLON{wcs_param[7]}"
    primary_hdu.header["CRVAL2"] = wcs_param[1]
    primary_hdu.header["CRPIX2"] = wcs_param[3]
    primary_hdu.header["CDELT2"] = wcs_param[5]
    primary_hdu.header["CROTA2"] = wcs_param[6]
    primary_hdu.header["CTYPE2"] = f"GLAT{wcs_param[7]}"
    primary_hdu.header["DATAMIN"] = np.min(dust_out)
    primary_hdu.header["DATAMAX"] = np.max(dust_out)
    primary_hdu.header["NPHOT"] = inp_par.num_photon
    primary_hdu.header["ALBEDO"] = inp_par.albedo
    primary_hdu.header["G"] = inp_par.g
    primary_hdu.header["WAVELENG"] = wavelength
    
    hdulist = fits.HDUList([primary_hdu])
    hdulist.writeto(filename, overwrite=True)

    # else:
    #     # Open the existing FITS file
    #     with fits.open(filename, mode='update' ) as hdulist:

    #         image_hdu = fits.ImageHDU(dust_out)
    #         image_hdu.header["CRVAL1"] = wcs_param[0]
    #         image_hdu.header["CRPIX1"] = wcs_param[2]
    #         image_hdu.header["CDELT1"] = wcs_param[4]
    #         image_hdu.header["CROTA1"] = wcs_param[6]
    #         image_hdu.header["CTYPE1"] = f"GLON{wcs_param[7]}"
    #         image_hdu.header["CRVAL2"] = wcs_param[1]
    #         image_hdu.header["CRPIX2"] = wcs_param[3]
    #         image_hdu.header["CDELT2"] = wcs_param[5]
    #         image_hdu.header["CROTA2"] = wcs_param[6]
    #         image_hdu.header["CTYPE2"] = f"GLAT{wcs_param[7]}"
    #         image_hdu.header["DATAMIN"] = np.min(dust_out)
    #         image_hdu.header["DATAMAX"] = np.max(dust_out)
    #         image_hdu.header["NPHOT"] = inp_par.num_photon
    #         image_hdu.header["ALBEDO"] = inp_par.albedo
    #         image_hdu.header["G"] = inp_par.g
    #         image_hdu.header["WAVELENG"] = wavelengths_list[i]
            
    #         hdulist.append(image_hdu)
    #         hdulist.flush()
            # hdulist.writeto(filename, overwrite=True)

    print(f'{filename} updated\n')
    hdulist.close()