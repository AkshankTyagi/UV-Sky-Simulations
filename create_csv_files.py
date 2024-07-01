# Create and save photons file + sorted file for all stars and save it
import pandas as pd
import csv
import time
from configparser import ConfigParser

from Params_configparser import *
# from Satellite_configparser import *
from Params_configparser import get_folder_loc
from dust_scatter import get_cords_from_ra_dec
from star_spectrum import *
from star_data import *
from read_dust_files import *

folder_loc = get_folder_loc()
params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global hipp_file, castelli_file
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    min_lim = float(config.get(param_set, 'limit_min'))
    max_lim = float(config.get(param_set, 'limit_max'))
    star_mag_threshold = float(config.get(param_set, 'star_mag_threshold'))
    
    return min_lim, max_lim

wave_min, wave_max = read_parameter_file()

def Get_hipstars():
    STARS = read_hipparcos_data(hipp_file)
    all_spectra = READ_CASTELLI_SPECTRA(castelli_file)
    hipstars = []
    for j in range(len(STARS['hip'])):
        star = stars()
        star.hip_no = STARS['hip'].iloc[j]
        star.index = j
        star.mag = STARS['mag'].iloc[j]
        star.sp_type = STARS['Spectral_type'].iloc[j]
        star.parallax = STARS['trig_parallax'].iloc[j]
        if star.parallax < 0:
            print(f"{j}) hip_no={star.hip_no} excluded ---- parallax is negative")
            continue
        star.ra = STARS['ra_deg'].iloc[j]
        star.dec = STARS['de_deg'].iloc[j]
        star.ebv = STARS['B-V'].iloc[j]
        if pd.isna(star.ebv):
            if pd.isna(STARS['B_mag'].iloc[j]) or pd.isna(STARS['V_mag'].iloc[j]):
                print(f"{j}) hip_no={star.hip_no} excluded ----- does not have B_mag, V_mag, or even B-V ")
                continue
            else:
                star.ebv = STARS['B_mag'].iloc[j] - STARS['V_mag'].iloc[j]
        star.distance= GET_DISTANCE(star.parallax)
        star.x, star.y, star.z = get_cords_from_ra_dec(star.ra, star.dec, star.distance )
        star.temperature= GET_STAR_TEMP(star.sp_type, star.hip_no)
        data = all_spectra[star.temperature]
        wavelengths = data['wavelength']
        Iflux = data['spectrum']
        star.scale, photons = GET_SCALE_FACTOR(star.mag, star.parallax, star.ebv, wavelengths, Iflux )
        star.wavelengths, star.Iflux, star.tot_photons = Trim_Spectral_data(data, photons)
        # print(star.scale, list(zip(star.wavelengths, star.Iflux, star.tot_photons)))
        # print(star.hip_no, list(zip(star.wavelengths,star.tot_photons)))
        star.min_phi = 0
        star.max_phi = 0
        star.min_theta = 0
        star.max_theta = 0
        hipstars.append(star)
    del STARS, data
    print("hipstars created")
    return hipstars

def Create_Allstars_flux(hipstars):
    data_dict = {'Wavelengths': hipstars[0].wavelengths}
    # Populate the dictionary with HIP numbers as keys and total photons as values
    for star in hipstars:
        data_dict[star.index] = star.tot_photons
    # Convert data to ensure it is little-endian
    for key in data_dict:
        if isinstance(data_dict[key], np.ndarray) and data_dict[key].dtype.byteorder == '>':
            data_dict[key] = data_dict[key].byteswap().newbyteorder()

    # Create a DataFrame from the dictionary and save it as csv file
    df = pd.DataFrame(data_dict)
    df.to_csv(f'diffused_data\\Allstars_flux_data[{wave_min},{wave_max}].csv', index=False)
    print(f'\ndiffused_data\\Allstars_flux_data[{wave_min},{wave_max}].csv created')
    return df

def sort_star_list(wavel_range, photon_data):
    num_stars = photon_data.shape[1]

    header = list(range(num_stars))
    header[0] ='wavelengths'
    df = pd.DataFrame(columns= header)

    for i in range(len(wavel_range)):
        print('sorting for ' ,wavel_range[i], i+1,len(wavel_range))
        photon_val = photon_data.iloc[i][1:]
        sorted_row = photon_val.sort_values(ascending =False)
        sorted_header = sorted_row.index.tolist()
        sorted_header.insert(0, wavel_range[i])
        sorted_header_df= pd.DataFrame([sorted_header], columns= df.columns)
        df = pd.concat([df, sorted_header_df], ignore_index = True)
    df.to_csv(f'diffused_data\\Sorted_star_list[{wave_min},{wave_max}].csv', index =False)
    print(f'\ndiffused_data\\Sorted_star_list[{wave_min},{wave_max}].csv created')
    return df

def create_sorted_list(wavelengths):
    try:
        photon_data = pd.read_csv(f'diffused_data\\Allstars_flux_data[{wave_min},{wave_max}].csv')
        # print('1\n',photon_data.iloc[366])
    except FileNotFoundError:
        print(f'\ndiffused_data\\Allstars_flux_data[{wave_min},{wave_max}].csv not_found')
        photon_data = Create_Allstars_flux(hipstars)
    df = sort_star_list(wavelengths, photon_data )
    return df, photon_data

# Calculate weighted probabilities for selecting stars
def calc_weigthed_probab(wavelengths):
    wave_min, wave_max = read_parameter_file()
    try: 
        weights = pd.read_csv(f'diffused_data\\Weighted_list[{wave_min},{wave_max}].csv')
        sorted_list = pd.read_csv(f'diffused_data\\Sorted_star_list[{wave_min},{wave_max}].csv')
        return weights, sorted_list
    except FileNotFoundError:
        try:
            sorted_list = pd.read_csv(f'diffused_data\\Sorted_star_list[{wave_min},{wave_max}].csv')
            photon_data = pd.read_csv(f'diffused_data\\Allstars_flux_data[{wave_min},{wave_max}].csv')
        except FileNotFoundError:
            print('sorted_list, photon_data  csv files not found')
            sorted_list, photon_data = create_sorted_list(wavelengths)
    
    print('sorted_list, photon_data obtained, creating weighted List\n')
    # print(photon_data[:3], sorted_list[:3])

    num_columns = photon_data.shape[1]
    
    header = list(range(num_columns))
    header[0] ='wavelengths'
    df = pd.DataFrame(columns= header)
    # print(df)
    df.to_csv(f'diffused_data\\Weighted_list[{wave_min},{wave_max}].csv', index =False)

    print(num_columns)
    for w in range(len(wavelengths)):
        time1 = int(time.time())
        star_wgt = [0] * num_columns
        max_phot = photon_data.iloc[w][int(sorted_list.iloc[w][1])]
        min_phot_limit = max_phot/1e14

        print("wavelength_Num",w, wavelengths[w], sorted_list.iloc[w][1], max_phot, 'at time:', time1)
        star_wgt[0] = wavelengths[w]
        star_wgt[1] = photon_data.iloc[w][int(sorted_list.iloc[w][1])]

        for i in range(2, num_columns):
            # print(f'{i}) {sorted_list.iloc[w][i]}, photons: {photon_data.iloc[w][int(sorted_list.iloc[w][i])]},  ')
            star_photon = photon_data.iloc[w][int(sorted_list.iloc[w][i])]

            if i%1000 ==2: #Checkpoint for time
                time2 = int(time.time())
                print('Checkpoints------- ',i,') star_no:', sorted_list.iloc[w][i], star_photon, '### duration:', time2 - time1)

            if star_photon < min_phot_limit:
                time3 = int(time.time())
                cumul_phot = star_wgt[i-1]
                print(i, "additions were done in time", time3-time1, "(s),---- Cumul_Photons:", cumul_phot)
                star_wgt[i:num_columns] = [cumul_phot] * (num_columns - i)
                break

            star_wgt[i] = star_wgt[i-1] + star_photon

        with open(f'diffused_data\\Weighted_list[{wave_min},{wave_max}].csv', mode='a') as file:
            writer = csv.writer(file)
            writer.writerow(star_wgt)

    print(f'diffused_data\\Weighted_list[{wave_min},{wave_max}].csv created')
    return pd.read_csv(f'diffused_data\\Weighted_list[{wave_min},{wave_max}].csv'), sorted_list

# hipstars = Get_hipstars()
# wavelengths = hipstars[0].wavelengths
# star_wgt, sort_list = calc_weigthed_probab()


def CHECKPOINT(dust_arr, inp_par, nphoton, tot_star, wcs, hipstars, starlog, misslog, totlog, distlog, scatlog):
    # Time related
    current_time = time()
    print(f"Checkpoint of {nphoton} photons at {ctime(current_time)}")

    # Add to cumulative grids. Scale by the number of photons from the star over
    # the number of photons in the simulation. Write them out to FITS files.

    write_fits_file(wcs, dust_arr, nphoton, tot_star, inp_par)

    with open("datalogger.txt", "w") as logfile:
        logfile.write("HIP_NO Dist.  Star_flux Star_phot Miss_phot Dist_sca scat_flux tot_flux\n")
        for i in range(len(hipstars)):
            logfile.write(f"{hipstars[i]['HIP_NO']} {hipstars[i]['distance']} {hipstars[i]['tot_photons']} "
                          f"{starlog[i]} {misslog[i]} {distlog[i]} {scatlog[i]} {totlog[i]}\n")

def write_fits_file(wcs_out, grid, nphoton, tot_star, inp_par):
    filename = "scattered.fits"
    axes = [wcs_out.nx, wcs_out.ny]
    dust_out = np.zeros((wcs_out.ny, wcs_out.nx), dtype=np.float32)

    for i in range(wcs_out.nx):
        for j in range(wcs_out.ny):
            ipixel = i + j * wcs_out.nx
            dust_out[j, i] = grid[ipixel] * tot_star / nphoton

    hdu = fits.PrimaryHDU(dust_out)
    hdu.header["CRVAL1"] = wcs_out.xrefval
    hdu.header["CRPIX1"] = wcs_out.xrefpix
    hdu.header["CDELT1"] = wcs_out.xinc
    hdu.header["CROTA1"] = wcs_out.rot
    hdu.header["CTYPE1"] = f"GLON{wcs_out.coordtype}"
    hdu.header["CRVAL2"] = wcs_out.yrefval
    hdu.header["CRPIX2"] = wcs_out.yrefpix
    hdu.header["CDELT2"] = wcs_out.yinc
    hdu.header["CROTA2"] = wcs_out.rot
    hdu.header["CTYPE2"] = f"GLAT{wcs_out.coordtype}"
    hdu.header["DATAMIN"] = np.min(dust_out)
    hdu.header["DATAMAX"] = np.max(dust_out)
    hdu.header["NPHOT"] = nphoton
    hdu.header["ALBEDO"] = inp_par.albedo
    hdu.header["G"] = inp_par.g
    hdu.header["WAVELENG"] = inp_par.wave
    hdu.header["COMMENT"] = f"Dust file: {inp_par.dust_file[:30]}"
    hdu.header["COMMENT"] = inp_par.version

    hdu.writeto(filename, overwrite=True)

############


   

