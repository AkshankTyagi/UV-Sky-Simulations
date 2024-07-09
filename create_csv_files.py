# Create and save photons file + sorted file for all stars and save it
import pandas as pd
import csv
import time
import ast
import os
from configparser import ConfigParser

from Params_configparser import *
# from Satellite_configparser import *
# from Params_configparser import get_folder_loc
from Coordinates import *
from star_spectrum import *
from star_data import *
from read_dust_files import *

folder_loc, params_file = get_folder_loc()
# params_file = f'{folder_loc}init_parameter.txt'

def read_parameter_file(filename= params_file, param_set = 'Params_1'):
    config = ConfigParser()
    config.read(filename)
    global hipp_file, castelli_file, exclude_stars_file, wave_min, wave_max, star_mag_threshold
    hipp_file = config.get(param_set, 'hipparcos_catalogue')
    castelli_file = config.get(param_set, 'Castelli_data')
    exclude_stars_file = config.get(param_set, 'exclude_stars')
    wave_min = float(config.get(param_set, 'limit_min'))
    wave_max = float(config.get(param_set, 'limit_max'))
    star_mag_threshold = float(config.get(param_set, 'star_mag_threshold'))

    return wave_min, wave_max

_, _ = read_parameter_file()

def read_excluded_stars(file =exclude_stars_file,):
    excluded_stars_list = np.loadtxt(file, dtype=int)
    return excluded_stars_list.tolist()

def Get_hipstars():
    try:
        hipstars = pd.read_csv(f'hipstars_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
        print('hipstars_df obtained')
    except FileNotFoundError:
        print(f'hipstars_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv not found')
        STARS = read_hipparcos_data(hipp_file)
        all_spectra = READ_CASTELLI_SPECTRA(castelli_file)
        excluded_stars = read_excluded_stars()

        star_data_list = []
        k = 0
        
        # Iterate over each star in the STARS DataFrame
        for j in range(len(STARS['hip'])):
            star_data = {}
            hip_no = STARS['hip'].iloc[j]
            
            if hip_no in excluded_stars:  # Trim excluded_stars after removing the star
                excluded_stars.remove(hip_no)
                print(f'hip_no={hip_no} excluded ')
                continue

            # Fill the star_data dictionary with star properties
            star_data['index'] = k
            star_data['hip_no'] = hip_no         
            star_data['mag'] = STARS['mag'].iloc[j]
            star_data['sp_type'] = STARS['Spectral_type'].iloc[j]
            star_data['parallax'] = STARS['trig_parallax'].iloc[j]
            star_data['ebv'] = STARS['B-V'].iloc[j]
            star_data['ra'] = STARS['ra_deg'].iloc[j]
            star_data['dec'] = STARS['de_deg'].iloc[j]
            star_data['distance'] = GET_DISTANCE(star_data['parallax'])
            star_data['x'], star_data['y'], star_data['z'] = conv_eq_to_cart(star_data['ra'], star_data['dec'], star_data['distance'])
            star_data['gl'], star_data['gb'] = conv_eq_to_gal(star_data['ra'], star_data['dec'])
            star_data['temperature'] = GET_STAR_TEMP(star_data['sp_type'], star_data['hip_no'])
            data = all_spectra[star_data['temperature']]
            wavelengths = data['wavelength']
            Iflux = data['spectrum']
            star_data['scale'], photons = GET_SCALE_FACTOR(star_data['mag'], star_data['parallax'], star_data['ebv'], wavelengths, Iflux)
            star_data['wavelengths'], star_data['Iflux'], star_data['tot_photons'] = Trim_Spectral_data(data, photons)
            star_data['min_phi'] = 0
            star_data['max_phi'] = 0
            star_data['min_theta'] = 0
            star_data['max_theta'] = 0
    
            k += 1
            star_data_list.append(star_data)
    
        # Convert the list of dictionaries to a DataFrame
        hipstars_df = pd.DataFrame(star_data_list)

        # Save the DataFrame to a CSV file
        print(excluded_stars)
        hipstars_df.to_csv(f'hipstars_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv', index=False)
        print(f"hipstars DataFrame created and saved to 'hipstars[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv'")
        hipstars = pd.read_csv(f'hipstars_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
    return hipstars

# def Get_hipstars():
#     STARS = read_hipparcos_data(hipp_file)
#     all_spectra = READ_CASTELLI_SPECTRA(castelli_file)
#     hipstars = []
#     excluded_stars = read_excluded_stars()
#     k = 0
#     for j in range(len(STARS['hip'])):
#         star = stars()
#         star.hip_no = STARS['hip'].iloc[j]
#         if star.hip_no in excluded_stars[:3]:  #since i am trimming excluded_stars after removing the star       # [0] or star.hip_no == excluded_stars[1]:
#             excluded_stars.remove(star.hip_no)
#             print(f'hip_no={star.hip_no} excluded ')  #  , next_exclusion{excluded_stars[0]}
#             continue

#         star.parallax = STARS['trig_parallax'].iloc[j]
#         # if star.parallax < 0:
#         #     print(f"{j}) hip_no={star.hip_no} excluded ----- parallax is negative")
#         #     excluded_stars.append(star.hip_no)
#         #     continue

#         star.ebv = STARS['B-V'].iloc[j]
#         # if pd.isna(star.ebv):
#         #     if pd.isna(STARS['B_mag'].iloc[j]) or pd.isna(STARS['V_mag'].iloc[j]):
#         #         print(f"{j}) hip_no={star.hip_no} excluded ----- does not have B_mag, V_mag, or even B-V ")
#         #         excluded_stars.append(star.hip_no)
#         #         continue
#         #     else:
#         #         star.ebv = STARS['B_mag'].iloc[j] - STARS['V_mag'].iloc[j]
        
#         star.index = k
#         star.mag = STARS['mag'].iloc[j]
#         star.sp_type = STARS['Spectral_type'].iloc[j]
#         star.ra = STARS['ra_deg'].iloc[j]
#         star.dec = STARS['de_deg'].iloc[j]
#         star.distance= GET_DISTANCE(star.parallax)
#         star.x, star.y, star.z = conv_eq_to_cart(star.ra, star.dec, star.distance )
#         star.gl, star.gb = conv_eq_to_gal(star.ra, star.dec)
#         star.temperature= GET_STAR_TEMP(star.sp_type, star.hip_no)
#         data = all_spectra[star.temperature]
#         wavelengths = data['wavelength']
#         Iflux = data['spectrum']
#         star.scale, photons = GET_SCALE_FACTOR(star.mag, star.parallax, star.ebv, wavelengths, Iflux )
#         star.wavelengths, star.Iflux, star.tot_photons = Trim_Spectral_data(data, photons)
#         # print(star.scale, list(zip(star.wavelengths, star.Iflux, star.tot_photons)))
#         # print(star.hip_no, list(zip(star.wavelengths,star.tot_photons)))
#         star.min_phi = 0
#         star.max_phi = 0
#         star.min_theta = 0
#         star.max_theta = 0
#         k += 1
#         hipstars.append(star)
#     del STARS, data
#     print("hipstars created")
#     return hipstars , excluded_stars

def Create_Allstars_flux(hipstars):
    wavelengths_str = hipstars.loc[0, 'wavelengths']
    wavelengths_str = wavelengths_str.replace(' ', ',')
    wavelengths_list = ast.literal_eval(wavelengths_str)
    print(wavelengths_list)

    data_dict = {'Wavelengths': wavelengths_list}
    for i, row in hipstars.iterrows():
        phot_str = row['tot_photons']
        # phot_str= phot_str.replace(' ', ',')
        
        # Handle missing or NaN values
        if phot_str[:4] == '[nan':
            phot = [0] * len(wavelengths_list)  # Or handle as you prefer
        else:
            phot = ast.literal_eval(phot_str)

        data_dict[f'{int(row["index"])}'] = phot
    
    # Create the DataFrame from the dictionary
    df = pd.DataFrame(data_dict)
    df = df.round(3)
    print(df.head())

    df.to_csv(f'diffused_data{os.sep}Allstars_flux_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv', index=False)
    print(f'\ndiffused_data{os.sep}Allstars_flux_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv created')
    return df

def sort_star_list(wavel_range, photon_data):
    num_column = photon_data.shape[1]

    header = list(range(num_column))
    header[0] ='wavelengths'
    df = pd.DataFrame(columns= header)

    for i in range(len(wavel_range)):
        print('\nsorting for ' ,wavel_range[i], i+1,len(wavel_range))
        photon_val = photon_data.iloc[i][1:]
        sorted_row = photon_val.sort_values(ascending =False)
        sorted_header = sorted_row.index.tolist()
        sorted_header.insert(0, wavel_range[i])
        sorted_header_df= pd.DataFrame([sorted_header], columns= df.columns)
        df = pd.concat([df, sorted_header_df], ignore_index = True)
    
    print(df.head())
    
    df.to_csv(f'diffused_data{os.sep}Sorted_star_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv', index =False)
    print(f'\ndiffused_data{os.sep}Sorted_star_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv created')
    return df

def create_sorted_list(wavelengths, hipstars):
    try:
        photon_data = pd.read_csv(f'diffused_data{os.sep}Allstars_flux_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
        # print('1\n',photon_data.iloc[366])
    except FileNotFoundError:
        print(f'\ndiffused_data{os.sep}Allstars_flux_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv not_found')
        photon_data = Create_Allstars_flux(hipstars)
    df = sort_star_list(wavelengths, photon_data )
    return df, photon_data

# Calculate weighted probabilities for selecting stars
def calc_weigthed_probab(wavelengths, hipstars):
    # wave_min, wave_max = read_parameter_file()
    try: 
        weights = pd.read_csv(f'diffused_data{os.sep}Weighted_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
        sorted_list = pd.read_csv(f'diffused_data{os.sep}Sorted_star_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
        print('weighted_list, sorted_stars_list obtained.')
        return weights, sorted_list
    except FileNotFoundError:
        try:
            sorted_list = pd.read_csv(f'diffused_data{os.sep}Sorted_star_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
            photon_data = pd.read_csv(f'diffused_data{os.sep}Allstars_flux_data[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv')
        except FileNotFoundError:
            print('sorted_list, photon_data  csv files not found')
            sorted_list, photon_data = create_sorted_list(wavelengths, hipstars)
    
        print('sorted_list, photon_data obtained, creating weighted List\n')
        # print(photon_data[:3], sorted_list[:3])

        num_columns = photon_data.shape[1]
        
        header = list(range(num_columns))
        header[0] ='wavelengths'
        df = pd.DataFrame(columns= header)
        # print(df)
        df.to_csv(f'diffused_data{os.sep}Weighted_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv', index =False)

        print('num_columns:', num_columns)

        for w in range(len(wavelengths)):
            time1 = int(time.time())
            star_wgt = [0] * num_columns
            max_phot = photon_data.iloc[w][int(sorted_list.iloc[w][1])+1]
            min_phot_limit = max_phot/1e10

            print("wavelength_Num",w, wavelengths[w], int(sorted_list.iloc[w][1]), max_phot, 'at time:', time1)
            star_wgt[0] = wavelengths[w]
            star_wgt[1] = photon_data.iloc[w][int(sorted_list.iloc[w][1])+1]

            for i in range(2, num_columns):
                # print(f'{i}) {sorted_list.iloc[w][i]}, photons: {photon_data.iloc[w][int(sorted_list.iloc[w][i])]},  ')
                star_photon = photon_data.iloc[w][int(sorted_list.iloc[w][i])+1]

                if i%1000 ==2: #Checkpoint for time
                    time2 = int(time.time())
                    print('Checkpoints------- ',i,') star_no:', int(sorted_list.iloc[w][i]), star_photon, '### duration:', time2 - time1)

                if star_photon < min_phot_limit or i > 3003:
                    time3 = int(time.time())
                    cumul_phot = star_wgt[i-1]
                    print(i, "additions were done in time", time3-time1, "(s),---- Cumul_Photons:", cumul_phot)
                    star_wgt[i:num_columns] = [cumul_phot] * (num_columns - i)
                    break

                star_wgt[i] = star_wgt[i-1] + star_photon

            with open(f'diffused_data{os.sep}Weighted_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv', mode='a') as file:
                writer = csv.writer(file)
                writer.writerow(star_wgt)

    print(f'\nWeighted_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv created')
    return pd.read_csv(f'diffused_data{os.sep}Weighted_list[{int(wave_min)},{int(wave_max)}]_mag{int(star_mag_threshold)}.csv'), sorted_list


#--------
# Example RUn

# hipstars, _ = Get_hipstars()
# wavelengths = hipstars[0].wavelengths
# star_wgt, sort_list = calc_weigthed_probab(wavelengths)

############


   

