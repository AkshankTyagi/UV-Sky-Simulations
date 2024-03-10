from configparser import ConfigParser
from view_orbit import get_folder_loc
folder_loc = get_folder_loc()

config = ConfigParser()

config['OBJECT A'] = {
    'sat_name' : 'OBJECT A',
    'line1' : '1 56308U 23057A   23125.81960475  .00000660  00000+0  00000+0 0  9992',
    'line2' : '2 56308   9.8856 212.7161 0023120  32.8725 327.3113 14.88280559  2007'
}
config['RISAT-2B'] = {
    'sat_name' : 'RISAT-2B',
    'line1' : '1 44233U 19028A   23134.10174731  .00006152  00000+0  50662-3 0  9998',
    'line2' : '2 44233  36.9998 199.2877 0014452 250.9853 108.9323 14.98747358218208'
}
config['ISS'] = {
    'sat_name' : 'ISS',
    'line1' : '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991',
    'line2' : '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482'
}
config['Astrosat'] = {
    'sat_name' : 'Astrosat',
    'line1' : '1 40930U 15052A   24006.16340880  .00001943  00000-0  17423-3 0  9998',
    'line2' : '2 40930   5.9997  85.2473 0008613 282.4744  77.4420 14.78268977447249'
}

with open(f'{folder_loc}Satellite_TLE.txt',"w") as f:
    config.write(f)