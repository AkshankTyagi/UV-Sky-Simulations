import numpy as np
import math
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from numba import njit

# from dust_extinction.grain_models import D03

NANGLE =  100000

@njit
def PHASE_FUNCTION(cangle, g):
    # Henyey-Greenstein phase function
    pscat = (1 - g **2) / (4. * math.pi* (pow(1. + g**2 - 2. * g * cangle, 1.5)))
    return pscat

@njit
def SETUP_ANGLE():
    angle = np.zeros(NANGLE) 
    # NANGLE = len(angle)  # Assuming NANGLE is defined outside the function
    angle[0] = 0.0
    for i in range(1, NANGLE):
        theta = float(i) / float(NANGLE) * math.pi
        angle[i] = angle[i-1] + math.sin(theta)
    
    # Normalize the angles
    max_angle = angle[NANGLE - 1]
    for i in range(1, NANGLE):
        angle[i] /= max_angle
    
    return angle

@njit
def SETUP_SCATTER( g):
    fscat = np.zeros(NANGLE)  # Assuming NANGLE is defined outside the function
    fscat[0] = 0.0
    
    for i in range(1, NANGLE):
        theta = float(i) / float(NANGLE) * math.pi
        fscat[i] = fscat[i-1] + PHASE_FUNCTION(math.cos(theta), g)
    
    # Normalize the scattering probabilities
    max_fscat = fscat[NANGLE - 1]
    for i in range(1, NANGLE):
        fscat[i] /= max_fscat

    return fscat

@njit
def CALC_THETA(angle, random):
    NANGLE = len(angle)
    
    # Binary Search
    min_i = 0
    max_i = NANGLE
    while max_i - min_i > 1:
        mid_i = (min_i + max_i) // 2
        if angle[mid_i] < random:
            min_i = mid_i
        else:
            max_i = mid_i

    i = min_i
    while i < NANGLE and angle[i] < random:
        i += 1

    if i == 0:
        index = 0.0
    else:
        index = float(i - 1) + (random - angle[i - 1]) / (angle[i] - angle[i - 1])

    theta = index / float(NANGLE) * math.pi  # Theta will be between 0 and PI in radians
    return theta

# @njit
def DETECT(x1, y1, z1, dx, dy, dz, dust):
    # Calculate distances
    d1 = math.sqrt(dx*dx + dy*dy + dz*dz)
    x2 = dust.sun_x - x1
    y2 = dust.sun_y - y1
    z2 = dust.sun_z - z1
    d2 = math.sqrt(x2*x2 + y2*y2 + z2*z2)

    # Calculate the angle (cosine)
    cangl = (dx*x2 + dy*y2 + dz*z2) / d1 / d2
    # print('cangle:', cangl)

    # Calculate the phase function
    pscat = PHASE_FUNCTION(cangl, dust.g)
    return pscat

# @njit
def SCATTER(delta_x, delta_y, delta_z, theta, phi):
    xnew = math.cos(phi) * math.sin(theta)
    ynew = math.sin(phi) * math.sin(theta)
    znew = math.cos(theta)

    dist = math.sqrt(delta_x**2 + delta_y**2 + delta_z**2)
    dxy = math.sqrt(delta_x**2 + delta_y**2)
    cbeta = delta_z / dist
    sbeta = dxy / dist
    calpha = delta_x / dxy
    salpha = delta_y / dxy

    rot_mat = [
        [calpha*cbeta, cbeta*salpha, -sbeta],
        [-salpha, calpha, 0],
        [calpha*sbeta, sbeta*salpha, cbeta]
    ]

    delta_x = rot_mat[0][0]*xnew + rot_mat[0][1]*ynew + rot_mat[0][2]*znew
    delta_y = rot_mat[1][0]*xnew + rot_mat[1][1]*ynew + rot_mat[1][2]*znew
    delta_z = rot_mat[2][0]*xnew + rot_mat[2][1]*ynew + rot_mat[2][2]*znew

    return delta_x, delta_y, delta_z

def CHECK_LIMITS(x, y, z, dust):
    if 0 < x < (dust.dust_xsize - 1) and \
       0 < y < (dust.dust_ysize - 1) and \
       0 < z < (dust.dust_zsize - 1):
        return True
    else:
        return False

def INCREMENT(x_n, y_n, z_n, delta_x, delta_y, delta_z):

    # if x_n[0] > 1e154 or y_n[0] > 1e154 or z_n[0] > 1e154 or delta_y>1e154 or delta_z>1e154 or delta_x>1e154:
    #     print("input has problem",x_n, y_n, z_n, delta_x, delta_y, delta_z)
    x_new = x_n[0] + delta_x
    y_new = y_n[0] + delta_y
    z_new = z_n[0] + delta_z
    # if x_new > 1e154 or y_new > 1e154 or z_new > 1e154:
    #     print("whattt!!!, heres the poblem",x_new,y_new,z_new)
    return x_new, y_new, z_new

def CALC_DELTA_X( theta, phi):
    delta_x = math.cos(phi) * math.sin(theta)
    delta_y = math.sin(phi) * math.sin(theta)
    delta_z = math.cos(theta)

    return delta_x, delta_y, delta_z

def GET_DUST_INDEX(x_new, y_new, z_new):
    dust_index = [int(x_new + 0.5), int(y_new + 0.5), int(z_new + 0.5)]
    return dust_index

@njit
def CALC_DIST(a, b, c, x, y, z):
    # dx = np.abs(a - x)
    # dy = np.abs(b - y)
    # dz = np.abs(c - z)

    # Ensure values are within a safe range before squaring

    # dist = dx**2 + dy**2 + dz**2

    return (a-x)**2 + (b - y)**2 + (c - z)**2

def FIRST_PHOTON(dust_par, star, angle):#,  ran_array, nrandom, ran_ctr, init_seed):
    NEW_PHOTON = True
    iter = 0
    step_delta = 0
    # print("NEW_PHOTON started")
    while NEW_PHOTON == True and iter <= dust_par.num_photon:
        iter += 1
        xp = star['x'] / dust_par.dust_binsize + dust_par.sun_x
        yp = star['y'] / dust_par.dust_binsize + dust_par.sun_y
        zp = star['z'] / dust_par.dust_binsize + dust_par.sun_z
        x_new = xp
        y_new = yp
        z_new = zp
        # if x_new > 1e154 or y_new > 1e154 or z_new > 1e154:
        #     print("1st calc_dist has problem",x_new,y_new,z_new)
        dp = CALC_DIST(xp, yp, zp, dust_par.sun_x, dust_par.sun_y, dust_par.sun_z)
        d_new = dp
        d_old = dp
        
        unif_rand = GET_RANDOM_ARRAY()# ran_array, nrandom, ran_ctr, init_seed)
        theta = CALC_THETA(angle, unif_rand)
        phi = GET_RANDOM_ARRAY()* 2 * math.pi# ran_array, nrandom, ran_ctr, init_seed) 

        delta_x_ptr, delta_y_ptr, delta_z_ptr = CALC_DELTA_X(theta, phi)
        # print("NEW_PHOTON working", iter)

        dtheta = math.degrees(theta)
        dphi = math.degrees(phi)

        # print(f"dtheta:{dtheta}, dphi:{dphi}, {star['min_theta'], star['max_theta'], star['min_phi'], star['max_phi']}   ")
        if dtheta < star['min_theta'] or dtheta > star['max_theta'] or dphi < star['min_phi'] or dphi > star['max_phi']:
            return x_new, y_new, z_new, delta_x_ptr, delta_y_ptr, delta_z_ptr, -1

        if not CHECK_LIMITS(x_new, y_new, z_new, dust_par):
            step_delta = 1000
            while step_delta >= 1:
                while not (CHECK_LIMITS(x_new, y_new, z_new, dust_par)) and d_new <= d_old:
                    delta_x_ptr *= step_delta
                    delta_y_ptr *= step_delta
                    delta_z_ptr *= step_delta
                    # print(delta_x_ptr, delta_y_ptr, delta_z_ptr) 
                    x_new, y_new, z_new = INCREMENT([x_new], [y_new], [z_new], delta_x_ptr, delta_y_ptr, delta_z_ptr)
                    d_old = d_new
                    d_new = CALC_DIST(x_new, y_new, z_new, dust_par.sun_x, dust_par.sun_y, dust_par.sun_z)
                
                x_new, y_new, z_new = INCREMENT([x_new], [y_new], [z_new], -delta_x_ptr, -delta_y_ptr, -delta_z_ptr)
                
                dp = math.sqrt(delta_x_ptr**2 + delta_y_ptr**2 + delta_z_ptr**2)
                delta_x_ptr /= dp
                delta_y_ptr /= dp
                delta_z_ptr /= dp
                d_new = CALC_DIST(x_new, y_new, z_new, dust_par.sun_x, dust_par.sun_y, dust_par.sun_z)
                d_old = d_new
                step_delta /= 2
            
            x_new, y_new, z_new = INCREMENT([x_new], [y_new], [z_new], delta_x_ptr, delta_y_ptr, delta_z_ptr)


        if not CHECK_LIMITS(x_new, y_new, z_new, dust_par):
            iter = -1

        NEW_PHOTON = False

    return x_new, y_new, z_new, delta_x_ptr, delta_y_ptr, delta_z_ptr, iter

# def MATCH_TAU(x_new, y_new, z_new, dust_par, delta_x_ptr, delta_y_ptr, delta_z_ptr, sigma, dust_dist, tau):

#     delta_x = delta_x_ptr
#     delta_y = delta_y_ptr
#     delta_z = delta_z_ptr
#     cum_tau = 0
#     step_size = 1.0  # We travel around in fractions of a bin
#     cont = True
    
#     while cont:
#         # print(f'tau:{tau},-- {CHECK_LIMITS(x_new + delta_x, y_new + delta_y, z_new + delta_z, dust_par)}')
#         while cum_tau <= tau and CHECK_LIMITS(x_new + delta_x, y_new + delta_y, z_new + delta_z, dust_par):
#             x_new, y_new, z_new = INCREMENT([x_new], [y_new], [z_new], delta_x, delta_y, delta_z)
#             # print(x_new, y_new, z_new)
#             dust_index = GET_DUST_INDEX(x_new, y_new, z_new)
#             # print(f'cum_tau= {cum_tau}, dust_val = {dust_dist[dust_index[0], dust_index[1], dust_index[2]]}, sigma:{sigma}, step_size:{step_size} ')
#             cum_tau += (dust_dist[dust_index[0], dust_index[1], dust_index[2]] * sigma[0] * step_size * dust_par.dust_binsize)
        
#         if step_size < 0.1:
#             cont = False
        
#         if cont and CHECK_LIMITS(x_new, y_new, z_new, dust_par):
#             dust_index = GET_DUST_INDEX(x_new, y_new, z_new)
#             cum_tau -= (dust_dist[dust_index[0], dust_index[1], dust_index[2]] * sigma[0] * step_size * dust_par.dust_binsize)
#             x_new, y_new, z_new = INCREMENT([x_new], [y_new], [z_new], -delta_x, -delta_y, -delta_z)
#             delta_x /= 10.0
#             delta_y /= 10.0
#             delta_z /= 10.0 
        
#         step_size /= 10.0
    
#     return cum_tau

def MATCH_TAU(x_new, y_new, z_new, dust_par, delta_x_ptr, delta_y_ptr, delta_z_ptr, sigma, dust_dist, tau):
    delta_x = delta_x_ptr
    delta_y = delta_y_ptr
    delta_z = delta_z_ptr
    cum_tau = 0
    step_size = 1.0
    cont = True
    
    while cont:
        while cum_tau <= tau and CHECK_LIMITS(x_new + delta_x, y_new + delta_y, z_new + delta_z, dust_par):
            x_new += delta_x
            y_new += delta_y
            z_new += delta_z
            dust_index = GET_DUST_INDEX(x_new, y_new, z_new)
            cum_tau += (dust_dist[dust_index[0], dust_index[1], dust_index[2]] * sigma * step_size * dust_par.dust_binsize)
        
        if step_size < 0.1:
            cont = False
        
        if cont and CHECK_LIMITS(x_new, y_new, z_new, dust_par):
            dust_index = GET_DUST_INDEX(x_new, y_new, z_new)
            cum_tau -= (dust_dist[dust_index[0], dust_index[1], dust_index[2]] * sigma * step_size * dust_par.dust_binsize)
            x_new -= delta_x
            y_new -= delta_y
            z_new -= delta_z
            delta_x /= 10.0
            delta_y /= 10.0
            delta_z /= 10.0 
        
        step_size /= 10.0
    
    return cum_tau


# def GET_RANDOM_ARRAY(ran_array, nrandom, ran_ctr, init_seed):
def GET_RANDOM_ARRAY(): 
    ran_number = np.random.rand()
    # ran_number = np.random.uniform(0, 1)

    # if ran_ctr < nrandom:
    #     ran_number = ran_array[ran_ctr]
    #     ran_ctr += 1
    # else:
    #     # np.random.seed(init_seed)
    #     ran_array = np.random.rand(nrandom)
    #     ran_ctr = 0
    #     ran_number = ran_array[ran_ctr]
    #     ran_ctr += 1
        # init_seed = np.random.randint(0, np.iinfo(np.uint32).max)
        # Uncomment the following line to print the seed and ran_number
        # print(init_seed, ran_number)
    
    return ran_number
